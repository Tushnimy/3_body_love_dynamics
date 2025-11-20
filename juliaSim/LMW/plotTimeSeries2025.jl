###############################
# Dependencies & system include
###############################

include("TimeSeriesLMW.jl")

using Plots
using Main.MySystem      # assumes rhs, integrate, rk live here
using LinearAlgebra: norm

plotly()

###############################
# Parameters / thresholds
###############################

upper_bound        = 1e6
fp_threshold       = 1e-2   # for state variation and derivative norm
periodic_threshold = 1e-2   # spatial closeness for maxima
time_threshold     = 1e-2   # temporal closeness for periods

# Classification codes:
# 0 = unclassified / other
# 1 = periodic
# 2 = fixed point
# 3 = blow-up / integration failure
# 4 = chaotic / quasi-periodic (heuristic)

###############################
# Low-allocation helpers
###############################

"""
    state_variation(time_series, size)

Max variation of the state vector y over all rows in `time_series`.

Rows are assumed to be [t, y1..ysize, f1..fsize].
"""
function state_variation(time_series, size::Int)
    n = length(time_series)
    n == 0 && return 0.0

    first = time_series[1]
    ymin = Vector{Float64}(undef, size)
    ymax = Vector{Float64}(undef, size)

    @inbounds for j in 1:size
        v = first[1 + j]
        ymin[j] = v
        ymax[j] = v
    end

    @inbounds for i in 2:n
        row = time_series[i]
        for j in 1:size
            v = row[1 + j]
            v < ymin[j] && (ymin[j] = v)
            v > ymax[j] && (ymax[j] = v)
        end
    end

    # Euclidean norm of (max - min)
    s = 0.0
    @inbounds for j in 1:size
        d = ymax[j] - ymin[j]
        s += d*d
    end
    return sqrt(s)
end

"""
    derivative_norm(row, size)

Euclidean norm of f(y) given one row [t, y..., f...].
"""
function derivative_norm(row, size::Int)
    s = 0.0
    first_f = size + 2
    last_f  = 1 + 2*size
    @inbounds for k in first_f:last_f
        v = row[k]
        s += v*v
    end
    return sqrt(s)
end

###############################
# Fixed-point detection
###############################

"""
    is_fixed_point(tail, size; eps_state, eps_deriv)

Check if the trajectory `tail` has small variation in y and small |f(y)|
at the final time.
"""
function is_fixed_point(
    tail,
    size::Int;
    eps_state::Float64,
    eps_deriv::Float64
)
    n = length(tail)
    n < 5 && return false  # not enough data

    # 1) State variation small
    if state_variation(tail, size) > eps_state
        return false
    end

    # 2) Derivative at last point small
    lastrow = tail[end]
    derivative_norm(lastrow, size) < eps_deriv
end

###############################
# Maxima detection (low allocation)
###############################

"""
    find_maxima(time_series, f, size)

Find approximate local maxima of y₁ based on sign change of f₁.
Returns Vector{Vector{Float64}} with each element [t_max, y...].
"""
function find_maxima(time_series, f, size::Int)
    maxima_data = Vector{Vector{Float64}}()
    n = length(time_series)
    n < 2 && return maxima_data

    @inbounds for i in 1:(n - 1)
        # derivative component f1 at rows i and i+1
        d1 = time_series[i][size + 2]
        d2 = time_series[i+1][size + 2]

        if d1 > 0.0 && d2 <= 0.0
            t1 = time_series[i][1]
            t2 = time_series[i+1][1]

            y  = @view time_series[i][2:size+1]       # state at i
            fy = @view time_series[i][size+2:end]     # derivative at i

            # First interpolation estimate for zero of f1
            te = (t2 - t1) * (0.0 - d1) / (d2 - d1)
            ye, fye, _ = rk(y, f, fy, te)

            # Refine using your original logic
            if fye[1] > 0.0
                te = te + (te - t2) * (0.0 - fye[1]) / (fye[1] - d2)
                ye, fye = rk(y, f, fy, te)
            else
                te = te * (0.0 - d1) / (fye[1] - d1)
                ye, fye, _ = rk(y, f, fy, te)
            end

            # Store [te, ye...]
            vec = Vector{Float64}(undef, 1 + size)
            vec[1] = te
            @inbounds for j in 1:size
                vec[1 + j] = ye[j]
            end
            push!(maxima_data, vec)
        end
    end

    return maxima_data
end

###############################
# Periodicity detection (no FFT)
###############################

"""
    is_periodic_from_maxima(maxima_data, size;
                            space_eps, time_eps, min_cycles)

Use recurrence of maxima in state space and nearly constant period
to detect periodic orbits.
"""
function is_periodic_from_maxima(
    maxima_data::Vector{Vector{Float64}},
    size::Int;
    space_eps::Float64,
    time_eps::Float64,
    min_cycles::Int = 3
)
    n = length(maxima_data)
    n < min_cycles && return false

    # Reference maximum
    yref = @view maxima_data[1][2:1+size]
    tref = maxima_data[1][1]

    # Collect times when maxima are close to yref
    t_hits = Float64[tref]

    @inbounds for k in 2:n
        yk = @view maxima_data[k][2:1+size]

        # Compute norm(yk - yref) without allocation
        s = 0.0
        for j in 1:size
            d = yk[j] - yref[j]
            s += d*d
        end
        if sqrt(s) < space_eps
            push!(t_hits, maxima_data[k][1])
        end
    end

    length(t_hits) < min_cycles && return false

    # Periods between successive hits
    m = length(t_hits) - 1
    periods = Vector{Float64}(undef, m)
    @inbounds for k in 1:m
        periods[k] = t_hits[k+1] - t_hits[k]
    end

    meanP = sum(periods) / m
    maxdev = 0.0
    @inbounds for k in 1:m
        d = abs(periods[k] - meanP)
        d > maxdev && (maxdev = d)
    end

    return maxdev < time_eps
end

###############################
# Tail classification
###############################

"""
    classify_tail(tail, f, size;
                  fp_eps_state, fp_eps_deriv,
                  space_eps, time_eps)

Return:
 2 = fixed point
 1 = periodic
 4 = chaotic / quasi-periodic (heuristic)
 0 = none of the above
"""
function classify_tail(
    tail,
    f,
    size::Int;
    fp_eps_state::Float64,
    fp_eps_deriv::Float64,
    space_eps::Float64,
    time_eps::Float64
)::UInt8
    n = length(tail)
    n < 10 && return UInt8(0) # too short to say much

    # 1. Fixed point?
    if is_fixed_point(tail, size; eps_state = fp_eps_state, eps_deriv = fp_eps_deriv)
        return UInt8(2)
    end

    # 2. Periodic?
    maxima = find_maxima(tail, f, size)
    if is_periodic_from_maxima(maxima, size; space_eps = space_eps, time_eps = time_eps, min_cycles = 3)
        return UInt8(1)
    end

    # 3. Chaotic / quasi-periodic? (heuristic)
    #    Non-trivial variation + some maxima but not periodic
    if state_variation(tail, size) > fp_eps_state && length(maxima) >= 3
        return UInt8(4)
    end

    # 4. Otherwise, unclassified
    return UInt8(0)
end

###############################
# Parameter sweep & matrix M
###############################

nd = 150
M  = fill(UInt8(0), nd, nd)  # smaller element type

X = range(-3, stop = 3, length = nd)
Y = range(-3, stop = 3, length = nd)

# System parameters: last two entries (x, y) will be varied
constant_params = (-1.0, 1.0, 1.0, -1.0, -1.8, -1.8)
size = 8  # dimension of state vector

for (i, x) in enumerate(X)
    for (j, y) in enumerate(Y)
        params = (constant_params..., x, y)
        f = rhs(params...)  # build vector field once for this (x, y)

        series, flag = integrate(f, size; upper_bound = upper_bound)

        # 3 = blow-up / integration failure
        if flag == 1
            M[j, i] = UInt8(3)
            continue
        end

        n = length(series)
        n < 10 && continue  # leave as 0 (unclassified) if too short

        tail_start = Int(floor(0.7 * n))
        tail = @view series[tail_start:end]

        cls = classify_tail(
            tail, f, size;
            fp_eps_state = fp_threshold,
            fp_eps_deriv = fp_threshold,
            space_eps    = periodic_threshold,
            time_eps     = time_threshold
        )

        M[j, i] = cls
    end
end

###############################
# Basic visualization (we'll refine later)
###############################

parameterPlot = heatmap(
    X,
    Y,
    M;
    xlabel = "j1",
    ylabel = "j2",
    size   = (600, 600)
)

display(parameterPlot)
