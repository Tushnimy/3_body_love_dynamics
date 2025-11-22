###############################
# Dependencies & system include
###############################

include("TimeSeriesLMW.jl")

using Plots
using Main.MySystem      # assumes rhs, integrate, rk live here
using LinearAlgebra: norm
using FFTW: fft

# You can comment this out if Plotly isn't installed:
plotly()

###############################
# Parameters / thresholds
###############################

upper_bound        = 1e6

fp_threshold       = 1e-2   # for state variation and derivative norm
periodic_threshold = 1e-2   # spatial closeness for maxima (state space)
time_threshold     = 1e-2   # temporal closeness for periods

# FFT-based thresholds (heuristic, you can tune)
fft_min_points       = 256       # minimal number of samples for FFT
fft_max_points       = 4096      # maximum number of samples (we downsample if longer)
fft_periodic_ratio   = 0.7       # dominant peak power / total power to call it periodic
fft_quasi_low_ratio  = 0.3       # lower bound for quasi-periodic

# Classification codes:
# 0 = unclassified / other
# 1 = periodic
# 2 = fixed point
# 3 = blow-up / integration failure
# 4 = quasi-periodic
# 5 = chaotic / irregular

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
# Periodicity from maxima (no FFT)
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
# FFT-based spectral analysis (no sort, downsample)
###############################

"""
    fft_metrics_component(tail, size, comp;
                          min_n, max_n)

Compute FFT of y_comp over the tail (downsampled if necessary) and return:

    (valid_fft, dominant_ratio, n_eff, total_power)

- dominant_ratio = (largest peak power excluding DC) / total_power
- n_eff = effective number of active frequencies (participation ratio)
          ≈ 1 for single-line spectrum, larger for spread-out spectra.
"""
function fft_metrics_component(
    tail,
    size::Int,
    comp::Int;
    min_n::Int = fft_min_points,
    max_n::Int = fft_max_points
)
    n = length(tail)
    if n < min_n
        return (false, 0.0, 0.0, 0.0)
    end

    # Decide sampling step if downsampling
    if n <= max_n
        step = 1
        m = n
    else
        step = max(1, fld(n, max_n))
        m = fld(n, step)
    end

    m < min_n && return (false, 0.0, 0.0, 0.0)

    # Extract component `comp` at stride `step`
    vals = Vector{Float64}(undef, m)
    idx = 1
    @inbounds for k in 1:m
        vals[k] = tail[idx][1 + comp]
        idx += step
    end

    # Remove mean (DC)
    meanv = sum(vals) / m
    @inbounds for k in 1:m
        vals[k] -= meanv
    end

    Y = fft(vals)

    # Use only positive frequencies (skip DC at index 1)
    kmax = fld(m, 2)
    if kmax < 2
        return (false, 0.0, 0.0, 0.0)
    end

    total_power = 0.0
    dominant_power = 0.0
    sum_p2 = 0.0

    @inbounds for k in 2:(kmax + 1)
        p = abs2(Y[k])
        total_power += p
        sum_p2 += p * p
        p > dominant_power && (dominant_power = p)
    end

    total_power <= 0 && return (false, 0.0, 0.0, 0.0)

    dominant_ratio = dominant_power / total_power
    # "effective" number of contributing frequencies
    n_eff = total_power^2 / sum_p2

    return (true, dominant_ratio, n_eff, total_power)
end

###############################
# Tail classification
###############################

"""
    classify_tail(tail, f, size;
                  fp_eps_state, fp_eps_deriv,
                  space_eps, time_eps,
                  fft_periodic_ratio, fft_quasi_low_ratio)

Return codes:
 2 = fixed point
 1 = periodic
 4 = quasi-periodic
 5 = chaotic / irregular
 0 = none of the above / unclassified
"""
function classify_tail(
    tail,
    f,
    size::Int;
    fp_eps_state::Float64,
    fp_eps_deriv::Float64,
    space_eps::Float64,
    time_eps::Float64,
    fft_periodic_ratio::Float64,
    fft_quasi_low_ratio::Float64,
    fft_comp::Int = 1,          # which y-component to FFT (1 = y1)
    fft_min_points::Int = fft_min_points,
    fft_max_points::Int = fft_max_points
)::UInt8
    n = length(tail)
    n < 10 && return UInt8(0) # too short to say much

    # 1. Fixed point?
    if is_fixed_point(tail, size; eps_state = fp_eps_state, eps_deriv = fp_eps_deriv)
        return UInt8(2)
    end

    # 2. Maxima-based info
    maxima = find_maxima(tail, f, size)
    periodic_by_maxima = is_periodic_from_maxima(
        maxima, size;
        space_eps = space_eps,
        time_eps  = time_eps,
        min_cycles = 3
    )

    # 3. FFT-based spectral info on chosen component
    valid_fft, dom_ratio, n_eff, total_power = fft_metrics_component(
        tail, size, fft_comp;
        min_n = fft_min_points,
        max_n = fft_max_points
    )

    # 4. Periodic classification (maxima AND/OR FFT)
    if periodic_by_maxima ||
       (valid_fft && dom_ratio >= fft_periodic_ratio && n_eff <= 2.0)
        return UInt8(1)  # periodic
    end

    # 5. Non-trivial dynamics?
    var_state = state_variation(tail, size)
    if var_state <= fp_eps_state
        # State almost flat but failed fixed-point test: ambiguous -> 0
        return UInt8(0)
    end

    # 6. Quasi-periodic vs chaotic / irregular
    # Heuristic:
    #  - valid FFT, moderate concentration of power over a handful of frequencies
    #    -> quasi-periodic (4)
    #  - else -> chaotic / irregular (5)
    if valid_fft
        if dom_ratio >= fft_quasi_low_ratio &&
           dom_ratio < fft_periodic_ratio &&
           2.0 <= n_eff <= 20.0
            return UInt8(4)  # quasi-periodic: a few active lines
        else
            return UInt8(5)  # chaotic / irregular spectrum
        end
    else
        # no reliable FFT but nontrivial variation: treat as chaotic/irregular
        return UInt8(5)
    end
end

###############################
# Parameter sweep & matrix M
###############################

nd = 300
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
            fp_eps_state        = fp_threshold,
            fp_eps_deriv        = fp_threshold,
            space_eps           = periodic_threshold,
            time_eps            = time_threshold,
            fft_periodic_ratio  = fft_periodic_ratio,
            fft_quasi_low_ratio = fft_quasi_low_ratio,
            fft_comp            = 1,               # use y1; change to 2,3,... if needed
            fft_min_points      = fft_min_points,
            fft_max_points      = fft_max_points
        )

        M[j, i] = cls
    end
end


using DelimitedFiles

# Directory of this script
save_dir = @__DIR__

# Convert M to plain Int matrix
M_int = Array{Int}(M)

# Save all files in the script directory
writedlm(joinpath(save_dir, "Data/M.csv"), M_int, ',')
writedlm(joinpath(save_dir, "Data/X.csv"), collect(X), ',')
writedlm(joinpath(save_dir, "Data/Y.csv"), collect(Y), ',')


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


