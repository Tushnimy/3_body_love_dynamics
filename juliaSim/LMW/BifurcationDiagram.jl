###############################
# Bifurcation-style diagram
# - vary j1, fix j2
# - plot local maxima of y₁ vs j1
# - skip blow-up trajectories
###############################

include("TimeSeriesLMW.jl")

using Plots
using Main.MySystem      # assumes rhs, integrate, rk live here
using LinearAlgebra: norm
using DelimitedFiles     # <-- NEW: for writing CSV/TSV data

# Use GR backend (no Plotly dependency)
plotly()

###############################
# User parameters
###############################

size         = 8           # state dimension
upper_bound  = 1e6         # same as your other scripts

# System parameters: last two are (j1, j2)
constant_params = (-1.0, 1.0, 1.0, -1.0, -1.8, -1.8)

# Bifurcation sweep: vary j1, fix j2
j1_min, j1_max = 2, 6
nj1            = 5000
j2_fixed       = -2.0       # keep this Float64

j1_values = range(j1_min, stop=j1_max, length=nj1)

# Observable to use for maxima (component index in y)
obs_comp = 1               # 1 = y₁, 2 = y₂, ...

# Transient fraction to discard
transient_frac = 0.5       # ignore first 50% of samples

###############################
# Blow-up detection (copied & slightly trimmed)
###############################

"""
    is_blowup_trajectory(series, size;
                         R_escape=50.0,
                         growth_factor=5.0,
                         tail_frac=0.15,
                         min_points_tail=20)

Heuristic blow-up detector:

- Compute r_k = ‖y_k‖₂ along the trajectory.
- If r_max < R_escape => not blow-up.
- Find first index where r_k ≥ R_escape.
- Look at the tail from max(that index, floor((1-tail_frac)*n)) to end.
- If the tail norms are mostly increasing and the final norm is at least
  `growth_factor * R_escape`, declare blow-up.
"""
function is_blowup_trajectory(
    series,
    size::Int;
    R_escape::Float64 = 50.0,
    growth_factor::Float64 = 5.0,
    tail_frac::Float64 = 0.15,
    min_points_tail::Int = 20,
)
    n = length(series)
    n == 0 && return false

    # Norm of y at each stored time
    norms = Vector{Float64}(undef, n)
    @inbounds for k in 1:n
        row = series[k]
        s = 0.0
        for j in 1:size
            v = row[1 + j]
            s += v*v
        end
        norms[k] = sqrt(s)
    end

    r_max = maximum(norms)
    # Never even reached the escape radius -> not blow-up
    r_max < R_escape && return false

    # First time we cross the escape radius
    idx_first = findfirst(>=(R_escape), norms)
    idx_first === nothing && return false

    # Tail region to inspect for monotone-ish escape
    start_tail = max(idx_first, Int(floor((1 - tail_frac) * n)))
    start_tail = min(start_tail, n)  # safety
    length(start_tail:n) < min_points_tail && return false

    tail_norms = @view norms[start_tail:end]

    # Simple monotonicity check: fraction of differences ≥ 0
    inc_count = 0
    total_diff = length(tail_norms) - 1
    @inbounds for k in 1:total_diff
        if tail_norms[k+1] >= tail_norms[k]
            inc_count += 1
        end
    end
    frac_increasing = total_diff > 0 ? inc_count / total_diff : 0.0

    # Require:
    # - norms mostly increasing in the tail
    # - final norm significantly beyond R_escape
    if frac_increasing ≥ 0.8 && tail_norms[end] ≥ growth_factor * R_escape
        return true
    else
        return false
    end
end

###############################
# Helpers for bifurcation
###############################

# Extract scalar observable y_comp from series of [t, y..., f...]
function extract_observable(series, size::Int, comp::Int)
    n = length(series)
    obs = Vector{Float64}(undef, n)
    idx = 1 + comp
    @inbounds for k in 1:n
        obs[k] = series[k][idx]
    end
    return obs
end

# Simple local maxima of a scalar series
function local_maxima_from_vector(x::Vector{Float64})
    n = length(x)
    maxima = Float64[]
    n < 3 && return maxima
    @inbounds for k in 2:(n-1)
        if x[k] >= x[k-1] && x[k] >= x[k+1]
            push!(maxima, x[k])
        end
    end
    return maxima
end

###############################
# Main loop: build bifurcation data
###############################

bx = Float64[]   # parameter values (j1) at maxima
by = Float64[]   # corresponding observable maxima

for (i, j1) in enumerate(j1_values)
    # Ensure Float64 for params
    params = (constant_params..., float(j1), float(j2_fixed))
    f = rhs(params...)  # vector field for this (j1, j2_fixed)

    series, flag = integrate(f, size; upper_bound = upper_bound)

    # 1) If your integrator already flags blow-up, skip directly
    if flag == 1
        @info "Integration flagged blow-up at j1 = $j1, j2 = $j2_fixed"
        continue
    end

    n = length(series)
    if n < 20
        @info "Too short series at j1 = $j1, skipping"
        continue
    end

    # 2) Extra blow-up sanity check
    if is_blowup_trajectory(series, size;
                            R_escape = 50.0,
                            growth_factor = 2.0,
                            tail_frac = 0.15,
                            min_points_tail = 20)
        @info "Heuristic blow-up detected at j1 = $j1, j2 = $j2_fixed"
        continue
    end

    # Discard transient
    start_idx = max(2, Int(floor(transient_frac * n)))
    tail_series = @view series[start_idx:end]

    # Extract observable on the tail
    obs_tail = extract_observable(tail_series, size, obs_comp)

    # Find local maxima (scalar series)
    maxvals = local_maxima_from_vector(obs_tail)

    if !isempty(maxvals)
        for mv in maxvals
            push!(bx, j1)
            push!(by, mv)
        end
    end

    if i % 20 == 0
        @info "Processed $(i) / $(nj1) j1-values"
    end
end

###############################
# Save bifurcation data for Python
###############################

save_dir = @__DIR__
plots_dir = joinpath(save_dir, "Plots")
Data_dir = joinpath(save_dir, "Data")

# CSV file with two columns: j1, max_y1
csv_path = joinpath(
    Data_dir,
    "bifurcation_data_j1_y1_j2_$(j2_fixed)_($(j1_min)-$(j1_max)).csv",
)

data = hcat(bx, by)  # n×2 matrix: [j1  max_y1]

open(csv_path, "w") do io
    # header
    writedlm(io, ["j1" "max_y1"], ',')
    # data
    writedlm(io, data, ',')
end

@info "Saved bifurcation data to $csv_path"

###############################
# Plot bifurcation diagram (optional)
###############################

fixedparambif_plot = scatter(
    bx,
    by;
    markersize = 0.5,
    markerstrokewidth = 0.0,
    xlabel = "j₁",
    ylabel = "maxima of y₁",
    title  = "Bifurcation-style diagram (j₂ = $j2_fixed, observable = y₁)",
    legend = false,
    size   = (800, 600),
    ylims  = (-2, 10),
)

display(fixedparambif_plot)

