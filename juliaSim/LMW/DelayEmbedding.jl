###############################
# Delay embedding demo
# - pick (j1, j2)
# - integrate once
# - delay-embed a scalar observable y_comp
###############################

include("TimeSeriesLMW.jl")

using Plots
using Main.MySystem      # assumes rhs, integrate, rk live here

# Plotly is nice for 3D interactive viewing
plotly()

###############################
# User parameters
###############################

size        = 8
upper_bound = 1e6

constant_params = (-1.0, 1.0, 1.0, -1.0, -1.8, -1.8)

# Choose a parameter pair (j1, j2) you suspect is chaotic
# (e.g. from your LLE map or classification)
j1 = 4.0   # <-- change these
j2 = -2.0

# Observable index (component of y)
obs_comp = 1  # 1 = y₁

# Fraction of trajectory to discard as transient
transient_frac = 0.3

# Delay embedding parameters
embed_dim = 3       # m = 2 or 3 is typical for plotting
lag_steps = 10      # τ in "index units" (careful: not physical time, but samples)

###############################
# Helpers
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

# Build delay embedding of a scalar series x:
# Y[i, :] = [x[i], x[i+τ], x[i+2τ], ..., x[i + (m-1)τ]]
function delay_embedding(x::Vector{Float64}, m::Int, τ::Int)
    N = length(x) - (m - 1)*τ
    if N <= 0
        error("Not enough points for delay embedding: length(x) = $(length(x)), m = $m, τ = $τ")
    end
    Y = Array{Float64}(undef, N, m)
    @inbounds for i in 1:N
        for j in 0:(m-1)
            Y[i, j+1] = x[i + j*τ]
        end
    end
    return Y
end

###############################
# Integrate system for chosen (j1, j2)
###############################

params = (constant_params..., j1, j2)
f = rhs(params...)

series, flag = integrate(f, size; upper_bound = upper_bound)

if flag == 1
    error("Integration blew up for (j1, j2) = ($j1, $j2)")
end

n = length(series)
if n < 50
    error("Time series too short (n = $n) for delay embedding.")
end

# Discard transient
start_idx = max(1, Int(floor(transient_frac * n)))
series_tail = @view series[start_idx:end]

# Extract observable on tail
obs_tail = extract_observable(series_tail, size, obs_comp)

###############################
# Build delay embedding
###############################

Y = delay_embedding(obs_tail, embed_dim, lag_steps)

###############################
# Plots
###############################

# 2D embedding: (x(t), x(t+τ))
save_dir = @__DIR__
if embed_dim == 2
    de2 = scatter(
        Y[:,1],
        Y[:,2];
        markersize = 1.5,
        markerstrokewidth = 0.0,
        xlabel = "x(t)",
        ylabel = "x(t + τ)",
        title  = "2D delay embedding (j₁ = $j1, j₂ = $j2, τ = $lag_steps, y_comp = $obs_comp)",
        legend = false,
        size   = (700, 600),
    )
    display(de2)
    savefig(de2, joinpath(save_dir, "Plots/delay_embedding_2D_j1_$(j1)_j2_$(j2).png"))
end

# 3D embedding: (x(t), x(t+τ), x(t+2τ))
if embed_dim == 3
    de3 = scatter(
        Y[:,1],
        Y[:,2],
        Y[:,3];
        markersize = 1.5,
        markerstrokewidth = 0.0,
        xlabel = "x(t)",
        ylabel = "x(t + τ)",
        zlabel = "x(t + 2τ)",
        title  = "3D delay embedding (j₁ = $j1, j₂ = $j2, τ = $lag_steps, y_comp = $obs_comp)",
        legend = false,
        size   = (700, 600),
    )
    display(de3)
    savefig(de3, joinpath(save_dir, "Plots/delay_embedding_3D_j1_$(j1)_j2_$(j2).png"))
end
