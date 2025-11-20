##########################
# Setup & imports
##########################

include("LyapTryFinal.jl")
using LaTeXStrings
using Main.lyapLMW
using LinearAlgebra
using Plots
using DelimitedFiles

# use default backend (GR) for speed; no plotly() to avoid warnings

##########################
# Parameters
##########################

nd   = 100          # grid resolution
size = 8            # state dimension

X = range(-3, stop = 3, length = nd)
Y = range(-3, stop = 3, length = nd)

constant_params3 = (-1.0, 1.0, 1.0, -1.0, -1.8, -1.8)
constant_params  = (-1.0, 1.0)   # currently unused here

# LLE_rk4 parameters
T_final    = 1e4   # total integration time
transient  = 10.0  # time to discard before averaging
dt         = 1e-3
nsteps_max = 10000

periodicTol = 1e-3  # |LLE| < periodicTol treated as "periodic-ish"

# Matrix of Lyapunov exponents: NaN = unbound / undefined
M = fill(NaN, nd, nd)

# Periodic candidates in parameter space (for later plotting if you like)
px = Float64[]
py = Float64[]

##########################
# LLE sweep
##########################

# Option 1: fixed initial tangent direction for all parameters
u0 = randn(size)
u0 ./= norm(u0)

for (i, x) in enumerate(X)
    for (j, y) in enumerate(Y)
        # copy initial direction so we don't mutate u0
        u = copy(u0)

        params3 = (constant_params3..., x, y)

        LLE = LLE_rk4(rhs3(params3), u, T_final, transient, dt, size, nsteps_max)

        if LLE == "unbound"
            # leave M[j,i] = NaN to mark unbounded/escape
            continue
        end

        # Store LLE
        M[j, i] = LLE

        # Mark "periodic-ish" based on |LLE| ~ 0
        if abs(LLE) < periodicTol
            push!(px, x)
            push!(py, y)
        end
    end

    # coarse progress output
    @info "Finished column $i / $nd"
end

##########################
# Quick Julia plot (optional)
##########################

finite_vals = filter(!isnan, vec(M))
minL = isempty(finite_vals) ? -0.1 : minimum(finite_vals)
maxL = isempty(finite_vals) ?  0.1 : maximum(finite_vals)

# symmetric color scale around 0
maxAbs = max(abs(minL), abs(maxL))
clims = (-maxAbs, maxAbs)

plt = heatmap(
    X, Y, M;
    xlabel  = L"j_1",
    ylabel  = L"j_2",
    title   = "Largest Lyapunov exponent",
    size    = (800, 700),
    colorbar = true,
    clims   = clims,
)

plot!(plt, X, X; color = :black, linestyle = :dash, label = L"j_1=j_2")
scatter!(plt, px, py; color = :orange, ms = 3, label = "|LLE| < $(periodicTol)")

display(plt)

##########################
# Export matrices to same directory as this script
##########################

save_dir = @__DIR__

# LLE matrix
writedlm(joinpath(save_dir, "M_LLE.csv"), M, ',')

# parameter grids
writedlm(joinpath(save_dir, "X_LLE.csv"), collect(X), ',')
writedlm(joinpath(save_dir, "Y_LLE.csv"), collect(Y), ',')

# periodic candidates (optional)
writedlm(joinpath(save_dir, "periodic_px.csv"), px, ',')
writedlm(joinpath(save_dir, "periodic_py.csv"), py, ',')
