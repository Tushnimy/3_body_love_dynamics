##########################
# Setup & imports
##########################

include("LyapJacLMW.jl")  

using LaTeXStrings
using .LyapJacLMW
using LinearAlgebra
using Plots
using DelimitedFiles

gr()  

##########################
# Parameters
##########################

size = 8
nd   = 300                

# Time parameters for Lyapunov computation
Ttr        = 3000.0      # transient time
T_LLE      = 800.0       # averaging time
Δt_renorm  = 0.05        # time between QR renormalizations

# Solver tolerances and bounds
reltol      = 1e-6
abstol      = 1e-8
upper_bound = 1e9

# base system parameters
constant_params = (-1.0, 1.0, 1.0, -1.0, 0.0, -1.0, -1.0)

# parameter grid in (j1, j2)
X = range(-3.0, stop = 3.0, length = nd)   # j1
Y = range(-3.0, stop = 3.0, length = nd)   # j2

# classification thresholds for data files
tol         = 5e-2
periodicTol = tol
chaosTol    = tol

# storage
M  = fill(NaN, nd, nd)   # largest Lyapunov exponent
px = Float64[]; py = Float64[]
cx = Float64[]; cy = Float64[]

# Diagnostics
# 0 = OK, 1 = lyapunov_spectrum returned nothing
diag_flag = fill(0, nd, nd)

##########################
#Initial condition 
##########################

u0 = randn(size)
u0 ./= norm(u0)

##########################
# LLE sweep
##########################

for (i, x) in enumerate(X)      # j1
    for (j, y) in enumerate(Y)  # j2

        params = (constant_params..., x, y)

        λ = lyapunov_spectrum(params;
                              u0          = u0,
                              Ttr         = Ttr,
                              T_LLE       = T_LLE,
                              Δt_renorm   = Δt_renorm,
                              pexp        = 1,
                              reltol      = reltol,
                              abstol      = abstol,
                              upper_bound = upper_bound)

        if λ === nothing
            diag_flag[j, i] = 1   # numerical failure / unbounded
            continue
        end

        LLE = λ[1]
        M[j, i] = LLE
        diag_flag[j, i] = 0

        # optional classification
        if abs(LLE) < periodicTol
            push!(px, x); push!(py, y)
        elseif LLE > chaosTol
            push!(cx, x); push!(cy, y)
        end

        if j % 20 == 0
            @info "Finished row $j / $nd for column $i"
        end
    end
    @info "Finished column $i / $nd"
end

##########################
# Plot
##########################

finite_vals = filter(!isnan, vec(M))
minL = isempty(finite_vals) ? -0.1 : minimum(finite_vals)
maxL = isempty(finite_vals) ?  0.1 : maximum(finite_vals)

maxAbs = max(abs(minL), abs(maxL))
clims = (-maxAbs, maxAbs)

plt = heatmap(
    X, Y, M;
    xlabel   = L"j_1",
    ylabel   = L"j_2",
    title    = "Largest Lyapunov exponent (Jacobian/QR method)",
    size     = (800, 700),
    colorbar = true,
    clims    = clims,
)

display(plt)

##########################
# Export
##########################

save_dir = @__DIR__
isdir(joinpath(save_dir, "Data")) || mkpath(joinpath(save_dir, "Data"))

writedlm(joinpath(save_dir, "Data/M_LLE_jacobian.csv"), M, ',')
writedlm(joinpath(save_dir, "Data/X_LLE_jacobian.csv"), collect(X), ',')
writedlm(joinpath(save_dir, "Data/Y_LLE_jacobian.csv"), collect(Y), ',')

# diagnostics
writedlm(joinpath(save_dir, "Data/diag_flag_jacobian.csv"), diag_flag, ',')

# classification (optional)
writedlm(joinpath(save_dir, "Data/periodic_px_jacobian.csv"), px, ',')
writedlm(joinpath(save_dir, "Data/periodic_py_jacobian.csv"), py, ',')
writedlm(joinpath(save_dir, "Data/chaotic_cx_jacobian.csv"),  cx, ',')
writedlm(joinpath(save_dir, "Data/chaotic_cy_jacobian.csv"),  cy, ',')
