##########################
# Setup & imports
##########################

include("LyapTryFinal.jl")  # defines module lyapLMW above
using LaTeXStrings
using Main.lyapLMW
using LinearAlgebra
using Plots
using DelimitedFiles
using StaticArrays
using Random

# use default backend (GR) for speed
gr()
Random.seed!(1234)   # for reproducibility (optional)

##########################
# Parameters
##########################

size = 8                   # state dimension
nd   = 80                  # 80x80 grid

Ttr        = 3000.0        # total transient time budget
T_LLE      = 500.0         # Lyapunov averaging time
dt         = 1e-4
steps      = 15
time_per_cycle = steps * dt
iter       = max(Int(floor(T_LLE / time_per_cycle)), 1)

@info "LLE parameters: iter=$iter, steps=$steps, dt=$dt, Ttr=$Ttr, T_LLE=$(iter*steps*dt)"

# cheap boundedness pre-check time per IC
T_check = 1000.0           # you can bump this if needed

# base system parameters (unchanged dynamics)
constant_params = (-1.0, 1.0, 1.0, -1.0, 0.0, -1.0, -1.0)
X = range(-3.0, stop = 3.0, length = nd)   # j1
Y = range(-3.0, stop = 3.0, length = nd)   # j2

tol         = 5e-2
periodicTol = tol
chaosTol    = tol

M  = fill(NaN, nd, nd)
px = Float64[]; py = Float64[]
cx = Float64[]; cy = Float64[]

upper_bound   = 1e9
max_ic_tries  = 8      # max random ICs per parameter

# diagnostics: 0 = OK, 1 = all ICs fail pre-check,
#              2 = some pass pre-check but all fail during LLE
diag_flag = fill(0, nd, nd)

##########################
# LLE sweep with IC pre-check & diagnostics
##########################

for (i, x) in enumerate(X)      # j1
    for (j, y) in enumerate(Y)  # j2

        params = (constant_params..., x, y)
        f = rhs(params)

        LLE_val = "unbound"
        cell_diag = 0
        any_passed_precheck = false

        # --- try several initial conditions ---
        for attempt in 1:max_ic_tries
            # random initial condition
            u0_vec = randn(size)
            u0_vec ./= norm(u0_vec)
            u0 = SVector{8,Float64}(u0_vec)

            # 1) cheap boundedness pre-check up to T_check
            u_check, flag_check = integrate(
                f, T_check, u0;
                tol         = 1e-6,
                upper_bound = upper_bound,
                h0          = dt,
                hmax        = 10dt,
            )

            if flag_check != 0 || !isfinite(norm(u_check))
                # this IC clearly escapes quickly -> try another
                cell_diag = 1          # so far, only pre-check failures
                continue
            end

            any_passed_precheck = true

            # 2) run LLE from the pre-checked state
            # remaining transient time after pre-check:
            #Ttr_eff = max(Ttr - T_check, 0.0)
            Ttr_eff    = 0.0

            LLE_val = LLE_rk4(
                f, u_check, iter, steps, dt, size, Ttr_eff;
                upper_bound = upper_bound,
            )

            if LLE_val != "unbound"
                # got a bounded orbit with a finite LLE → keep it
                cell_diag = 0
                break
            else
                # passed pre-check but blew up later (transient/LLE)
                cell_diag = 2
            end
        end

        if LLE_val == "unbound"
            # every IC we tried escaped / failed at some stage
            diag_flag[j, i] = cell_diag
            continue
        end

        # success: store LLE
        LLE = LLE_val::Float64
        M[j, i] = LLE
        diag_flag[j, i] = 0

        # classification just for overlays if you want
        if abs(LLE) < periodicTol
            push!(px, x); push!(py, y)
        elseif LLE > chaosTol
            push!(cx, x); push!(cy, y)
        end

        if j % 40 == 0
            @info "Finished row $j / $nd for column $i"
        end
    end
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
    xlabel   = L"j_1",
    ylabel   = L"j_2",
    title    = "Largest Lyapunov exponent",
    size     = (800, 700),
    colorbar = true,
    clims    = clims,
)

display(plt)

##########################
# Export matrices to same directory as this script
##########################

save_dir = @__DIR__
isdir(joinpath(save_dir, "Data")) || mkpath(joinpath(save_dir, "Data"))

# LLE matrix and parameter grids
writedlm(joinpath(save_dir, "Data/M_LLE.csv"), M, ',')
writedlm(joinpath(save_dir, "Data/X_LLE.csv"), collect(X), ',')
writedlm(joinpath(save_dir, "Data/Y_LLE.csv"), collect(Y), ',')

# Diagnostics: why cells are grey
writedlm(joinpath(save_dir, "Data/diag_flag.csv"), diag_flag, ',')

# periodic / chaotic candidate lists (optional)
writedlm(joinpath(save_dir, "Data/periodic_px.csv"), px, ',')
writedlm(joinpath(save_dir, "Data/periodic_py.csv"), py, ',')
writedlm(joinpath(save_dir, "Data/chaotic_cx_.csv"),  cx, ',')
writedlm(joinpath(save_dir, "Data/chaotic_cy.csv"),  cy, ',')
