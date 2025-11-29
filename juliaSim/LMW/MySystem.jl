module MySystem

export rhs, rk, integrate

using LinearAlgebra: norm
using StaticArrays

###############################
# Vector field
###############################

"""
    rhs(a1, a2, a3, a4, b1, b2, b3, j1, j2)

Return a closure `f(x)` implementing the 8D vector field.
"""
function rhs(a1::Float64, a2::Float64, a3::Float64, a4::Float64,
             b1::Float64, b2::Float64, b3::Float64,
             j1::Float64, j2::Float64)
    function f(x)
        @SVector [
            a1 + (x[3]^2 - x[4]^2) + b1*(x[1] + x[4]),
            2*x[3]*x[4]           + b1*(x[2] - x[3]),
            a2 + (x[1]^2 - x[2]^2) + b2*(x[3] + x[2]) + j1*x[5],
            2*x[1]*x[2]           + b2*(x[4] - x[1]) + j1*x[6],
            a4 + (x[7]^2 - x[8]^2) + b1*(x[5] + x[8]),
            2*x[7]*x[8]           + b1*(x[6] - x[7]),
            a3 + (x[5]^2 - x[6]^2) + b3*(x[7] + x[6]) + j2*x[1],
            2*x[5]*x[6]           + b3*(x[8] - x[5]) + j2*x[2]
        ]
    end
end

###############################
# Dormand–Prince RK(5,4) step
###############################

"""
    rk(y, f, fy, h) -> (z, fz, est)

One adaptive RK(5,4) step (Dormand–Prince) from state `y` with derivative
`fy = f(y)` and step size `h`.

Returns:
  z   – 5th order solution
  fz  – f(z)
  est – error estimate ‖z - w‖₂
"""
function rk(y, f, fy, h)
    k1 = fy
    k2 = f(y + h * (1//5) * k1)
    k3 = f(y + h * (3//40 * k1 + 9//40 * k2))
    k4 = f(y + h * (44//45 * k1 - 56//15 * k2 + 32//9 * k3))
    k5 = f(y + h * (19372//6561 * k1 - 25360//2187 * k2 +
                    64448//6561 * k3 - 212//729 * k4))
    k6 = f(y + h * (9017//3168 * k1 - 355//33 * k2 + 46732//5247 * k3 +
                    49//176 * k4 - 5103//18656 * k5))

    # 5th-order solution
    z = y + h * (35//384 * k1 + 500//1113 * k3 + 125//192 * k4 -
                 2187//6784 * k5 + 11//84 * k6)

    k7 = f(z)

    # 4th-order solution
    w = y + h * (5179//57600 * k1 + 7571//16695 * k3 +
                 393//640 * k4 - 92097//339200 * k5 +
                 187//2100 * k6 + 1//40 * k7)

    est = norm(z - w)
    return z, k7, est
end

###############################
# Adaptive integrator with
#   - nsteps cap
#   - downsampled storage
#   - rolling tail buffer
###############################

"""
    integrate(f, size; kwargs...) -> (time_series, flag)

Integrate `y' = f(y)` in ℝ^size with an adaptive RK(5,4) method.

Returns:
  time_series :: Vector{Vector{Float64}} with rows [t, y..., f...]
  flag        :: Int
                0 → success
                1 → step too small / blow-up / max steps reached

Keyword arguments (defaults chosen to be safe for large 2D scans):

  tol          = 1e-6      – error tolerance
  tend         = 1e4       – final time
  upper_bound  = 1e5       – blow-up radius for state norm
  h0           = 1e-3      – initial step
  h_min        = 1e-8      – minimum step
  h_max        = 1e-1      – maximum step
  nsteps_max   = 2_000_000 – hard cap on accepted steps
  save_every   = 10        – store every N-th accepted step
  max_store    = 5000      – keep at most this many stored rows (rolling tail)
"""
function integrate(
    f,
    size::Int;
    tol::Float64        = 1e-6,
    tend::Float64       = 1e4,
    upper_bound::Float64 = 1e5,
    h0::Float64         = 1e-3,
    h_min::Float64      = 1e-8,
    h_max::Float64      = 1e-1,
    nsteps_max::Int     = 2_000_000,
    save_every::Int     = 10,
    max_store::Int      = 5000,
)

    # Initial condition (you can change this if you want deterministic ICs)
    y = randn(size)

    t  = 0.0
    h  = h0
    fy = f(y)
    flag = 0
    nsteps = 0
    save_counter = 0

    # time_series rows: [t, y..., fy...]
    time_series = Vector{Vector{Float64}}()
    # store initial row
    row0 = Vector{Float64}(undef, 1 + 2*size)
    row0[1] = t
    @inbounds for j in 1:size
        row0[1 + j]      = y[j]
        row0[1 + size+j] = fy[j]
    end
    push!(time_series, row0)

    while t < tend
        # if step size gets too small, give up gracefully
        if h < h_min
            @info "Step size too small at t = $t"
            flag = 1
            break
        end

        if nsteps >= nsteps_max
            @info "Reached nsteps_max = $nsteps_max at t = $t, aborting."
            flag = 1
            break
        end

        # don't step past tend
        h = min(h, tend - t)

        # trial step
        z, fz, est = rk(y, f, fy, h)

        # scale error
        sc = max(norm(y), 1.0)
        err = est / sc

        if err <= tol || est == 0.0
            # ---- ACCEPT step ----
            t += h
            y, fy = z, fz
            nsteps += 1
            save_counter += 1

            # blow-up / non-finite detection
            ny = norm(y)
            if !isfinite(ny) || ny > upper_bound
                @info "Blow-up / escape detected: ‖y‖ = $ny at t = $t"
                flag = 1
                break
            end

            # store every `save_every`-th accepted step
            if save_counter % save_every == 0
                row = Vector{Float64}(undef, 1 + 2*size)
                row[1] = t
                @inbounds for j in 1:size
                    row[1 + j]      = y[j]
                    row[1 + size+j] = fy[j]
                end

                if length(time_series) < max_store
                    push!(time_series, row)
                else
                    # rolling tail buffer: drop oldest, keep latest max_store points
                    popfirst!(time_series)
                    push!(time_series, row)
                end
            end

            # choose next h (growth/shrink factors capped)
            fac = (err == 0.0) ? 2.0 : (tol / err)^(1/5)
            h   = h * min(2.0, max(0.3, 0.8 * fac))
            h   = min(h, h_max)

        else
            # ---- REJECT step: shrink h and retry ----
            fac = (tol / err)^(1/5)
            h   = h * max(0.1, 0.8 * fac)
        end
    end

    return time_series, flag
end

end # module
