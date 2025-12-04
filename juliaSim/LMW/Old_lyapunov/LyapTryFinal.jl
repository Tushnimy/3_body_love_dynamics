module lyapLMW

export rhs, rk, rk4, integrate_steps, integrate, LLE_rk4

using LinearAlgebra
using StaticArrays

###########################################################
# Right-hand side: 8D system with 9 parameters
# p = (a1, a2, a3, a4, b1, b2, b3, j1, j2)
###########################################################

function rhs(p)
    @assert length(p) == 9 "rhs expects 9 parameters (a1,a2,a3,a4,b1,b2,b3,j1,j2)"

    a1, a2, a3, a4, b1, b2, b3, j1, j2 = p

    return function (u::SVector{8,Float64})
        @SVector [
            # oscillator 1  (L1, M)
            a1 + (u[3]^2 - u[4]^2) + b1*(u[1] + u[4]),
            2*u[3]*u[4]           + b1*(u[2] - u[3]),

            # oscillator 2 (Majnu, coupled via j1 to L2)
            a2 + (u[1]^2 - u[2]^2) + b2*(u[3] + u[2]) + j1*u[5],
            2*u[1]*u[2]           + b2*(u[4] - u[1]) + j1*u[6],

            # oscillator 3 (Ward, coupled via j2 to L1)
            a4 + (u[7]^2 - u[8]^2) + b1*(u[5] + u[8]),
            2*u[7]*u[8]           + b1*(u[6] - u[7]),

            a3 + (u[5]^2 - u[6]^2) + b3*(u[7] + u[6]) + j2*u[1],
            2*u[5]*u[6]           + b3*(u[8] - u[5]) + j2*u[2],
        ]
    end
end

###########################################################
# Dormand–Prince 5(4) embedded Runge–Kutta step
###########################################################

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

###########################################################
# Adaptive integrator using DOPRI (with rejection)
#
# Returns:
#   u    : state at t = tend (or earlier if stopped)
#   flag : 0 = OK, 1 = hit upper_bound/NaN, 2 = h < hmin,
#          3 = exceeded max_steps
#
# Suggested: tol ~ 1e-4–1e-6, upper_bound ~ 1e6–1e10
###########################################################

function integrate(f, tend, u0;
                   tol::Float64 = 1e-4,
                   upper_bound::Float64 = 1e10,
                   h0::Float64 = 1e-3,
                   hmin::Float64 = 1e-8,
                   hmax::Float64 = 1e-2)

    u = u0
    t = 0.0
    h = h0
    fu = f(u)
    flag = 0

    while t < tend
        if h < hmin
            flag = 2               # step too small
            break
        end

        h = min(h, tend - t)

        # trial step
        z, fz, est = rk(u, f, fu, h)

        # scale error by state size
        sc = max(norm(u), 1.0)
        err = est / sc

        if err <= tol || est == 0.0
            # accept step
            t += h
            u, fu = z, fz

            if !isfinite(norm(u)) || norm(u) > upper_bound
                flag = 1          # blew up / NaN
                break
            end

            # choose next h
            fac = (err == 0.0) ? 2.0 : (tol/err)^(1/5)
            h = h * min(2.0, max(0.3, 0.8*fac))
            h = min(h, hmax)
        else
            # reject step, shrink h and retry
            fac = (tol/err)^(1/5)
            h = h * max(0.1, 0.8*fac)
        end
    end
    return u, flag
end

###########################################################
# Fixed-step RK4
###########################################################

function rk4(y, f, fy, h)
    k1 = fy
    k2 = f(y + h * k1 / 2)
    k3 = f(y + h * k2 / 2)
    k4 = f(y + h * k3)
    ynew = y + (h/6) * (k1 + 2k2 + 2k3 + k4)
    return ynew, f(ynew)
end

###########################################################
# Integrate for a fixed number of RK4 steps
###########################################################

function integrate_steps(f, steps, h, u0)
    u = u0
    fu = f(u0)
    nsteps = Int(steps)
    for _ = 1:nsteps
        u, fu = rk4(u, f, fu, h)
    end
    return u, fu
end

###########################################################
# Largest Lyapunov exponent via Benettin method (RK4)
#
# Arguments:
#   f     : rhs(u) from rhs(p)
#   u0    : initial condition (SVector{8,Float64})
#   iter  : number of Benettin renormalisation loops
#   steps : RK4 steps per loop
#   h     : RK4 step size
#   size  : dimension (8 for this system)
#   Ttr   : transient time (integrated with adaptive DOPRI)
#
# Recommended ranges for this system:
#   h    ~ 5e-4 – 1e-3
#   steps ~ 100 – 200
#   iter  ~ 3000 – 8000
#   Ttr   ~ 200 – 600
###########################################################

function LLE_rk4(f, u0, iter, steps, h, size, Ttr;
                 upper_bound::Float64 = 1e8,
                 max_steps::Int = 2_000_000)

    # 1. Transient with adaptive integrator
    ur, flag = integrate(f, Ttr, u0;
                         tol         = 1e-5,
                         upper_bound = upper_bound,
                         h0          = h,
                         hmax        = 10h)

    # any nonzero flag → treat as unbounded/unreliable
    if flag != 0 || !isfinite(norm(ur))
        return "unbound"
    end

    # 2. Initial small perturbation
    dev_vec = SVector{size,Float64}(randn(size))
    ur_vec  = ur
    ut      = ur_vec + (dev_vec / norm(dev_vec)) * 1e-8

    log_factors = Float64[]

    # 3. Benettin loops (fixed-step RK4)
    for _ in 1:iter
        ur_vec, _ = integrate_steps(f, steps, h, ur_vec)
        ut,     _ = integrate_steps(f, steps, h, ut)

        # blow-up / NaN guard
        if !isfinite(norm(ur_vec)) || !isfinite(norm(ut)) ||
           norm(ur_vec) > upper_bound || norm(ut) > upper_bound
            return "unbound"
        end

        diff = ut - ur_vec
        sep  = norm(diff)
        if sep == 0.0
            continue
        end

        push!(log_factors, log(sep / 1e-8))

        # renormalise perturbation
        dir = diff / sep
        ut  = ur_vec + dir * 1e-8
    end

    if isempty(log_factors)
        return 0.0
    end

    # total physical time covered by the RK4 evolution
    T = iter * steps * h
    return sum(log_factors) / T
end

end # module
