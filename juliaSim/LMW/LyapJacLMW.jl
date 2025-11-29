module LyapJacLMW

using LinearAlgebra
using OrdinaryDiffEq
using OrdinaryDiffEq: Rodas4
const LYAP_ALG = Rodas4()


const N = 8  # state dimension

############################
# 1. System RHS (in-place)
############################

"""
    rhs!(du, u, p, t)

Right-hand side of the 8D system.

Arguments:
- du :: AbstractVector  (output)
- u  :: AbstractVector  (state of length 8)
- p  :: NTuple{9,Float64} = (a1,a2,a3,a4,b1,b2,b3,j1,j2)
- t  :: Float64 (ignored; autonomous system)
"""
function rhs!(du, u, p, t)
    @inbounds begin
        a1, a2, a3, a4, b1, b2, b3, j1, j2 = p
        x1, x2, x3, x4, x5, x6, x7, x8 = u

        # oscillator 1
        du[1] = a1 + (x3^2 - x4^2) + b1*(x1 + x4)
        du[2] = 2*x3*x4 + b1*(x2 - x3)

        # oscillator 2
        du[3] = a2 + (x1^2 - x2^2) + b2*(x3 + x2) + j1*x5
        du[4] = 2*x1*x2 + b2*(x4 - x1) + j1*x6

        # oscillator 3
        du[5] = a4 + (x7^2 - x8^2) + b1*(x5 + x8)
        du[6] = 2*x7*x8 + b1*(x6 - x7)

        du[7] = a3 + (x5^2 - x6^2) + b3*(x7 + x6) + j2*x1
        du[8] = 2*x5*x6 + b3*(x8 - x5) + j2*x2
    end
    return nothing
end


function jacobian!(J::AbstractMatrix, u::AbstractVector, p)
    @assert length(u) == 8 "jacobian! expects u of length 8"
    @assert size(J) == (8, 8) "jacobian! expects an 8×8 J"

    a1, a2, a3, a4, b1, b2, b3, j1, j2 = p
    x1, x2, x3, x4, x5, x6, x7, x8 = u

    # start from zero to be safe
    @inbounds begin
        J[1,1] = 0.0; J[1,2] = 0.0; J[1,3] = 0.0; J[1,4] = 0.0;
        J[1,5] = 0.0; J[1,6] = 0.0; J[1,7] = 0.0; J[1,8] = 0.0;
        J[2,1] = 0.0; J[2,2] = 0.0; J[2,3] = 0.0; J[2,4] = 0.0;
        J[2,5] = 0.0; J[2,6] = 0.0; J[2,7] = 0.0; J[2,8] = 0.0;
        J[3,1] = 0.0; J[3,2] = 0.0; J[3,3] = 0.0; J[3,4] = 0.0;
        J[3,5] = 0.0; J[3,6] = 0.0; J[3,7] = 0.0; J[3,8] = 0.0;
        J[4,1] = 0.0; J[4,2] = 0.0; J[4,3] = 0.0; J[4,4] = 0.0;
        J[4,5] = 0.0; J[4,6] = 0.0; J[4,7] = 0.0; J[4,8] = 0.0;
        J[5,1] = 0.0; J[5,2] = 0.0; J[5,3] = 0.0; J[5,4] = 0.0;
        J[5,5] = 0.0; J[5,6] = 0.0; J[5,7] = 0.0; J[5,8] = 0.0;
        J[6,1] = 0.0; J[6,2] = 0.0; J[6,3] = 0.0; J[6,4] = 0.0;
        J[6,5] = 0.0; J[6,6] = 0.0; J[6,7] = 0.0; J[6,8] = 0.0;
        J[7,1] = 0.0; J[7,2] = 0.0; J[7,3] = 0.0; J[7,4] = 0.0;
        J[7,5] = 0.0; J[7,6] = 0.0; J[7,7] = 0.0; J[7,8] = 0.0;
        J[8,1] = 0.0; J[8,2] = 0.0; J[8,3] = 0.0; J[8,4] = 0.0;
        J[8,5] = 0.0; J[8,6] = 0.0; J[8,7] = 0.0; J[8,8] = 0.0;
    end

    @inbounds begin
        # Row 1: f1 = a1 + x3^2 - x4^2 + b1*(x1 + x4)
        J[1,1] = b1
        J[1,3] = 2*x3
        J[1,4] = b1 - 2*x4

        # Row 2: f2 = 2*x3*x4 + b1*(x2 - x3)
        J[2,2] = b1
        J[2,3] = -b1 + 2*x4
        J[2,4] = 2*x3

        # Row 3: f3 = a2 + x1^2 - x2^2 + b2*(x3 + x2) + j1*x5
        J[3,1] = 2*x1
        J[3,2] = b2 - 2*x2
        J[3,3] = b2
        J[3,5] = j1

        # Row 4: f4 = 2*x1*x2 + b2*(x4 - x1) + j1*x6
        J[4,1] = -b2 + 2*x2
        J[4,2] = 2*x1
        J[4,4] = b2
        J[4,6] = j1

        # Row 5: f5 = a4 + x7^2 - x8^2 + b1*(x5 + x8)
        J[5,5] = b1
        J[5,7] = 2*x7
        J[5,8] = b1 - 2*x8

        # Row 6: f6 = 2*x7*x8 + b1*(x6 - x7)
        J[6,6] = b1
        J[6,7] = -b1 + 2*x8
        J[6,8] = 2*x7

        # Row 7: f7 = a3 + x5^2 - x6^2 + b3*(x7 + x6) + j2*x1
        J[7,1] = j2
        J[7,5] = 2*x5
        J[7,6] = b3 - 2*x6
        J[7,7] = b3

        # Row 8: f8 = 2*x5*x6 + b3*(x8 - x5) + j2*x2
        J[8,2] = j2
        J[8,5] = -b3 + 2*x6
        J[8,6] = 2*x5
        J[8,8] = b3
    end

    return nothing
end

function jacobian(u::AbstractVector, p)
    J = zeros(8, 8)
    jacobian!(J, u, p)
    return J
end


############################
# 3. Augmented system: u + Q
############################

"""
    flow_tangent!(dY, Y, p, t)

Augmented ODE for base state u(t) and tangent matrix Q(t):

    dot(u) = f(u)
    dot(Q) = J(u) * Q

State layout:
- Y[1:N]        = u
- Y[N+1:end]    = vec(Q) in column-major order (length N * pexp)
"""
function flow_tangent!(dY, Y, p, t)
    n = 8
    total_len = length(Y)
    pexp = (total_len - n) ÷ n

    u   = @view Y[1:n]
    du  = @view dY[1:n]
    Qv  = @view Y[n+1:total_len]
    dQv = @view dY[n+1:total_len]

    Q  = reshape(Qv,  n, pexp)
    dQ = reshape(dQv, n, pexp)

    # base dynamics
    rhs!(du, u, p, t)

    # Jacobian * Q
    J = zeros(n, n)          
    jacobian!(J, u, p)
    mul!(dQ, J, Q)

    return nothing
end


############################
# 4. Lyapunov spectrum via Jacobian + QR
############################

"""
    lyapunov_spectrum(p;
                      u0=nothing,
                      Ttr=500.0,
                      T_LLE=500.0,
                      Δt_renorm=0.1,
                      pexp=N,
                      reltol=1e-8,
                      abstol=1e-10,
                      upper_bound=1e8)

Compute up to `pexp` Lyapunov exponents for parameters `p`.

Returns:
- λ :: Vector{Float64} of length `pexp` on success
- `nothing` if trajectory blows up (norm(u) > upper_bound or NaN)

Arguments:
- p: parameter tuple (a1,a2,a3,a4,b1,b2,b3,j1,j2)
- u0: optional initial condition (Vector{Float64} of length 8). If `nothing`,
      a random normalized vector is used.
- Ttr: transient time before LLE accumulation
- T_LLE: time over which exponents are averaged
- Δt_renorm: physical time between QR renormalizations
- pexp: number of exponents to compute (1 ≤ pexp ≤ N)
"""
function lyapunov_spectrum(p;
                           u0=nothing,
                           Ttr=500.0,
                           T_LLE=500.0,
                           Δt_renorm=0.1,
                           pexp::Int=N,
                           reltol=1e-8,
                           abstol=1e-10,
                           upper_bound=1e8)

    n = N
    @assert 1 ≤ pexp ≤ n "pexp must be between 1 and $n"

    # -------------------------
    # (1) Initial condition
    # -------------------------
    if u0 === nothing
        u0 = randn(n)
        u0 ./= norm(u0)
    else
        u0 = copy(u0)
    end

    # -------------------------
    # (2) Transient integration of base system
    # -------------------------    
    prob_tr = ODEProblem(rhs!, u0, (0.0, Ttr), p)
    sol_tr = solve(prob_tr, Tsit5();
               reltol        = reltol,
               abstol        = abstol,
               maxiters      = 5_000_000,    # <-- this is the key
               save_everystep = false)
 #   
    if sol_tr.retcode != :Success
        return nothing  # treat as failure for this (j1, j2)
    end

    u_tr = Array(sol_tr.u[end])
    if !all(isfinite, u_tr) || norm(u_tr) > upper_bound
        return nothing
    end
 
    u_tr = Array(sol_tr.u[end])
    if !all(isfinite, u_tr) || norm(u_tr) > upper_bound
        return nothing
    end

    # -------------------------
    # (3) Initialize tangent matrix Q0
    # -------------------------
    # Use first pexp columns of identity as initial basis
    Q0 = Matrix(I, n, n)[:, 1:pexp]

    # Augmented initial state Y0 = [u_tr; vec(Q0)]
    Y = [u_tr; vec(Q0)]

    # -------------------------
    # (4) Main LLE accumulation loop
    # -------------------------
    S = zeros(pexp)     # running sum of log stretch
    total_time = 0.0

    while total_time < T_LLE
        dt_seg = min(Δt_renorm, T_LLE - total_time)
        
        prob = ODEProblem(flow_tangent!, Y, (0.0, dt_seg), p)
        sol  = solve(prob, Tsit5();
                 reltol        = reltol,
                 abstol        = abstol,
                 maxiters      = 5_000_000,  # <-- and here
                 save_everystep = false)

        Y = Array(sol.u[end])
        u = @view Y[1:n]

        if sol.retcode != :Success
            return nothing  # this segment failed -> LLE not reliable
        end

        # blow-up guard
        if !all(isfinite, u) || norm(u) > upper_bound
            return nothing
        end

        # Extract Q and perform QR
        Qv = @view Y[n+1:end]
        Q  = reshape(Qv, n, pexp)

        F = qr(Q)  # thin QR
        Q_orth = Matrix(F.Q)[:, 1:pexp]
        R      = F.R[1:pexp, 1:pexp]

        diagR = abs.(diag(R))
        # Guard against degenerate directions
        diagR .= max.(diagR, eps())

        S .+= log.(diagR)

        # Rebuild Y with orthonormal Q
        Qv .= vec(Q_orth)
        Y[1:n] .= u

        total_time += dt_seg
    end

    λ = S ./ total_time
    return λ
end

export rhs!, jacobian, flow_tangent!, lyapunov_spectrum

end # module
