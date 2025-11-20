# Computing bifurcation diagram in Julia
# written by W. Marszalek, H. Podhaisky, and J. Sadecki (CC BY-SA 3.0)
# June 2019 

module bif
export integrate, rhs
using LinearAlgebra: norm
using StaticArrays

# The function below is the RHS function
function rhs(a1::Float64, a2::Float64, a3::Float64, b1::Float64, b2::Float64, b3::Float64, j1::Float64, j2::Float64)
    function f(x)
	@SVector [a1+(x[3]^2-x[4]^2)+b1*(x[1]+x[4]),2*x[3]*x[4]+b1*(x[2]-x[3]),a2+(x[1]^2-x[2]^2)+b2*(x[3]+x[2])+j1*x[5],2*x[1]*x[2]+b2*(x[4]-x[1])+j1*x[6],a1+(x[7]^2-x[8]^2)+b1*(x[5]+x[8]),2*x[7]*x[8]+b1*(x[6]-x[7]),a3+(x[5]^2-x[6]^2)+b3*(x[7]+x[6])+j2*x[1],2*x[5]*x[6]+b3*(x[8]-x[5])+j2*x[2]]
    end
end

#= Runge-Kutta 4th order; takes in input, output at input, function and step size; 
returns new input, error at new input and value 
=#
function rk(y, f, fy, h)
    k1 = fy
    k2 = f(y + h * 1 / 5 * k1)
    k3 = f(y + h * (3 / 40 * k1 + 9 / 40 * k2))
    k4 = f(y + h * (44 / 45 * k1 - 56 / 15 * k2 + 32 / 9 * k3))
    k5 = f(y + h * (19372 / 6561 * k1 - 25360 / 2187 * k2 + 
         64448 / 6561 * k3 - 212 / 729 * k4))
    k6 = f(y + h * (9017 / 3168 * k1 - 355 / 33 * k2 + 46732 / 5247 * k3 + 
          49 / 176 * k4 - 5103 / 18656 * k5))
    z = y + h * (35 / 384 * k1 + 500 / 1113 * k3 + 125 / 192 * k4 - 
          2187 / 6784 * k5 + 11 / 84 * k6)
    k7 = f(z)
    w = y + h * (5179 / 57600 * k1 + 7571 / 16695 * k3 + 
           393 / 640 * k4 - 92097 / 339200 * k5 + 187 / 2100 * k6 + 1 / 40 * k7)
    est = norm(z - w)
    return (z, k7, est)
end 

#=
takes in function, input, output, runge-kutta output, step size, icomp(?)
returns new step size, new input, value at new input
=#

function findevent(f, y, fy, fz, h, icomp)
    he = h * fy[icomp] / (fy[icomp] - fz[icomp])
    ye, fye, _ = rk(y, f, fy, he)
    he = he * fy[icomp] / (fy[icomp] - fye[icomp])
    ye, fye, _ = rk(y, f, fy, he)
    (he, ye, fye)
end

#=
Checks period; sets error(?) to 10^-3, 
=#

function checkperiod(E, ip, px, px21)
    eps = 1e-3
    for p = 1:px
        p1 = mod(ip - p - 1, px21) + 1 
        p2 = mod(ip - 2 * p - 1, px21) + 1
        d1 = norm(E[ip,2:end] - E[p1,2:end])
        d2 = norm(E[p1,2:end] - E[p2,2:end])
        t1 = E[ip,1] - E[p1,1]
        t2 = E[p1,1] - E[p2,1]
        if (d1 < eps) && (d2 < eps) && (abs(t1 - t2) < eps)
            return p
        end
    end
    return 0
end

function integrate(f; tol = 1e-6, tend = 2000, px = 16, plotting = false)
# returns the type of orbit, the maxima of the second component, the solution and 
# the event data 
    y = [0.4,0.9,0.3,0.1,0.3,0.3,0.8,0.24]
    t = 0; h = 1e-6; hnew = h; fy = f(y);
    icomp = 2
    px21 = 2 * px + 1
    E = zeros((px21, 9))
    ip = 1
    period = 0
    stat = Dict(:accept => 0, :reject => 0)
    L = []
    while (t < tend) && (norm(y) < 1e4)
        h = hnew
        (z, fz, est) = rk(y, f, fy, h)
        fac = (tol / est)^(1 / 5)
        hnew = min(2, max(0.5, 0.8 * fac)) * h
        if fac < 1 
            stat[:reject] += 1
            continue
        else
            if plotting
                push!(L, [t, y[1], y[2], y[3], y[4], y[5], y[6], y[7], y[8]])
            end
            if  (fy[icomp] > 0) && (fz[icomp] < 0)
                he, ye, fye = findevent(f, y, fy, fz, h, icomp)
                E[ip,1] = t + he
                E[ip,2:9] = ye
                period = checkperiod(E, ip, px, px21)
                if period > 0
                    break
                end
                ip = mod(ip, px21) + 1
            end 
            y, fy = z, fz
            t += h
            stat[:accept] += 1
        end
    end
    y2x = Array{Float64}([E[mod(ip-k-1,px21)+1,icomp+1] for k=0:max(0,period-1)])
    if norm(y)>=1e4
        period = -1
    end
    period, y2x , L, E
end
end