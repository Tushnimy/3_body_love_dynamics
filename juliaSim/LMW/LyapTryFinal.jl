module lyapLMW

export rhs3, rk, integrate_steps, LLE_rk4, rhs

using Plots
using LinearAlgebra
using StaticArrays 


function rhs3(p)
    a1 =p[1]
    a2 =p[2]
    a3 =p[3]
    b1 =p[4]
    b2 =p[5]
    b3 =p[6]
    j1 =p[7]
    j2 =p[8]

    function f(u)
    @SVector [a1+(u[3]^2-u[4]^2)+b1*(u[1]+u[4]),2*u[3]*u[4]+b1*(u[2]-u[3]),a2+(u[1]^2-u[2]^2)+b2*(u[3]+u[2])+j1*u[5],2*u[1]*u[2]+b2*(u[4]-u[1])+j1*u[6],a1+(u[7]^2-u[8]^2)+b1*(u[5]+u[8]),2*u[7]*u[8]+b1*(u[6]-u[7]),a3+(u[5]^2-u[6]^2)+b3*(u[7]+u[6])+j2*u[1],2*u[5]*u[6]+b3*(u[8]-u[5])+j2*u[2]]
    end
end

function rhs(a::Float64, b::Float64, c1::Float64, c2::Float64)
    function f(x)
	@SVector [a+(x[3]^2-x[4]^2)+c1*(x[1]+x[4]),2*x[3]*x[4]+c1*(x[2]-x[3]),b+(x[1]^2-x[2]^2)+c2*(x[3]+x[2]),2*x[1]*x[2]+c2*(x[4]-x[1])]
    end
end

function rhs2(p)
    sigma = p[1]
    rho = p[2]
    beta = p[3]
    function lorenz(u)
    @SVector [sigma*(u[2]-u[1]),u[1]*(rho-u[3])-u[2],u[1]*u[2]-beta*u[3]]
    end
end

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

function rk4(y,f,fy,h)
    k1 = fy
    k2 = f(y + h* k1 / 2)
    k3 = f(y + h * k2 /2)
    k4 = f(y + h*k3)
    ynew = y + (h/6)*(k1 + 2*k2 + 2*k3 + k4)
    return (ynew, f(ynew))
end

function integrate_steps(f, steps, h, u0)
    u = u0
    fu = f(u0)
    for step = 1:steps
        u, fu = rk4(u,f,fu,h)
    end
    return u, fu
end


function integrate(f, tend, u0; tol=1e-3, upper_bound=1e8)
    u = u0
    t=0
    h = 1e-6
    fu = f(u)
    flag=0
    # adaptive step size rk integration
    while t < tend
        h = min(h, tend - t)
        z,fz,est = rk(u,f,fu,h)

        t+= h
        u,fu = z, fz
        if norm(u) > upper_bound
            flag=1
            break
        end

        fac = (tol/est) ^ (1/5)
        h = min(2,max(0.5,0.8*fac))*h
    end
    return u,flag
end

function LLE_rk4(f, u0, iter, steps,h, size, Ttr)
    ur, flag = integrate(f,Ttr,u0)
    if flag==1
        return "unbound"
    end
    dev = randn(size)
    ut = ur+ (dev/norm(dev))*1e-8
    l = []
    for step=1:iter
        ur, _ = integrate_steps(f,steps,h,ur)
        ut, _ = integrate_steps(f,steps,h,ut)
        dir = (ut-ur)/ norm(ut-ur)
        push!(l,norm(ut-ur)/1e-8)
        ut = ur + dir*(1e-8)
    end
    l = log.(l)
    return sum(l)/(iter*h*steps)
end

end