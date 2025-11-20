module MySystem

export integrate, rhs, rhs2, rk
using LinearAlgebra: norm
using StaticArrays

function rhs(a1::Float64, a2::Float64, a3::Float64, b1::Float64, b2::Float64, b3::Float64, j1::Float64, j2::Float64)
    function f(x)
	@SVector [a1+(x[3]^2-x[4]^2)+b1*(x[1]+x[4]),2*x[3]*x[4]+b1*(x[2]-x[3]),a2+(x[1]^2-x[2]^2)+b2*(x[3]+x[2])+j1*x[5],2*x[1]*x[2]+b2*(x[4]-x[1])+j1*x[6],a1+(x[7]^2-x[8]^2)+b1*(x[5]+x[8]),2*x[7]*x[8]+b1*(x[6]-x[7]),a3+(x[5]^2-x[6]^2)+b3*(x[7]+x[6])+j2*x[1],2*x[5]*x[6]+b3*(x[8]-x[5])+j2*x[2]]
    end
end

function rhs2(a::Float64, b::Float64, c1::Float64, c2::Float64)
    function f(x)
	@SVector [a+(x[3]^2-x[4]^2)+c1*(x[1]+x[4]),2*x[3]*x[4]+c1*(x[2]-x[3]),b+(x[1]^2-x[2]^2)+c2*(x[3]+x[2]),2*x[1]*x[2]+c2*(x[4]-x[1])]
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


# Define the modified integrate function
function integrate(f,size; tol = 1e-3, tend = 2000, upper_bound = 1e4)
    # Initialization
    y = randn(size)
    t = 0
    h = 1e-6
    fy = f(y)
    flag = 0
    # Initialize time series
    time_series = [[t, y...,zeros(size)...]]
    
    # Adaptive step size RK integration
    while t < tend
        h = min(h, tend - t)  # Adjust step size for the last step
        
        # RK integration step
        (z, fz, est) = rk(y, f, fy, h)
        
        # Update time and state
        t += h
        y, fy = z, fz
         
        # Store the time and state in the time series
        push!(time_series, [t, y..., fy...])
        
        # Check if the norm of the state exceeds the upper bound
        if norm(y) > upper_bound
            flag = 1
            break
        end

        # Adjust step size for the next iteration
        fac = (tol / est)^(1 / 5)
        h = min(2, max(0.5, 0.8 * fac)) * h
    end
    
    return time_series, flag
end

end
