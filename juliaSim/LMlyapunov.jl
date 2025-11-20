using DynamicalSystems
using Plots

function LMW_rule!(du, u, p, t)
    a1 =p[1]
    a2 =p[2]
    b1 =p[3]
    b2 =p[4]

    du[1] = a1+(u[3]^2-u[4]^2)+b1*(u[1]+u[4])
    du[2] = 2*u[3]*u[4]+b1*(u[2]-u[3])
    du[3] = a2+(u[1]^2-u[2]^2)+b2*(u[3]+u[2])
    du[4] = 2*u[1]*u[2]+b2*(u[4]-u[1])
    return nothing
end


nd = 50
M = zeros(nd,nd)
X = range(-1, stop=1, length=nd)
Y = range(-1, stop=1, length=nd)

u1 = [0.0,0.0,0.0,0.0]
p1 = [-1.0, 1.0,-1,-1]
LMW = ContinuousDynamicalSystem(LMW_rule!,u1,p1)
for (i,x) = enumerate(X)
    for (j,y) = enumerate(Y)
        if (x+y >=0)
            break
        end
        print(i,' ',j,'\n')
        set_parameter!(LMW,3,x)
        set_parameter!(LMW,4,y)
        M[j,i] = lyapunov(LMW,3000;Ttr=50000,Δt = 0.1)
    end
end

p = heatmap(X,Y,M, xlabel="j1",ylabel="j2", size=(800,800))
