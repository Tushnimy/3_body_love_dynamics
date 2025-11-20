using DynamicalSystems
using Plots

function LMW_rule!(du, u, p, t)
    a1 =p[1]
    a2 =p[2]
    a3 =p[3]
    b1 =p[4]
    b2 =p[5]
    b3 =p[6]
    j1 =p[7]
    j2 =p[8]

    du[1] = a1+(u[3]^2-u[4]^2)+b1*(u[1]+u[4])
    du[2] = 2*u[3]*u[4]+b1*(u[2]-u[3])

    du[3] = a2+(u[1]^2-u[2]^2)+b2*(u[3]+u[2])+j1*u[5]
    du[4] = 2*u[1]*u[2]+b2*(u[4]-u[1])+j1*u[6]

    du[5] = a1+(u[7]^2-u[8]^2)+b1*(u[5]+u[8])
    du[6] = 2*u[7]*u[8]+b1*(u[6]-u[7])

    du[7] = a3+(u[5]^2-u[6]^2)+b3*(u[7]+u[6])+j2*u[1]
    du[8] = 2*u[5]*u[6]+b3*(u[8]-u[5])+j2*u[2]
    
    return nothing
end

thresh = 1e-7
nd = 100
M = zeros(nd,nd)
X = range(-3, stop=3, length=nd)
Y = range(-3, stop=3, length=nd)

u1 = [0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0]
p1 = [-1.0, 1.0, 1.0, -1.0, -1.8, -1.8,-3,-3]
LMW = ContinuousDynamicalSystem(LMW_rule!,u1,p1)
for (i,x) = enumerate(X)
    for (j,y) = enumerate(Y)
        if j==nd || i ==nd
            break
        ends
        print(i,' ',j,'\n')
        set_parameter!(LMW,7,x)
        set_parameter!(LMW,8,y)
        lam = lyapunov(LMW,10000;Ttr=90000,Δt = 0.01)

        if lam > thresh
            M[j,i] = 3 #Chaotic
        end
        if (abs(lam) <=thresh)
            M[j,i] = 1 #Periodic
        end
        if (lam < -1*thresh)
            M[j,i] = 0 #FP
        end
    end
end

p = heatmap(X,Y,M, xlabel="j1",ylabel="j2", size=(800,800))
print("Done")
display(p)
