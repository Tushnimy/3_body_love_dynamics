include("LyapTryFinal.jl")
using Main.lyapLMW
using LinearAlgebra
using Plots
using StaticArrays

plotly()

nd =100
M = zeros(nd,nd)
X = range(-3, stop=3, length=nd)
Y = range(-3, stop=3, length=nd)
constant_params3 = (-1.0, 1.0, 1.0, -1.0, -1.8, -1.8)
constant_params = (-1.0,1.0)
size = 8
px = []
py = []
periodicTol = 0.001
for (i,x) = enumerate(X)
    for (j,y) = enumerate(Y)
        u = randn(size)
        u = u/norm(u)
        params = (constant_params...,x,y)
        params3 = [constant_params3...,x,y]
        #LLE = LLE_rk4(rhs(params...),u,1e4,10,1e-2,size,10000)
        LLE = LLE_rk4(rhs3(params3),u,1e4,10,1e-3,size,10000)
        print(i)

        if LLE =="unbound"
            M[j,i] = 3
            continue
        end
        M[j,i] = LLE
        if abs(LLE) < periodicTol
            push!(px,x)
            push!(py,y)
        end
        #if LLE > 5*1e-2
        #    M[j,i] = 0
        #end
        #if abs(LLE) < 5*1e-2
        #    M[j,i] = 1
        #end
        #if LLE < -1*5*1e-2
        #    M[j,i] = 2
        #end
    end
end

parameterPlot = heatmap(X,Y,M, xlabel="j1",ylabel="j2", size=(600,600))
plot!(parameterPlot,X,X,color=:red,label="j_1=j_2")
scatter!(parameterPlot,px,py, color=:orange, label="periodic points")
display(parameterPlot)
#1e4 and 1e-2 give accuracy upto 1st decimal digit (0.90ish)
#1e5 and 1e-2 give 0.905ish

#Lyapunov essentially