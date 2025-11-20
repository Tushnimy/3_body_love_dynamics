include("LaylaMajnuWardAll.jl")
using Plots, Main.bif, NPZ
plotly()

nd = 200 #200
M = zeros(nd,nd)
X = range(-3, stop=3, length=nd)
Y = range(-3, stop=3, length=nd)
constant_params = (-1.0, 1.0, 1.0, -1.0, -1.8, -1.8)
constp2= (-2, -0.35)
for (i,x) = enumerate(X)
    for (j,y) = enumerate(Y)
        params = (constant_params...,x,y)
        print(i)
        M[j,i] = integrate(rhs(params...))[1]
    end
end
p = heatmap(X,Y,M, xlabel="j1",ylabel="j2", size=(600,600))
plot!(p,X,X,color=:white,lw=2,label="j_1=j_2")
display(p)

#The weird dotty plot in orange purple and black