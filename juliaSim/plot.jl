include("bif.jl")
using Plots, Main.bif
plotly()

nd = 600
M = zeros(nd,nd)
X = range(2.36, stop=2.46, length=nd)
Y = range(2.29, stop=2.39, length=nd)
for (i,x) = enumerate(X)
    for (j,y) = enumerate(Y)
        a,b,c = (1.665, x, y)
        M[j,i] = integrate(rhs(a,b,c))[1]
    end
end
M[findall(x->x==0,M)].=17
M[findall(x->x==-1,M)].=NaN
heatmap(X,Y,M, xlabel="b",ylabel="c", size=(800,800))
