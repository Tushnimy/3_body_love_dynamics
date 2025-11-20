using BifurcationKit, Parameters

function LMW(X,p)
    @unpack a1, a2, a3, b1, b2, b3, j1, j2 = p
    x1, x2, x3, x4, x5, x6, x7, x8 = X
    [
        a1+(x3^2-x4^2)+b1*(x1+x4)
        2*x3*x4+b1*(x2-x3)
        a2+(x1^2-x2^2)+b2*(x3+x2)+j1*x5
        2*x1*x2+b2*(x4-x1)+j1*x6
        a1+(x7^2-x8^2)+b1*(x5+x8)
        2*x7*x8+b1*(x6-x7)
        a3+(x5^2-x6^2)+b3*(x7+x6)+j2*x1
        2*x5*x6+b3*(x8-x5)+j2*x2
    ]
end

par_sl = (a1=-1,a2=1,a3=1,b1=-1,b2=-1.8,b3=-1.8,j1=0,j2=0)
u0 = [0.5419,-1.6723,1.5112,-1.0533,-0.5419,-1.6723,1.5112,-1.0533]
prob = BifurcationProblem(LMW,u0,par_sl,(@lens _.j1))
opts = ContinuationPar()
br = continuation(prob, PALC(),opts, bothside=true)

