function dx = LaylaMajnuWard(t,x,param)
dx = [
    param(1)+(x(3)^2-x(4)^2)+param(4)*(x(1)+x(4));
    2*x(3)*x(4)+param(4)*(x(2)-x(3));
    param(2)+(x(1)^2-x(2)^2)+param(5)*(x(3)+x(2)) +param(7)*x(5);
    2*x(1)*x(2)+param(5)*(x(4)-x(1))+param(7)*x(6);
    param(1)+(x(7)^2-x(8)^2)+param(4)*(x(5)+x(8));
    2*x(7)*x(8)+param(4)*(x(6)-x(7));
    param(3)+(x(5)^2-x(6)^2)+param(6)*(x(7)+x(6))+param(8)*x(1);
    2*x(5)*x(6)+param(6)*(x(8)-x(5))+param(8)*x(2);
];