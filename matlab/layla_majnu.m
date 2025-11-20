function dx = layla_majnu(t,x,param)
dx = [
    param(1)+x(3)^2-x(4)^2+param(3)*(x(1)+x(4));
    2*x(3)*x(4)+param(3)*(x(2)-x(3));
    param(2)+x(1)^2-x(2)^2+param(4)*(x(3)+x(2));
    2*x(1)*x(2)+param(4)*(x(4)-x(1));
];