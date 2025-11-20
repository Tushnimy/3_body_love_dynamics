function dx = three_lovers(t, x, param)
% Environmental factors
a1=param(1);
a2=param(2);
a3=param(3);
% Memory factors
b1=param(4);
b2=param(5);
b3=param(6);
% Jealousy factors
j1=param(7); 
j2=param(8);
% DEs
dx =[
    %Layla: Same as the memory dynamics equations
    a1+(x(3)^2-x(4)^2)+b1*(x(1)+x(4));
    2*x(3)*x(4) + b1*(x(2)-x(3));
    %Majnu: Same as memory dynamics but with added jealousy terms
    a2+(x(1)^2-x(2)^2)+b2*(x(3)+x(2))+j1*x(5);
    2*x(1)*x(2)+b2*(x(4)-x(1))+j1*x(6);
    % Layla to Ward: Same as memory dynamics but with added jealousy terms
    a1+(x(7)^2-x(8)^2) + b1*(x(5)+x(8))+ j2*(x(1));
    2*x(7)*x(8)+ b1*(x(6)-x(7)) + j2*(x(2));
    % Ward to Layla
    a3+(x(5)^2-x(6)^2)+b3*(x(7)+x(6))+j2*x(5);
    2*x(5)*x(6)+b3*(x(8)-x(5))+j2*x(6);

];