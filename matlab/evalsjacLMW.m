function eigenvalues = evalsjacLMW(values,param)
    % values is a vector of size 8
    a1 = param(1);
    a2 = param(2);
    a3 = param(3);
    a4 = param(4);
    b1 = param(5);
    b2 = param(6);
    b3 = param(7);
    j1 = param(8);
    j2 = param(9);

    A = [b1,0,2*values(3),-2*values(4)+b1,0,0,0,0;
         0,b1,2*values(4)-b1,2*values(3),0,0,0,0;
         2*values(1),-2*values(2)+b2,b2,0,j1,0,0,0;
         2*values(2)-b2,2*values(1),0,b2,0,j1,0,0;
         0,0,0,0,b1,0,2*values(7),-2*values(8)+b1;
         0,0,0,0,0,b1,2*values(8)-b1,2*values(7);
         j2,0,0,0,2*values(5),-2*values(6)+b3,b3,0;
         0,j2,0,0,2*values(6)-b3,2*values(5),0,b3];

    % Calculate eigenvalues
    eigenvalues = eig(A);
end