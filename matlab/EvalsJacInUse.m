function eigenvalues = EvalsJac(values,c)
    % values is a vector [x, y, z, w]

    % Example 4x4 matrix with elements as functions of x, y, z, w
    A = [c,        0,        2*values(3),        -2*values(4)+c;
         0,        c,      2*values(4)-c,           2*values(3);
         2*values(1), -2*values(2)-1.8, -1.8,                 0;
         2*values(2)+1.8,    2*values(1),  0,              -1.8];

    % Calculate eigenvalues
    eigenvalues = eig(A);
end