function eigenvalues = evalsjacLMW(values,c,a)
    % values is a vector of size 8

    A = [-1,0,2*values(3),-2*values(4)-1,0,0,0,0;
         0,-1,2*values(4)+1,2*values(3),0,0,0,0;
         2*values(1),-2*values(2)-1.8,-1.8,0,c,0,0,0;
         2*values(2)+1.8,2*values(1),0,-1.8,0,c,0,0;
         0,0,0,0,-1,0,2*values(7),-2*values(8)-1;
         0,0,0,0,0,-1,2*values(8)+1,2*values(7);
         a,0,0,0,2*values(5),-2*values(6)-1.8,-1.8,0;
         0,a,0,0,2*values(6)+1.8,2*values(5),0,-1.8];

    % Calculate eigenvalues
    eigenvalues = eig(A);
end