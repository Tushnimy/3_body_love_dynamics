iclear all; clc;close all;
%x0 = [100, 0, 0, -10];
t=500;
r = linspace(-1,1,t);
p = zeros(4,t);
z=zeros(t);
syms x y 
sols=[];

for i=1:length(r)
    c=r(i);
    eq1 =-1+y^2+c*(x-j*y);
    eq2 =1+x^2-1.8*(y-j*x);
    [solx,soly] = solve([eq1,eq2],[x,y]);
    solutions = [vpa(solx), vpa(soly)];
    sols=[real(solutions(:,1)),imag(solutions(:,1)),real(solutions(:,2)),imag(solutions(:,2))];
    for k=1:4
        eigenvalues = EvalsJacInUse(sols(k,:),c);
        Real = findMaxRealVec(eigenvalues,0);
        p(k,i)=Real;
    end
end

mat=[];
for l=1:t
    mat(l)=min(p(:,l));
end

plot(r,mat);
hold on
plot(r,z);
hold on
plot(z,r);
hold on;
plot([-1,0.206413],[0,0],'LineWidth',3,'Color','blue');
hold on;
plot([0.206413,1],[0,0],'LineWidth',3,'Color','yellow');
ylabel('max(Re(\lambda)');
xlabel('c_1')
title('Eigenvalue plot');
