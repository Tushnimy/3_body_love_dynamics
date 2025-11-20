close all; clc;close all;
t=50;
r = linspace(-3,3,t);
p = zeros(16,t);
Z=zeros(t);
syms x y z w 
sols=[];
a=2.5;
for i=1:length(r)
    c=r(i);
    eq1 =-1+y^2-1*(x-1i*y);
    eq2 =1+x^2-1.8*(y-1i*x)+c*z;
    eq3 =-1+w^2-1*(z-1i*w);
    eq4 =1+z^2-1.8*(w-1i*z)+a*x; 
    [solx,soly,solz,solw] = solve([eq1,eq2,eq3,eq4],[x,y,z,w]);
    solutions = [vpa(solx), vpa(soly), vpa(solz), vpa(solw)];
    sols=[real(solutions(:,1)),imag(solutions(:,1)),real(solutions(:,2)),imag(solutions(:,2)),real(solutions(:,3)),imag(solutions(:,3)),real(solutions(:,4)),imag(solutions(:,4))];
    for k=1:16
        eigenvalues = evalsjacLMW(sols(k,:),c,a);
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
plot(r,Z);
hold on
plot(Z,r);
hold on;
plot([-3,-2.02041],[0,0],'LineWidth',1,'Color','green');
hold on;
plot([-2.02041,1.49],[0,0],'LineWidth',1,'Color','blue');
hold on
plot([1.49,2.632],[0,0],'LineWidth',1,'Color','green');
hold on
plot([2.632,3],[0,0],'LineWidth',1,'Color','blue');
hold on
ylabel('max(Re(\lambda))');
xlabel('j_1')
title('Eigenvalue plot');