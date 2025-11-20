clear all, close all, clc
param1 = [-1;1;1.0;0.75];
param2 = [-1;1;0;-1.8];
param3 = [-1;1;0.2;-1.8];
param4 = [-1;1;0.4;-1.8];

x0 = -1+2.*rand(8,4);

for k = 1:8
    [t,x] = ode45(@(t,x)layla_majnu(t,x,param1),[1 100000],x0(k,:));
    figure(1) 
    subplot(2,1,1) 
    %plot(t(t>90000),x((t>90000),1));
    plot(t,x(:,1));
    title('Laylas love') 
    xlabel 'time'; 
    ylabel 'L_r';
    hold on
    subplot(2,1,2) 
    %plot(x((t>90000),1), x((t>90000),2));
    plot(x(:,1), x(:,2));
    title('Trajectory') 
    xlabel 'L_r'; 
    ylabel 'L_i'; 
    hold on
end
hold off

