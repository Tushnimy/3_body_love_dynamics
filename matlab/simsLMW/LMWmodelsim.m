close all; clear all; clc;
param = [-1.0;1.0;1.0;-1.0;-1.8;-1.8;2.6;2.75];
arr = 0:0.1:100000;
x0 = -1+2.*rand(1,8);

    [t,x] = ode45(@(t,x)LMWtest2(t,x,param),[0 100000],x0);
    %[t,x] = ode45(@(t,x)LMWtest2(t,x,param),[0 100000],x0);
    figure(1) 
    subplot(2,1,1)
    %plot(t,x(:,1));
    plot(t(t>90000),x(t>90000,1)); 
    title('Laylas love for Majnu') 
    xlabel 'time'; 
    ylabel 'L_11';
    hold on
    subplot(2,1,2) 
    plot(x(:,1), x(:,2)); 
    title('Trajectory') 
    xlabel 'L_1'; 
    ylabel 'L_2'; 
    hold on
    saveas(figure(1),'MajnuFails-1782.png');
hold off
