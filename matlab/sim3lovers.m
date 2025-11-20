clear all, close all, clc
% Choosing parameter values;
%{
We keep the values of a,b,c2 fixed as in the last paper (-1,1,-1.8)
We fix d = 0.5, he has the tendency to love Layla as she is his wife but
it's not as strong as Majnu.
c1 fixed at 0; no real reason, showed interesting dynamics last time
NOW; we come across the new parameters added to the system;
c3,g1,g2,e,f
e and f shall be fixed at 100,100(?) as they are there to merely act as
constraints
%}
%J=-1;
%J=-0.9;
%J=-0.8;
%J=-0.7;
%J=-0.6;
%J=-0.5;
%J=-0.4;
%J=-0.3;
%J=-0.2;
%J=-0.1;
%J=0;
%J=0.1;
%J=0.2;
%J=0.3;
%J=0.4;
%J=0.5;
%J=0.6;
%J=0.7;
%J=0.8;
%J=0.9;
%J=1.0;
%Seems to oscillate constantly for J[-2,2], I=-0.6
%I=-1; %Several different circles
%I=-0.2 OSCILLATIOONMS SKAFKA:FKAE
%I=0; VERY interesting
%I=0.2 circles, dies
%I=0.5 kinda interesting
%I=0.75 also interesting?
%I=1, also kind of interesting?
%M=0.2; is interesting
%M=0; interesting
%M=0.2; Blows up Further also blows up
%M=-1.5; circles?
%M=-1; circles again
%SEEMS LIKE WARD WINS ALL THE TIME
%I=0.979899;
%J=1.311558;

%I=0; %FP
%J=0; %FP

I=-2.5; %Unb
J=-1; %Unb

param = [-1;1;1;-1;-1.8;-1.8;I;J];

x0 = -1+2.*rand(6,8);

for k=1:6

    [t,x] = ode45(@(t,x)three_lovers(t,x,param),[0 100000],x0(k,:));
%    [t,x] = ode45(@(t,x)three_lovers(t,x,param),[0 100000],x0);
    figure(1) 
    subplot(2,1,1) 
    plot(t((t>97300)),x((t>97300),1)); 
    title('Laylas love') 
    xlabel 'time'; 
    ylabel 'L_r';
    hold on
    subplot(2,1,2) 
    plot(x(:,1), x(:,2)); 
    title('Trajectory') 
    xlabel 'L_r'; 
    ylabel 'L_i'; 
    hold on
    saveas(figure(1),'I06J16.png');
end
hold off
