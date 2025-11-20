clc; close all; clear all;

x0 = 2*rand(3,1)-1; % Initial condition
tdisp = linspace(0,200,10000);
[t, x] = ode45(@mySystem, tdisp,x0);

figure(1);
plot3(x(1:end-20,1),x(11:end-10,1),x(21:end,1));

lyapunovExponent(x(:,1),1, 'Dimension',3,'lag',10, 'ExpansionRange',200);

%disp(['Maximum Lyapunov exponent: ', num2str(l)]);

% Define the multivariate dynamical system (modify this part according to your system)
function dxdt = mySystem(t, x)
    % Example system: Lorenz system
    sigma = 10;
    rho = 28;
    beta = 8/3;
    dxdt = [sigma * (x(2) - x(1));
            x(1) * (rho - x(3)) - x(2);
            x(1) * x(2) - beta * x(3)];
end