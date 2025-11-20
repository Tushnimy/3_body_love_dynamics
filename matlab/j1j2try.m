close all; clc; clear all;
x0 = -1+2.*rand(1,8);

I=0;
J=1;

param = [-1;1;1;-1;-1.8;-1.8;I;J];
[t,x] = ode45(@(t,x)three_lovers(t,x,param),[0 100000],x0);

% Calculate the maximum Lyapunov exponent
n = size(x, 1); % Number of data points
m = size(x, 2); % Dimension of the system
eps = 1e-7; % Small perturbation
d0 = eps * eye(m); % Initial perturbation vector
l = 0; % Initialize Lyapunov exponent

for i = 1:n
    [~, d] = ode45(@(t,x)three_lovers(t,x,param), [0 1], x(i, :)' + d0(:, 1));
    d = d(end, :)'; % Final state after integrating for 1 unit of time
    d = d - mean(d); % Subtract the mean
    d = d / norm(d); % Normalize the vector
    d0 = d; % Update the perturbation vector
    l = l + log(norm(d) / eps); % Update the Lyapunov exponent
end
l = l / n; % Average over all data points

% Display the maximum Lyapunov exponent
disp(['Maximum Lyapunov exponent: ', num2str(l)]);




