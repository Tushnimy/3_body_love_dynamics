clear; close all; clc 
%% -------- User choices -------- 
% Choose your parameter pair (j1, j2): 
j1 = 3; 
j2 = 3; 
% Integration settings 
Tfinal = 1e6; % total integration time 
t_span = [0, Tfinal]; 
% For FFT: choose a uniform sampling step 
dt = 5; 
% Initial condition (8D) 
rng(1); 
x0 = -1 + 2.*rand(1, 8); 

%% -------- Base parameters -------- 
param_base = [-1; 1; 1; -1; -1.8; -1.8; 0; 0]; 
param = param_base; 
param(7) = j1; 
param(8) = j2; 
%% -------- Integrate system -------- 
fprintf('Integrating three_lovers with j1 = %.4f, j2 = %.4f ...\n', j1, j2); 
[t, x] = ode45(@(t,x) three_lovers(t, x, param), t_span, x0); 
% Extract Layla's love (first two components) 
Lr = x(:,1); Li = x(:,2); 
%% -------- Resample onto uniform grid for FFT -------- 
t_uniform = (t(1):dt:t(end)).'; 
Lr_uniform = interp1(t, Lr, t_uniform, 'pchip'); 
% Remove mean for FFT clarity 
Lr_uniform = Lr_uniform - mean(Lr_uniform);
%% -------- FFT -------- 
N = numel(t_uniform); 
fs = 1/dt; Y = fft(Lr_uniform); 
N_half = floor(N/2); 
Y_pos = Y(1:N_half); 
f_pos = (0:N_half-1).' * (fs/N);
%% -------- Plotting -------- 
figure(1); clf
set(gcf, 'Color', 'w', 'Position', [100 100 1200 450]); 
%% --- (1) Trajectory in phase space 
subplot(1,3,1); 
plot(Lr, Li, 'LineWidth', 0.8); 
xlabel('L_r'); 
ylabel('L_i'); 
title(sprintf('Trajectory (L_r,L_i), j_1 = %.3f, j_2 = %.3f', j1, j2)); 
grid on; 
box on; 
%% --- (2) L_r(t) vs time 
N_u = numel(t_uniform); 
start_idx = ceil(0.9995 * N);
% Corrected to use ceil for proper indexing 
subplot(1,3,2); 
plot(t_uniform(start_idx:end), Lr_uniform(start_idx:end), 'LineWidth', 0.8); 
xlabel('t'); 
ylabel('L_r'); 
title('L_r(t) vs time'); 
grid on; 
box on; 
%% --- (3) FFT magnitude spectrum 
subplot(1,3,3); 
plot(f_pos, abs(Y_pos), 'LineWidth', 0.8); 
xlabel('Frequency'); ylabel('|FFT(L_r)|'); 
title('FFT magnitude spectrum'); grid on; box on; 

%AP: (-1.173913, -1.173913), (-2.979933110367893,-0.5719063545150501); 
%C: (-1, -3) 
%P: (3, 1.4548494983277591)