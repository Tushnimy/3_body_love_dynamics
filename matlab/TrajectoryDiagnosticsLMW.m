clear; close all; clc

%% -------- User choices --------
% Choose your parameter pair (j1, j2):
j1 = -1.5569620253164558;   
j2 = 0.6455696202531646;

% Integration settings
Tfinal = 1e5;           % total integration time
t_span = [0, Tfinal];

% For FFT: choose a uniform sampling step
dt = 5;                 % adjust as needed (tradeoff: resolution vs speed)

% Initial condition (8D)
rng(1);
x0 = -1 + 2.*rand(1, 8);


%% -------- Base parameters --------
% param = [a; b; ???; ???; c1; c2; j1; j2]  (your original ordering)
param_base = [-1; 1.0; 1.0; -1.0; 0.0; -1.0; -1.0; 0; 0];
param = param_base;
param(8) = j1;
param(9) = j2;

%% -------- Integrate system --------
fprintf('Integrating three_lovers with j1 = %.6f, j2 = %.6f ...\n', j1, j2);
opts   = odeset('RelTol',1e-6,'AbsTol',1e-8,'MaxStep',1);
[t, x] = ode45(@(t,x) three_lovers(t, x, param), t_span, x0);

% Extract Layla's and Majnun's loves:
Lr = x(:,1);   % Layla real
Li = x(:,2);   % Layla imag
Mr = x(:,3);   % Majnun real (assumed)
Mi = x(:,4);   % Majnun imag (assumed)

%% -------- Trim transient (for diagnostics) --------
fracTransient = 0.9;    
N = numel(t);
idx_tr = ceil(fracTransient * N) : N;

t_tr  = t(idx_tr);
Lr_tr = Lr(idx_tr);
Li_tr = Li(idx_tr);
Mr_tr = Mr(idx_tr);
Mi_tr = Mi(idx_tr);

N_tr = numel(t_tr);

%% -------- Resample onto uniform grid for FFT --------
t_uniform = (t_tr(1):dt:t_tr(end)).';   % uniform grid on trimmed interval
Lr_uniform = interp1(t_tr, Lr_tr, t_uniform, 'pchip');
Li_uniform = interp1(t_tr, Li_tr, t_uniform, 'pchip');

% Remove mean for FFT clarity
Lr_uniform_m = Lr_uniform - mean(Lr_uniform);

% Sampling frequency
fs = 1/dt;

%% -------- Power spectrum via Welch --------
% Use Welch for smooth spectrum
nfft  = 2^nextpow2(numel(Lr_uniform_m));
win   = floor(numel(Lr_uniform_m) / 8);
if mod(win, 2) == 1
    win = win + 1;
end
if win < 32
    win = min(numel(Lr_uniform_m), 256);
end
noverlap = floor(0.5 * win);

[Pxx, f] = pwelch(Lr_uniform_m, win, noverlap, nfft, fs);

%% -------- Poincaré section (Mr vs Mi at Lr = const) --------
[xP, yP] = poincare_section(Lr_tr, Mr_tr, Mi_tr);

%% -------- Plotting: publication-style diagnostics --------
fig = figure(1); clf(fig);
set(fig, 'Color', 'w', 'Position', [100 100 1200 800]);

tiledlayout(2,2, 'TileSpacing', 'compact', 'Padding', 'compact');

paramStr = sprintf('j_1 = %.4f, j_2 = %.4f', j1, j2);
sgtitle(sprintf('Diagnostics for three lovers: %s', paramStr), ...
    'FontName', 'Times New Roman', 'FontWeight', 'bold', 'FontSize', 14);

%% --- (1) L_r(t) vs time (tail only)
nexttile;
plot(t_tr(ceil(0.9 * N_tr):end), Lr_tr(ceil(0.9 * N_tr):end), 'LineWidth', 0.8);
xlabel('t', 'FontName', 'Times New Roman');
ylabel('L_r(t)', 'FontName', 'Times New Roman');
title('Time series of L_r(t)', 'FontName', 'Times New Roman');
grid on; box on;
set(gca, 'LineWidth', 0.75, 'FontName', 'Times New Roman', 'FontSize', 11);

L = numel(t);
frac = 0.0;
start = ceil(frac*N);

%% --- (2) Trajectory in (L_r, L_i)
nexttile;
plot(Lr(start+1:end), Li(start+1:end), 'LineWidth', 0.3);
xlabel('L_r', 'FontName', 'Times New Roman');
ylabel('L_i', 'FontName', 'Times New Roman');
title('Trajectory in (L_r, L_i)', 'FontName', 'Times New Roman');
grid on; box on;
axis tight;
set(gca, 'LineWidth', 0.75, 'FontName', 'Times New Roman', 'FontSize', 11);

%% --- (3) Poincaré section (M_r vs M_i at L_r = const)
nexttile;
if ~isempty(xP)
    scatter(xP, yP, 10, 'filled');
    xlabel('M_r (normalised) at L_r = const', 'FontName', 'Times New Roman');
    ylabel('M_i (normalised) at L_r = const', 'FontName', 'Times New Roman');
    title('Poincaré section (M_r vs M_i)', 'FontName', 'Times New Roman');
    grid on; box on;
    axis tight;
else
    text(0.5, 0.5, 'No valid Poincaré crossings found', ...
        'HorizontalAlignment', 'center', 'FontName', 'Times New Roman');
    axis off;
end
set(gca, 'LineWidth', 0.75, 'FontName', 'Times New Roman', 'FontSize', 11);

%% --- (4) Power spectrum of L_r
nexttile;
plot(f, Pxx, 'LineWidth', 0.3);
xlim([0, fs/2]);
set(gca, 'YScale', 'log');  % log scale usually looks nicer
xlabel('Frequency', 'FontName', 'Times New Roman');
ylabel('Power', 'FontName', 'Times New Roman');
title('Power spectrum of L_r(t)', 'FontName', 'Times New Roman');
grid on; box on;
set(gca, 'LineWidth', 0.75, 'FontName', 'Times New Roman', 'FontSize', 11);

% Optional: save as vector PDF for publication
% exportgraphics(fig, sprintf('diagnostics_j1_%.3f_j2_%.3f.pdf', j1, j2), ...
%     'ContentType', 'vector');

max_vals = max(x, [], 1);
min_vals = min(x, [], 1);
for k = 1:8
    fprintf("Observable y%d: min = %g, max = %g\n", ...
            k, min_vals(k), max_vals(k));
end

%% ------------- Local function(s) ------------- %%
function [xP, yP] = poincare_section(xObs, x2, x3)
%POINCARE_SECTION Simple Poincaré section on xObs = const, d/dt > 0.
%
%   xObs : observable used to define the section (vector)
%   x2,x3: coordinates to be reported at section crossings
%
%   Returns:
%   xP, yP: points (x2, x3) at Poincaré crossings.

    xObs = xObs(:);
    x2   = x2(:);
    x3   = x3(:);

    Nloc = length(xObs);
    if any([length(x2), length(x3)] ~= Nloc)
        error('xObs, x2, x3 must have same length.');
    end

    % Define section at mean of xObs in this window
    xSec = mean(xObs);

    % Approximate derivative sign using forward difference
    dx   = diff(xObs);

    xP = [];
    yP = [];

    for k = 1:Nloc-1
        % Crossing of xObs = xSec with positive slope
        if (xObs(k) < xSec && xObs(k+1) >= xSec && dx(k) > 0)
            % Linear interpolation for better accuracy
            alpha = (xSec - xObs(k)) / (xObs(k+1) - xObs(k));  % in [0,1]
            x2_cross = x2(k) + alpha * (x2(k+1) - x2(k));
            x3_cross = x3(k) + alpha * (x3(k+1) - x3(k));

            xP(end+1,1) = x2_cross; %#ok<AGROW>
            yP(end+1,1) = x3_cross; %#ok<AGROW>
        end
    end
end


%P: (4, -2)/ (4, -1.9966555183946488)/ (-2,2.5)