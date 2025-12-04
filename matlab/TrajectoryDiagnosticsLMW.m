clear; close all; clc

%% -------- User choices --------
% Chooseparameter pair (j1, j2):
j1 = 3.0;  
j2 = -2.7;   %Poincare for 1.27 is interesting

% Integration settings
Tfinal = 1e5;           
t_span = [0, Tfinal];
f = 0.5;

% Initial condition (8D)
%rng(1);
x0 = -1 + 2.*rand(1, 8);


%% -------- Base parameters --------
% param = [a1; a2; a3; a4; b1; b2; j1; j2]  
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
plot(t(t>0), Lr(t>0), 'LineWidth', 0.8);
xlabel('t', 'FontName', 'Times New Roman');
ylabel('L_r(t)', 'FontName', 'Times New Roman');
title('Time series of L_r(t)', 'FontName', 'Times New Roman');
grid on; box on;
set(gca, 'LineWidth', 0.75, 'FontName', 'Times New Roman', 'FontSize', 11);

%% --- (2) Trajectory in (L_r, L_i)
nexttile;
plot(Lr(t>0), Li(t>0), 'LineWidth', 0.3);
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

            xP(end+1,1) = x2_cross; 
            yP(end+1,1) = x3_cross; 
        end
    end
end


%P: (1.25, 3), (3,1.6), (3,1.625) 
%FP: (3, 3): x0 = [-3, -1, 1, 0, -3, -1, 1, 0] + f*(-1 + 2.*rand(1,8));
%C: (3, -1.0)
%UB: (3, 3)

% FP: (3, -2.8)  -0.8102   -0.2492    0.0920   -0.7766    0.8089    0.2666    0.8108    0.2611
% FP: (3, -2.7)   0.7345    0.5325   -0.7301   -0.8702    0.3449   -0.0356   -0.0092   -0.3730

