clear; close all; clc

%% -------- Parameter pairs (j1, j2) to test --------
paramPairs = [
    2.678929765886290   2.49832775919732;  % class 0 (aperiodic)
    3                   1.4548494983277591;  % class 1 (periodic)
    0.05                0;                 % class 2 (fixed point)
   -1                  -3;                 % class 3 (blow-up / unbounded)
];

nPairs  = size(paramPairs, 1);    % should be 4
nIC     = 1;                      % number of random initial conditions per (j1,j2)

%% -------- Per-pair integration times --------
Tfinal_vec = [
    1e6;      % for pair 1
    1e7;      % for pair 2
    1e2;      % for pair 3
    1e6;      % for pair 4
];

%% -------- Per-pair plotting windows (different per subplot) --------
% Left column: Layla's love vs time
timeWindows_L = [
    9e5     1e6;     % pair 1, time-series window
    9.99999e6     1e7;     % pair 2
    0       15;      % pair 3
    9.9995e5     1e6;     % pair 4
];

% Right column: trajectories in (L_r, L_i)
% (Here you can choose *different* windows if you like.)
timeWindows_traj = [
    1e5   1e6;     % pair 1, phase-portrait window
    9.99999e6   1e7;     % pair 2
    0      20;     % pair 3
    9.99e5     1e6;     % pair 4
];

%% -------- Base parameters (fixed part of 'param') --------
param_base = [-1; 1; 1; -1; -1.8; -1.8; 0; 0];

%% -------- Initial conditions --------
rng(1);                            % for reproducibility
x0 = -1 + 2.*rand(nIC, 8);         % random ICs in [-1,1]^8

%% -------- Figure & layout (4 rows x 2 columns) --------
figure(1); clf
set(gcf, 'Color', 'w', 'Position', [100 100 900 1000]);

set(groot, 'DefaultAxesFontName', 'Times New Roman');
set(groot, 'DefaultAxesFontSize', 10);
set(groot, 'DefaultLineLineWidth', 0.8);

tiledlayout(nPairs, 2, 'TileSpacing', 'compact', 'Padding', 'compact');

for idx = 1:nPairs
    % ----- current parameter pair -----
    I = paramPairs(idx, 1);    % j1
    J = paramPairs(idx, 2);    % j2

    param = param_base;
    param(7) = I;
    param(8) = J;

    % ----- current integration time -----
    Tfinal = Tfinal_vec(idx);
    t_span = [0, Tfinal];

    % ----- current plotting windows -----
    tL_min = timeWindows_L(idx, 1);
    tL_max = timeWindows_L(idx, 2);

    tT_min = timeWindows_traj(idx, 1);
    tT_max = timeWindows_traj(idx, 2);

    % ----- Left subplot: Layla's love vs time -----
    ax1 = nexttile;
    hold(ax1, 'on');

    % ----- Right subplot: trajectory in (L_r, L_i) -----
    ax2 = nexttile;
    hold(ax2, 'on');

    for k = 1:nIC
        [t, x] = ode45(@(t,x) three_lovers(t, x, param), t_span, x0(k,:));

        % Pair & subplot-specific masks
        mask_L   = (t >= tL_min) & (t <= tL_max);   % for time-series
        mask_trj = (t >= tT_min) & (t <= tT_max);   % for trajectories

        % Layla's love (real part) vs time
        plot(ax1, t(mask_L), x(mask_L, 1));

        % Trajectory in (L_r, L_i), using *different* mask if desired
        plot(ax2, x(mask_trj, 1), x(mask_trj, 2));
    end

    % ----- Style left axis (Layla's love) -----
    ylabel(ax1, 'L_r', 'Interpreter', 'tex');
    if idx == nPairs
        xlabel(ax1, 'time');
    else
        set(ax1, 'XTickLabel', []);
    end
    grid(ax1, 'on'); box(ax1, 'on');

    % ----- Style right axis (trajectory) -----
    ylabel(ax2, 'L_i', 'Interpreter', 'tex');
    if idx == nPairs
        xlabel(ax2, 'L_r', 'Interpreter', 'tex');
    else
        set(ax2, 'XTickLabel', []);
    end
    grid(ax2, 'on'); box(ax2, 'on');

    set([ax1, ax2], 'LineWidth', 0.8, 'FontSize', 10);
end

%% -------- Global column headings --------
annotation('textbox', [0.15 0.96 0.3 0.03], ...
    'String', 'Layla''s love', ...
    'HorizontalAlignment', 'center', ...
    'VerticalAlignment', 'bottom', ...
    'FontSize', 12, ...
    'FontWeight', 'bold', ...
    'LineStyle', 'none', ...
    'Interpreter', 'tex');

annotation('textbox', [0.55 0.96 0.3 0.03], ...
    'String', 'Trajectories', ...
    'HorizontalAlignment', 'center', ...
    'VerticalAlignment', 'bottom', ...
    'FontSize', 12, ...
    'FontWeight', 'bold', ...
    'LineStyle', 'none', ...
    'Interpreter', 'tex');

%% -------- High-resolution export --------
exportgraphics(gcf, 'Plots/three_lovers_4params_pub.png', 'Resolution', 300);
