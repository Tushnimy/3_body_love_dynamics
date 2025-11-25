clear; close all; clc

%% -------- Parameter pairs (j1, j2) to test --------
% Each row is [j1, j2].
% Replace these with whatever 4 pairs you want to showcase.
paramPairs = [
    3  -2;   % example: class 0 (aperiodic)
    -1.5  -2.8;   % example: class 1 (periodic)
    0.1  0.1;   % example: class 2 (fixed point)
     -2  -2;   % example: class 3 (blow-up / unbounded)
];

nPairs  = size(paramPairs, 1);    % should be 4
nIC     = 1;                      % number of random initial conditions per (j1,j2)
Tfinal  = 2e5;                    % final integration time
Twindow = 1e2;                    % only plot the tail in time (remove transient)

%% -------- Base parameters (fixed part of 'param') --------
% Original code: param = [-1;1;1;-1;-1.8;-1.8;I;J];
% Keep first 6 fixed, last two will be set per pair.
param_base = [-1; 1; 1; -1; -1.8; -1.8; 0; 0];

%% -------- Initial conditions --------
rng(1);                            % for reproducibility
x0 = -1 + 2.*rand(nIC, 8);         % 6 random ICs in [-1,1]^8

%% -------- Figure & layout (4 rows x 2 columns) --------
figure(1); clf
set(gcf, 'Color', 'w', 'Position', [100 100 900 1000]);  % big white figure

% Publication-ish defaults
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

    % ----- Left subplot: Layla's love vs time -----
    ax1 = nexttile;
    hold(ax1, 'on');

    % ----- Right subplot: trajectory in (L_r, L_i) -----
    ax2 = nexttile;
    hold(ax2, 'on');

    for k = 1:nIC
        [t, x] = ode45(@(t,x) three_lovers(t, x, param), [0 Tfinal], x0(k,:));

        % Tail of the time series (to avoid transient)
        mask = (t > (Tfinal - Twindow));

        % Layla's love (real part) vs time
        plot(ax1, t(mask), x(mask, 2));

        % Trajectory in (L_r, L_i)
        plot(ax2, x(:,1), x(:,2));
    end

    % ----- Style left axis (Layla's love) -----
    %title(ax1, sprintf('j_1 = %.2f, j_2 = %.2f', I, J), ...
    %      'Interpreter', 'tex', 'FontWeight', 'normal');
    ylabel(ax1, 'L_i', 'Interpreter', 'tex');

    if idx == nPairs
        xlabel(ax1, 'time');
    else
        set(ax1, 'XTickLabel', []);   % hide x-labels for upper rows
    end

    grid(ax1, 'on');
    box(ax1, 'on');

    % ----- Style right axis (trajectory) -----
    %title(ax2, 'Trajectory', 'FontWeight', 'normal');
    ylabel(ax2, 'L_i', 'Interpreter', 'tex');

    if idx == nPairs
        xlabel(ax2, 'L_r', 'Interpreter', 'tex');
    else
        set(ax2, 'XTickLabel', []);
    end

    grid(ax2, 'on');
    box(ax2, 'on');

    set([ax1, ax2], 'LineWidth', 0.8, 'FontSize', 10);
end

%% -------- Global column headings (no tiles consumed) --------
% Normalized figure coordinates: [x y width height]

% Left column: "Layla's love"
annotation('textbox', [0.15 0.96 0.3 0.03], ...
    'String', 'Layla''s love', ...
    'HorizontalAlignment', 'center', ...
    'VerticalAlignment', 'bottom', ...
    'FontSize', 12, ...
    'FontWeight', 'bold', ...
    'LineStyle', 'none', ...
    'Interpreter', 'tex');

% Right column: "Trajectories"
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
