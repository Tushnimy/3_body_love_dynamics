clear; close all; clc

%% -------- Parameter pairs (j1, j2) to test --------
paramPairs = [
    1.25   3;
    -1.0   3;
    3.00   3
];

nPairs  = size(paramPairs, 1);

%% -------- Per-pair integration times --------
Tfinal_vec = [
    1e5;
    1e5;
    0.7e2
];

%% -------- Per-pair plotting thresholds --------
T_L_start_vec = [
    9.99e4;
    9.99e4;
    0.00
];

T_traj_start_vec = [
    9e4;
    0.00;
    0.00
];

%% -------- Base parameters --------
param_base = [-1; 1; 1; -1; 0; -1; -1; 0; 0];

%% -------- Initial conditions --------
IC_list = cell(nPairs, 1);
IC_list{1} = -3 + 6.*rand(4, 8);
IC_list{2} = -1 + 2.*rand(1,8) + 1e-15*(-1 + 2.*rand(4, 8));
IC_list{3} =  [-3,-1,1,0,-3,-1,1,0] + 1e-1*(-1 + 2.*rand(4,8));

%% -------- Figure setup --------
figure(1); clf
set(gcf, 'Color', 'w', 'Position', [100 100 900 1000]);

tiledlayout(nPairs, 2, 'TileSpacing', 'compact', 'Padding', 'compact');

%% ======== MAIN LOOP ========
for idx = 1:nPairs

    % ----- current parameters -----
    j1 = paramPairs(idx, 1);
    j2 = paramPairs(idx, 2);

    param = param_base;
    param(8) = j1;
    param(9) = j2;

    % ----- time and thresholds -----
    Tfinal       = Tfinal_vec(idx);
    t_span       = [0, Tfinal];
    T_L_start    = T_L_start_vec(idx);
    T_traj_start = T_traj_start_vec(idx);

    % ----- ICs -----
    x0_mat = IC_list{idx};
    nIC_idx = size(x0_mat, 1);

    %% ----- Left subplot: L_r vs time -----
    ax1 = nexttile;
    hold(ax1, 'on');

    %% ----- Right subplot: trajectory in (L_r, L_i) -----
    ax2 = nexttile;
    hold(ax2, 'on');

    for k = 1:nIC_idx
        x0 = x0_mat(k, :);   % 1×8 IC

        [t, x] = ode45(@(t,x) three_lovers(t, x, param), t_span, x0);

        % masks
        mask_L   = (t > T_L_start);
        mask_trj = (t > T_traj_start);

        % Left: L_r(t)
        plot(ax1, t(mask_L), x(mask_L, 1));

        % Right: (L_r, L_i)
        plot(ax2, x(mask_trj, 1), x(mask_trj, 2));
    end

    %% ------ Style left axis ------
    xlabel(ax1, 'time');
    ylabel(ax1, 'L_r');
    grid(ax1, 'on'); box(ax1, 'on');

    %% ------ Style right axis ------
    xlabel(ax2, 'L_r');
    ylabel(ax2, 'L_i');
    grid(ax2, 'on'); box(ax2, 'on');

end

%% -------- Export --------
exportgraphics(gcf, 'Plots/three_lovers_4params_pub.png', 'Resolution', 300);
