%% Clean up environment
close all; clc; clear;

%% Parameters
t = 50;
r = linspace(-3, 3, t);
p = zeros(16, t);
syms x y z w
a = 2.5;
sols = [];

%% Main loop: solve equations and compute eigenvalues
for i = 1:length(r)
    c = r(i);

    % System of equations
    eq1 = -1 + y^2 - 1*(x - 1i*y);
    eq2 =  1 + x^2 - 1.8*(y - 1i*x) + c*z;
    eq3 = -1 + w^2 - 1*(z - 1i*w);
    eq4 =  1 + z^2 - 1.8*(w - 1i*z) + a*x;

    % Solve symbolic system
    [solx, soly, solz, solw] = solve([eq1, eq2, eq3, eq4], [x, y, z, w]);
    solutions = [vpa(solx), vpa(soly), vpa(solz), vpa(solw)];

    % Extract real/imag parts
    sols = [real(solutions(:,1)), imag(solutions(:,1)), ...
            real(solutions(:,2)), imag(solutions(:,2)), ...
            real(solutions(:,3)), imag(solutions(:,3)), ...
            real(solutions(:,4)), imag(solutions(:,4))];

    % Evaluate Jacobian and maximum real eigenvalue
    for k = 1:16
        eigenvalues = evalsjacLMW(sols(k,:), c, a);  % user-defined function
        Real = findMaxRealVec(eigenvalues, 0);       % user-defined function
        p(k,i) = Real;
    end
end

%% Postprocess: minimum across 16 branches
mat = min(p, [], 1);

%% ========== PLOTTING (publication quality) ==========
fig = figure('Color','w','Units','inches','Position',[1 1 8 4]);
set(groot, ...
    'defaultTextInterpreter','latex', ...
    'defaultAxesTickLabelInterpreter','latex', ...
    'defaultLegendInterpreter','latex', ...
    'defaultAxesFontName','Times', ...
    'defaultAxesFontSize',11);

hold on

% 1) Plot the curve AS-IS (single neutral color)
hCurve = plot(r, mat, 'LineWidth', 1.8, 'Color', [0.1 0.1 0.1]);

% 2) Color the x-axis (y=0) by the sign of mat
green = [0.20 0.60 0.20];   % one stable fixed point (below x-axis)
blue  = [0.10 0.40 0.90];   % stable LC / multistable periodic/aperiodic (above x-axis)

isStable = mat <= 0;  % true where curve is below or on x-axis
edges = [true, diff(isStable)~=0, true];
breaks = find(edges); % indices that bound constant-sign segments

for s = 1:numel(breaks)-1
    i1 = breaks(s);
    i2 = breaks(s+1)-1;
    if i1 <= i2
        rseg = r(i1:i2);
        yseg = zeros(size(rseg));
        col  = green;
        if ~isStable(i1)
            col = blue;
        end
        plot(rseg, yseg, '-', 'Color', col, 'LineWidth', 3, 'Clipping', 'on');
    end
end

% 3) Draw y-axis only (no default y=0 line; we just colored it)
xline(0, '-', 'Color', [0 0 0], 'LineWidth', 0.8, 'Alpha', 0.5);

% 4) Labels and title
xlabel('$j_1$', 'FontSize', 12);
ylabel('$\max\!\big(\mathrm{Re}\,\lambda\big)$', 'FontSize', 12);
title('Eigenvalue spectrum vs.\ $j_1$', 'FontSize', 13, 'FontWeight', 'normal');

% 5) Axes style
grid on;
set(gca, 'GridAlpha', 0.2, ...
         'LineWidth', 0.75, ...
         'TickDir', 'out', ...
         'Box', 'off');
xlim([min(r) max(r)]);
ylim('padded');

% 6) Legend: only color-condition entries (NOT the curve)
hLegStable   = plot(nan, nan, '-', 'Color', green, 'LineWidth', 3);
hLegUnstable = plot(nan, nan, '-', 'Color', blue,  'LineWidth', 3);
legend([hLegStable, hLegUnstable], ...
       {'one stable F.P.', ...
        'stable L.C. or multistable periodic/aperiodic'}, ...
       'Location','northwest', 'Box','off');

% 7) Export
exportgraphics(fig, 'eigenvalue_plot.pdf', 'ContentType','vector');
exportgraphics(fig, 'eigenvalue_plot.png', 'Resolution', 300);

disp('Plot saved as eigenvalue_plot.pdf and eigenvalue_plot.png');
