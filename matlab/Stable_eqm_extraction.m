close all; clear all; clc

j1 = 3.0;
j2 = -1.0;

[u_star, sol_star, branchIdx, lambda_star, sols_all, u_all, real_eig_all] = ...
    mostStableEquilibrium(j2, j1, true);

fprintf('Most stable branch index: %d\n', branchIdx);
fprintf('Most stable max Re(eig):  %g\n', lambda_star);

for k = 1:numel(real_eig_all)
    fprintf('\nBranch %d equilibrium (8D real):\n', k);
    disp(u_all(k,:));
    fprintf('Real parts of eigenvalues at this FP:\n');
    disp(real_eig_all{k}.');   % row
end
