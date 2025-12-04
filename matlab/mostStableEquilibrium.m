function [u_star, sol_star, branchIdx, lambda_star, ...
          sols_all, u_all, real_eig_all] = mostStableEquilibrium(a, c, get_all_real)
    if nargin < 3
        get_all_real = false;
    end

    % --- default outputs in case we error early ---
    u_star      = [];
    sol_star    = [];
    branchIdx   = [];
    lambda_star = [];
    sols_all    = [];
    u_all       = [];
    real_eig_all = [];

    % ----- symbolic solve for equilibria -----
    syms x y z w

    % Equilibrium equations
    % a1=-1, a2=1, a3=1, a4=-1, b1=0, b2=-1, b3=-1)
    eq1 = -1 + y^2 + 0.0*(x - 1i*y);
    eq2 =  1 + x^2 - 1.0*(y - 1i*x) + c*z;
    eq3 = -1 + w^2 + 0.0*(z - 1i*w);
    eq4 =  1 + z^2 - 1.0*(w - 1i*z) + a*x;

    [solx, soly, solz, solw] = solve([eq1, eq2, eq3, eq4], [x, y, z, w]);

    sols_all = vpa([solx, soly, solz, solw]);
    nSol = size(sols_all, 1);

    if nSol == 0
        error('No equilibria found for a = %g, c = %g.', a, c);
    end

    % storage
    u_all = zeros(nSol, 8);   % each row = 8D real equilibrium
    pvals = zeros(nSol, 1);   % max Re(eig) per branch

    if get_all_real
        real_eig_all = cell(nSol, 1); % each cell: real parts of eigenvalues
    else
        real_eig_all = [];
    end

    % ----- loop over all fixed points -----
    for k = 1:nSol
        xk = sols_all(k, 1);
        yk = sols_all(k, 2);
        zk = sols_all(k, 3);
        wk = sols_all(k, 4);

        % 8D real representation
        u = [ real(xk), imag(xk), ...
              real(yk), imag(yk), ...
              real(zk), imag(zk), ...
              real(wk), imag(wk) ];

        u_all(k, :) = u;

        % all eigenvalues at this equilibrium 
        params = [-1, 1, 1, -1, 0, -1, -1, c, a];
        eigenvalues = evalsjacLMW(u, params);

        % max real part for this branch
        pvals(k) = max(real(eigenvalues));

        % optionally store real parts of all eigenvalues
        if get_all_real
            real_eig_all{k} = real(eigenvalues);
        end
    end

    % ----- pick the "most stable" branch (smallest max Re \lambda) -----
    [lambda_star, branchIdx] = min(pvals);

    sol_star = sols_all(branchIdx, :);
    u_star   = u_all(branchIdx, :);
end
