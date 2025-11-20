function fps = find_fixed_points(params)
    options = optimoptions('fsolve', 'Display', 'off');
    num_vars = 8; % Number of variables in the system
    fps = [];

    % Try different initial guesses if the system has multiple fixed points
    for i = 1:10
        initial_guess = rand(num_vars, 1) * 4 - 2; % Random initial guess between -2 and 2
        [fp, ~, exitflag, ~] = fsolve(@(x) LMWtest2(0, x, params), initial_guess, options);
        if exitflag > 0 % Solution converged
            fps = [fps; fp']; % Transpose to make it a row vector for each solution
        end
    end
    
    % Remove duplicates within a tolerance
    fps = unique(fps, 'rows', 'stable');
end
