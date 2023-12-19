function y = y_profile(sort_x, n_problems, n_runs)
    n_solvers = size(sort_x, 2);
    y = zeros(n_problems * n_runs, n_solvers);
    for i_run = 1:n_runs
        for i_solver = 1:n_solvers
            y_partial = NaN(n_problems * n_runs, 1);
            y_partial((i_run - 1) * n_problems + 1:i_run * n_problems) = linspace(1 / n_problems, 1.0, n_problems)';
            y_partial = y_partial(sort_x(:, i_solver));
            for i_problem = 1:n_problems * n_runs
                if isnan(y_partial(i_problem))
                    if i_problem > 1
                        y_partial(i_problem) = y_partial(i_problem - 1);
                    else
                        y_partial(i_problem) = 0.0;
                    end
                end
            end
            y(:, i_solver) = y(:, i_solver) + y_partial;
        end
    end
    y = y / n_runs;
end