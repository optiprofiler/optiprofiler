function [x, y, ratio_max] = getPerformanceDataProfileAxes(work, denominator)
    % Calculate the axes of the performance and data profiles.

    [n_problems, n_solvers, n_runs] = size(work);

    % Calculate the x-axis values.
    x = NaN(n_solvers, n_problems, n_runs);
    for i_run = 1:n_runs
        for i_problem = 1:n_problems
            x(:, i_problem, i_run) = work(i_problem, :, i_run) / denominator(i_problem, i_run);
        end
    end
    if all(isnan(x(:)))
        ratio_max = eps;
    else
        ratio_max = max(x(:), [], 'omitnan');
    end
    x(isnan(x)) = Inf;
    x = sort(x, 2);
    x = reshape(x, [n_solvers, n_problems * n_runs]);
    x = x';
    [x, sort_index_x] = sort(x);

    % Calculate the y-axis values.
    y = NaN(n_problems * n_runs, n_solvers, n_runs);
    for i_solver = 1:n_solvers
        for i_run = 1:n_runs
            y((i_run - 1) * n_problems + 1:i_run * n_problems, i_solver, i_run) = linspace(1 / n_problems, 1.0, n_problems);
            y_partial = y(:, i_solver, i_run);
            y(:, i_solver, i_run) = y_partial(sort_index_x(:, i_solver));
            for i_problem = 1:n_problems * n_runs
                if isnan(y(i_problem, i_solver, i_run))
                    if i_problem > 1
                        y(i_problem, i_solver, i_run) = y(i_problem - 1, i_solver, i_run);
                    else
                        y(i_problem, i_solver, i_run) = 0;
                    end
                end
            end
        end
    end
    
end