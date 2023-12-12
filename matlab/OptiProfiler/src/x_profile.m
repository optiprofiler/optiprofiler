function [sorted_x, ratio_max, sort_x] = x_profile(work, denominator, ratio_max_init)
    [n_problems, n_solvers, n_runs] = size(work);
    x = NaN(n_solvers, n_problems, n_runs);
    
    for i_run = 1:n_runs
        for i_problem = 1:n_problems
            if ~all(isnan(work(i_problem, :, i_run)))
                x(:, i_problem, i_run) = work(i_problem, :, i_run) / denominator(i_problem, i_run);
            end
        end
    end
    
    ratio_max = max(x(:), [], 'omitnan');
    if isempty(ratio_max) || isnan(ratio_max)
        ratio_max = ratio_max_init;
    end
    ratio_max = max(ratio_max, ratio_max_init);
    
    x(isnan(x)) = 2.0 * ratio_max;
    
    % for i_run = 1:n_runs
    %     x(i_run, :, :) = sort(x(i_run, :, :), 2);
    % end
    x = sort(x, 2);
    
    x = reshape(x, [n_solvers, n_problems * n_runs]);
    x = x';
    
    [sorted_x, sort_x] = sort(x);
    
end