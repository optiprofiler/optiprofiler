function value_histories_processed = processHistYaxes(value_histories, value_inits)
%PROCESSHISTYAXES Process the value_histories of the y-axis data.
%   value_histories has the size of (n_solvers, n_runs, n_evals).
%   value_inits is a vector of size (n_runs, 1).
%   This function processes the value_histories by replacing NaN and Inf values
%   with a value larger than the maximum finite value in value_histories.
%

    n_runs = size(value_histories, 2);
    value_histories_processed = value_histories;
    for i_run = 1:n_runs
        mask_hist_nan_inf = ~isfinite(value_histories(:, i_run, :));
        value_histories(mask_hist_nan_inf) = value_inits(i_run);
        value_max = max(value_histories(:, i_run, :), [], 'all', 'omitnan');
        value_min = min(value_histories(:, i_run, :), [], 'all', 'omitnan');
        value_histories_processed(mask_hist_nan_inf) = value_min + 1.5 * (value_max - value_min);
    end

    % Following code is used to truncate the value_histories according to the last evaluation where the value decreases. But at this moment, we do not use this strategy.

    % n_solvers = size(value_histories, 1);
    % n_runs = size(value_histories, 2);

    % % Find the last evaluation where the value decreases.
    % mask_diff = diff(value_histories, 1, 3) < 0;
    % if ~any(mask_diff(:))
    %     max_lastd = 2;
    % else
    %     for i_solver = 1:n_solvers
    %         for i_run = 1:n_runs
    %             lastd_tmp = find(mask_diff(i_solver, i_run, :), 1, 'last');
    %             if isempty(lastd_tmp)
    %                 lastd(i_solver, i_run) = 2;
    %             else
    %                 lastd(i_solver, i_run) = lastd_tmp;
    %             end
    %         end
    %     end
    %     max_lastd = max(lastd(:));
    % end

    % % Truncate the value_histories of the function values according to `max_last_true_fun_eval'.
    % value_histories_processed = value_histories_processed(:, :, 1:max_lastd);

end