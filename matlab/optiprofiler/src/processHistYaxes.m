function value_histories_processed = processHistYaxes(value_histories, value_init, is_cum)
%PROCESSHISTYAXES Process the value_histories of the y-axis data.

    n_solvers = size(value_histories, 1);
    n_runs = size(value_histories, 2);

    mask_hist_nan_inf = ~isfinite(value_histories);
    value_histories(mask_hist_nan_inf) = value_init;
    value_max = max(value_histories(:));
    value_min = min(value_histories(:));
    if is_cum
        value_histories_processed = cummin(value_histories, 3);
    else
        value_histories_processed = value_histories;
        value_histories_processed(mask_hist_nan_inf) = value_min + 1.5 * (value_max - value_min);
    end

    % Find the last evaluation where the value decreases.
    mask_diff = diff(value_histories, 1, 3) < 0;
    if ~any(mask_diff(:))
        max_lastd = 2;
    else
        for i_solver = 1:n_solvers
            for i_run = 1:n_runs
                lastd_tmp = find(mask_diff(i_solver, i_run, :), 1, 'last');
                if isempty(lastd_tmp)
                    lastd(i_solver, i_run) = 2;
                else
                    lastd(i_solver, i_run) = lastd_tmp;
                end
            end
        end
        max_lastd = max(lastd(:));
    end

    % Truncate the value_histories of the function values according to `max_last_true_fun_eval'.
    value_histories_processed = value_histories_processed(:, :, 1:max_lastd);

end