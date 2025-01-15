function drawHist(fun_histories, maxcv_histories, merit_histories, fun_init, maxcv_init, merit_init, solver_names, cell_axs_summary, is_cum, p_type, problem_n, n_eval, profile_options, default_height)
%DRAWHIST draws the history plots of the function values, the maximum constraint violation, and the merit function values.

    fun_histories = processHistYaxes(fun_histories, fun_init);
    maxcv_histories = processHistYaxes(maxcv_histories, maxcv_init);
    merit_histories = processHistYaxes(merit_histories, merit_init);
    
    % Define the shift of the y-axis. Shift the y-axis if there is value that is too close to zero.
    y_shift_fun = computeYShift(fun_histories, profile_options);

    % First, draw the histories of function values.
    drawFunMaxcvMeritHist(cell_axs_summary{1}, fun_histories, solver_names, is_cum, problem_n, y_shift_fun, n_eval, profile_options);
    [~, formatted_fun_shift] = formatFloatScientificLatex(y_shift_fun, 3);

    maxlength_fun = length(['Cummin of function values shifted above by $', formatted_fun_shift, '$']);
    label_fontsize_fun = min(12, 1.5 * default_height / maxlength_fun);

    if is_cum
        if y_shift_fun > 0
            y_label = ['Cummin of function values shifted above by $', formatted_fun_shift, '$'];
            ylabel(cell_axs_summary{1}, y_label, 'Interpreter', 'latex', 'FontSize', label_fontsize_fun);
        else
            ylabel(cell_axs_summary{1}, 'Cummin of function values', 'Interpreter', 'latex', 'FontSize', label_fontsize_fun);
        end
    else
        if y_shift_fun > 0
            y_label = ['Function values shifted above by $', formatted_fun_shift, '$'];
            ylabel(cell_axs_summary{1}, y_label, 'Interpreter', 'latex', 'FontSize', label_fontsize_fun);
        else
            ylabel(cell_axs_summary{1}, 'Function values', 'Interpreter', 'latex', 'FontSize', label_fontsize_fun);
        end
    end

    % If the problem is unconstrained, do not draw the histories of maximum constraint violations and merit function values.
    if strcmp(p_type, 'unconstrained')
        return;
    end

    % Do the same for the maximum constraint violations and the merit function values.
    y_shift_maxcv = computeYShift(maxcv_histories, profile_options);
    y_shift_merit = computeYShift(merit_histories, profile_options);

    % Second, draw the histories of maximum constraint violations and merit function values.
    drawFunMaxcvMeritHist(cell_axs_summary{2}, maxcv_histories, solver_names, is_cum, problem_n, y_shift_maxcv, n_eval, profile_options);
    [~, formatted_maxcv_shift] = formatFloatScientificLatex(y_shift_maxcv, 3);
    drawFunMaxcvMeritHist(cell_axs_summary{3}, merit_histories, solver_names, is_cum, problem_n, y_shift_merit, n_eval, profile_options);
    [~, formatted_merit_shift] = formatFloatScientificLatex(y_shift_merit, 3);

    maxlength_maxcv = length(['Cummin of maximum constraint violations shifted above by $', formatted_maxcv_shift, '$']);
    label_fontsize_maxcv = min(12, 1.5 * default_height / maxlength_maxcv);
    maxlength_merit = length(['Cummin of merit function values shifted above by $', formatted_merit_shift, '$']);
    label_fontsize_merit = min(12, 1.5 * default_height / maxlength_merit);
    if is_cum
        if y_shift_maxcv > 0
            y_label = ['Cummin of maximum constraint violations shifted above by $', formatted_maxcv_shift, '$'];
            ylabel(cell_axs_summary{2}, y_label, 'Interpreter', 'latex', 'FontSize', label_fontsize_maxcv);
        else
            ylabel(cell_axs_summary{2}, 'Cummin of maximum constraint violations', 'Interpreter', 'latex', 'FontSize', label_fontsize_maxcv);
        end
        if y_shift_merit > 0
            y_label = ['Cummin of merit function values shifted above by $', formatted_merit_shift, '$'];
            ylabel(cell_axs_summary{3}, y_label, 'Interpreter', 'latex', 'FontSize', label_fontsize_merit);
        else
            ylabel(cell_axs_summary{3}, 'Cummin of merit function values', 'Interpreter', 'latex', 'FontSize', label_fontsize_merit);
        end
    else
        if y_shift_maxcv > 0
            y_label = ['Maximum constraint violations shifted above by $', formatted_maxcv_shift, '$'];
            ylabel(cell_axs_summary{2}, y_label, 'Interpreter', 'latex', 'FontSize', label_fontsize_maxcv);
        else
            ylabel(cell_axs_summary{2}, 'Maximum constraint violations', 'Interpreter', 'latex', 'FontSize', label_fontsize_maxcv);
        end
        if y_shift_merit > 0
            y_label = ['Merit function values shifted above by $', formatted_merit_shift, '$'];
            ylabel(cell_axs_summary{3}, y_label, 'Interpreter', 'latex', 'FontSize', label_fontsize_merit);
        else
            ylabel(cell_axs_summary{3}, 'Merit function values', 'Interpreter', 'latex', 'FontSize', label_fontsize_merit);
        end
    end

end

function y_shift = computeYShift(history, profile_options)

    y_shift = 0;
    if strcmp(profile_options.(ProfileOptionKey.RANGE_TYPE.value), 'meanstd')
        y_mean = squeeze(mean(history, 2));
        y_std = squeeze(std(history, 0, 2));
        y_lower = y_mean - y_std;
        y_min = min(y_lower(:));
    else
        y_min = min(history(:));
    end
    
    % Shift the y-axis if there is value that is smaller than 1e-12 and the values are not all the same.
    if any(diff(history(:))) && y_min <= 1e-12
        y_shift = 1e-12 - y_min;
    end
end