function drawHist(fun_histories, maxcv_histories, merit_histories, fun_init, maxcv_init, merit_init, labels, cell_axs_summary, is_cum, problem_type, problem_n, n_eval)
%DRAWHIST draws the history plots of the function values, the maximum constraint violation, and the merit function values.

    fun_histories = processHistYaxes(fun_histories, fun_init);
    
    % Define the shift of the y-axis. Shift the y-axis if there is value that is too close to zero.
    y_shift_fun = 0;
    if min(fun_histories(:)) <= 1e-9 * (max(fun_histories(:)) - min(fun_histories(:)))
        y_shift_fun = max(1e-9 * (max(fun_histories(:)) - min(fun_histories(:))) - min(fun_histories(:)), 1e-12);
    end

    % First, draw the histories of function values.
    drawFunMaxcvMeritHist(cell_axs_summary{1}, fun_histories, labels, is_cum, problem_n, y_shift_fun, n_eval);
    [~, formatted_fun_shift] = formatFloatScientificLatex(y_shift_fun, 3);

    if is_cum
        if y_shift_fun > 0
            ylabel(cell_axs_summary{1}, ['Cummin of function values shifted above by $', formatted_fun_shift, '$'], 'Interpreter', 'latex');
        else
            ylabel(cell_axs_summary{1}, 'Cummin of function values', 'Interpreter', 'latex');
        end
    else
        if y_shift_fun > 0
            ylabel(cell_axs_summary{1}, ['Function values shifted above by $', formatted_fun_shift, '$'], 'Interpreter', 'latex');
        else
            ylabel(cell_axs_summary{1}, 'Function values', 'Interpreter', 'latex');
        end
    end

    % If the problem is unconstrained, do not draw the histories of maximum constraint violations and merit function values.
    if strcmp(problem_type, 'unconstrained')
        return;
    end

    % Do the same for the maximum constraint violations and the merit function values.
    maxcv_histories = processHistYaxes(maxcv_histories, maxcv_init);
    merit_histories = processHistYaxes(merit_histories, merit_init);
    y_shift_maxcv = 0;
    y_shift_merit = 0;
    if min(maxcv_histories(:)) <= 1e-9 * (max(maxcv_histories(:)) - min(maxcv_histories(:))) && any(diff(maxcv_histories(:)))
        % Note that if all maxcv_histories are the same (actually we care the case where it is all 0),
        % we do not need to shift the y-axis.
        y_shift_maxcv = max(1e-9 * (max(maxcv_histories(:)) - min(maxcv_histories(:))) - min(maxcv_histories(:)), 1e-12);
    end
    if min(merit_histories(:)) <= 1e-9 * (max(merit_histories(:)) - min(merit_histories(:)))
        y_shift_merit = max(1e-9 * (max(merit_histories(:)) - min(merit_histories(:))) - min(merit_histories(:)), 1e-12);
    end

    % Second, draw the histories of maximum constraint violations and merit function values.

    drawFunMaxcvMeritHist(cell_axs_summary{2}, maxcv_histories, labels, is_cum, problem_n, y_shift_maxcv, n_eval);
    [~, formatted_maxcv_shift] = formatFloatScientificLatex(y_shift_maxcv, 3);
    drawFunMaxcvMeritHist(cell_axs_summary{3}, merit_histories, labels, is_cum, problem_n, y_shift_merit, n_eval);
    [~, formatted_merit_shift] = formatFloatScientificLatex(y_shift_merit, 3);
    if is_cum
        if y_shift_maxcv > 0
            ylabel(cell_axs_summary{2}, 'Cummin of maximum constraint violations shifted above by $' + formatted_maxcv_shift + '$', 'Interpreter', 'latex');
        else
            ylabel(cell_axs_summary{2}, 'Cummin of maximum constraint violations', 'Interpreter', 'latex');
        end
        if y_shift_merit > 0
            ylabel(cell_axs_summary{3}, 'Cummin of merit function values shifted above by $' + formatted_merit_shift + '$', 'Interpreter', 'latex');
        else
            ylabel(cell_axs_summary{3}, 'Cummin of merit function values', 'Interpreter', 'latex');
        end
    else
        if y_shift_maxcv > 0
            ylabel(cell_axs_summary{2}, 'Maximum constraint violations shifted above by $' + formatted_maxcv_shift + '$', 'Interpreter', 'latex');
        else
            ylabel(cell_axs_summary{2}, 'Maximum constraint violations', 'Interpreter', 'latex');
        end
        if y_shift_merit > 0
            ylabel(cell_axs_summary{3}, 'Merit function values shifted above by $' + formatted_merit_shift + '$', 'Interpreter', 'latex');
        else
            ylabel(cell_axs_summary{3}, 'Merit function values', 'Interpreter', 'latex');
        end
    end

    
end