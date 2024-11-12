function drawHist(fun_histories, maxcv_histories, merit_histories, fun_init, maxcv_init, merit_init, labels, cell_axs_summary, is_cum, profile_options, problem_type, problem_n)
%DRAWHIST draws the history plots of the function values, the maximum constraint violation, and the merit function values.

    fun_histories = processHistYaxes(fun_histories, fun_init, is_cum);
    maxcv_histories = processHistYaxes(maxcv_histories, maxcv_init, is_cum);
    merit_histories = processHistYaxes(merit_histories, merit_init, is_cum);

    % First, draw the histories of function values.
    y_shift_fhist = drawFunMaxcvMeritHist(cell_axs_summary{1}, fun_histories, labels, false, profile_options, problem_n);
    [~, formatted_fhist] = formatFloatScientificLatex(y_shift_fhist, 3);

    if is_cum
        if y_shift_fhist > 0
            ylabel(cell_axs_summary{1}, ['Cummin of function values shifted above by $', formatted_fhist, '$'], 'Interpreter', 'latex');
        else
            ylabel(cell_axs_summary{1}, 'Cummin of function values', 'Interpreter', 'latex');
        end
    else
        if y_shift_fhist > 0
            ylabel(cell_axs_summary{1}, ['Function values shifted above by $', formatted_fhist, '$'], 'Interpreter', 'latex');
        else
            ylabel(cell_axs_summary{1}, 'Function values', 'Interpreter', 'latex');
        end
    end

    % If the problem is unconstrained, do not draw the histories of maximum constraint violations and merit function values.
    if strcmp(problem_type, 'unconstrained')
        return;
    end

    % Second, draw the histories of maximum constraint violations and merit function values. Note that there is no need to shift the y-axis since it is always nonnegative.

    drawFunMaxcvMeritHist(cell_axs_summary{2}, maxcv_histories, labels, true, profile_options, problem_n);
    y_shift_merithist = drawFunMaxcvMeritHist(cell_axs_summary{3}, merit_histories, labels, false, profile_options, problem_n);
    [~, formatted_merithist] = formatFloatScientificLatex(y_shift_merithist, 3);
    if is_cum
        ylabel(cell_axs_summary{2}, 'Cummin of maximum constraint violations', 'Interpreter', 'latex');
        if y_shift_merithist > 0
            ylabel(cell_axs_summary{3}, 'Cummin of merit function values shifted above by $' + formatted_merithist + '$', 'Interpreter', 'latex');
        else
            ylabel(cell_axs_summary{3}, 'Cummin of merit function values', 'Interpreter', 'latex');
        end
    else
        ylabel(cell_axs_summary{2}, 'Maximum constraint violations', 'Interpreter', 'latex');
        if y_shift_merithist > 0
            ylabel(cell_axs_summary{3}, 'Merit function values shifted above by $' + formatted_merithist + '$', 'Interpreter', 'latex');
        else
            ylabel(cell_axs_summary{3}, 'Merit function values', 'Interpreter', 'latex');
        end
    end

    
end