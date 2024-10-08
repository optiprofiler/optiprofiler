function drawHist(fun_histories, maxcv_histories, merit_histories, fun_init, maxcv_init, merit_init, labels, cell_axs_summary, is_cum, profile_options, problem_type, problem_n)
%DRAWHIST draws the history plots of the function values, the maximum constraint violation, and the merit function values.

    fun_histories = processHistYaxes(fun_histories, fun_init, is_cum);
    maxcv_histories = processHistYaxes(maxcv_histories, maxcv_init, is_cum);
    merit_histories = processHistYaxes(merit_histories, merit_init, is_cum);

    if strcmp(problem_type, 'unconstrained')
        drawFunMaxcvMeritHist(cell_axs_summary{1}, fun_histories, labels, false, profile_options, problem_n);
        if is_cum
            ylabel(cell_axs_summary{1}, 'Cummin of function values', 'Interpreter', 'latex');
        else
            ylabel(cell_axs_summary{1}, 'Function values', 'Interpreter', 'latex');
        end
    else
        drawFunMaxcvMeritHist(cell_axs_summary{1}, fun_histories, labels, false, profile_options, problem_n);
        drawFunMaxcvMeritHist(cell_axs_summary{2}, maxcv_histories, labels, true, profile_options, problem_n);
        drawFunMaxcvMeritHist(cell_axs_summary{3}, merit_histories, labels, false, profile_options, problem_n);
        if is_cum
            ylabel(cell_axs_summary{1}, 'Cummin of function values', 'Interpreter', 'latex');
            ylabel(cell_axs_summary{2}, 'Cummin of maximum constraint violations', 'Interpreter', 'latex');
            ylabel(cell_axs_summary{3}, 'Cummin of merit function values', 'Interpreter', 'latex');
        else
            ylabel(cell_axs_summary{1}, 'Function values', 'Interpreter', 'latex');
            ylabel(cell_axs_summary{2}, 'Maximum constraint violations', 'Interpreter', 'latex');
            ylabel(cell_axs_summary{3}, 'Merit function values', 'Interpreter', 'latex');
        end
    end
    
end