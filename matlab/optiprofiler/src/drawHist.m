function [fig_fun, fig_maxcv, fig_merit] = drawHist(fun_histories, maxcv_histories, merit_histories, fun_init, maxcv_init, merit_init, labels, cell_axs_summary, is_summary, is_fun, is_maxcv, is_merit, construct_cum, profile_options)
%DRAWHIST draws the profiles of the histories of the function values, the maximum constraint violation, and the merit function values when the whole test suite only contains one problem.

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Create the individual figures.
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    fig_fun = figure('visible', 'off');
    t_fun = tiledlayout(fig_fun, 1, 1, 'Padding', 'compact', 'TileSpacing', 'compact');
    ax_fun = nexttile(t_fun);
    fig_maxcv = figure('visible', 'off');
    t_maxcv = tiledlayout(fig_maxcv, 1, 1, 'Padding', 'compact', 'TileSpacing', 'compact');
    ax_maxcv = nexttile(t_maxcv);
    fig_merit = figure('visible', 'off');
    t_merit = tiledlayout(fig_merit, 1, 1, 'Padding', 'compact', 'TileSpacing', 'compact');
    ax_merit = nexttile(t_merit);
    % Draw the figures of the histories of the function values, the maximum constraint violation, and the merit function values.
    fun_histories = processHistYaxes(fun_histories, fun_init, construct_cum);
    maxcv_histories = processHistYaxes(maxcv_histories, maxcv_init, construct_cum);
    merit_histories = processHistYaxes(merit_histories, merit_init, construct_cum);
    drawFunMaxcvMeritHist(ax_fun, fun_histories, labels, false, profile_options);
    drawFunMaxcvMeritHist(ax_maxcv, maxcv_histories, labels, true, profile_options);
    drawFunMaxcvMeritHist(ax_merit, merit_histories, labels, false, profile_options);
    % Set y-axis labels.
    if construct_cum
        ylabel(ax_fun, 'Cummin of function values', 'Interpreter', 'latex');
        ylabel(ax_maxcv, 'Cummin of maximum constraint violations', 'Interpreter', 'latex');
        ylabel(ax_merit, 'Cummin of merit function values', 'Interpreter', 'latex');
    else
        ylabel(ax_fun, 'Function values', 'Interpreter', 'latex');
        ylabel(ax_maxcv, 'Maximum constraint violations', 'Interpreter', 'latex');
        ylabel(ax_merit, 'Merit function values', 'Interpreter', 'latex');
    end

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Create the figures in summary.
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    if is_summary
        if is_fun && is_maxcv && is_merit
            drawFunMaxcvMeritHist(cell_axs_summary{1}, fun_histories, labels, false, profile_options);
            drawFunMaxcvMeritHist(cell_axs_summary{2}, maxcv_histories, labels, true, profile_options);
            drawFunMaxcvMeritHist(cell_axs_summary{3}, merit_histories, labels, false, profile_options);
            if construct_cum
                ylabel(cell_axs_summary{1}, 'Cummin of function values', 'Interpreter', 'latex');
                ylabel(cell_axs_summary{2}, 'Cummin of maximum constraint violations', 'Interpreter', 'latex');
                ylabel(cell_axs_summary{3}, 'Cummin of merit function values', 'Interpreter', 'latex');
            else
                ylabel(cell_axs_summary{1}, 'Function values', 'Interpreter', 'latex');
                ylabel(cell_axs_summary{2}, 'Maximum constraint violations', 'Interpreter', 'latex');
                ylabel(cell_axs_summary{3}, 'Merit function values', 'Interpreter', 'latex');
            end
        elseif is_fun && is_maxcv
            drawFunMaxcvMeritHist(cell_axs_summary{1}, fun_histories, labels, false, profile_options);
            drawFunMaxcvMeritHist(cell_axs_summary{2}, maxcv_histories, labels, true, profile_options);
            if construct_cum
                ylabel(cell_axs_summary{1}, 'Cummin of function values', 'Interpreter', 'latex');
                ylabel(cell_axs_summary{2}, 'Cummin of maximum constraint violations', 'Interpreter', 'latex');
            else
                ylabel(cell_axs_summary{1}, 'Function values', 'Interpreter', 'latex');
                ylabel(cell_axs_summary{2}, 'Maximum constraint violations', 'Interpreter', 'latex');
            end
        elseif is_fun && is_merit
            drawFunMaxcvMeritHist(cell_axs_summary{1}, fun_histories, labels, false, profile_options);
            drawFunMaxcvMeritHist(cell_axs_summary{2}, merit_histories, labels, false, profile_options);
            if construct_cum
                ylabel(cell_axs_summary{1}, 'Cummin of function values', 'Interpreter', 'latex');
                ylabel(cell_axs_summary{2}, 'Cummin of merit function values', 'Interpreter', 'latex');
            else
                ylabel(cell_axs_summary{1}, 'Function values', 'Interpreter', 'latex');
                ylabel(cell_axs_summary{2}, 'Merit function values', 'Interpreter', 'latex');
            end
        elseif is_maxcv && is_merit
            drawFunMaxcvMeritHist(cell_axs_summary{1}, maxcv_histories, labels, true, profile_options);
            drawFunMaxcvMeritHist(cell_axs_summary{2}, merit_histories, labels, false, profile_options);
            if construct_cum
                ylabel(cell_axs_summary{1}, 'Cummin of maximum constraint violations', 'Interpreter', 'latex');
                ylabel(cell_axs_summary{2}, 'Cummin of merit function values', 'Interpreter', 'latex');
            else
                ylabel(cell_axs_summary{1}, 'Maximum constraint violations', 'Interpreter', 'latex');
                ylabel(cell_axs_summary{2}, 'Merit function values', 'Interpreter', 'latex');
            end
        elseif is_fun
            drawFunMaxcvMeritHist(cell_axs_summary{1}, fun_histories, labels, false, profile_options);
            if construct_cum
                ylabel(cell_axs_summary{1}, 'Cummin of function values', 'Interpreter', 'latex');
            else
                ylabel(cell_axs_summary{1}, 'Function values', 'Interpreter', 'latex');
            end
        elseif is_maxcv
            drawFunMaxcvMeritHist(cell_axs_summary{1}, maxcv_histories, labels, true, profile_options);
            if construct_cum
                ylabel(cell_axs_summary{1}, 'Cummin of maximum constraint violations', 'Interpreter', 'latex');
            else
                ylabel(cell_axs_summary{1}, 'Maximum constraint violations', 'Interpreter', 'latex');
            end
        elseif is_merit
            drawFunMaxcvMeritHist(cell_axs_summary{1}, merit_histories, labels, false, profile_options);
            if construct_cum
                ylabel(cell_axs_summary{1}, 'Cummin of merit function values', 'Interpreter', 'latex');
            else
                ylabel(cell_axs_summary{1}, 'Merit function values', 'Interpreter', 'latex');
            end
        end

    end
end