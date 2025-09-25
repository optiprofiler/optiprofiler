function exportHist(problem_name, problem_type, problem_dim, solver_names, solvers_success, fun_history, maxcv_history, fun_inits, maxcv_inits, n_eval, profile_options, path_hist_plots)
    if profile_options.(ProfileOptionKey.SCORE_ONLY.value) || all(~solvers_success(:))
        return;
    end

    try
        merit_fun = profile_options.(ProfileOptionKey.MERIT_FUN.value);
        try
            merit_history = meritFunCompute(merit_fun, fun_history, maxcv_history, maxcv_inits, 'single');
            merit_inits = meritFunCompute(merit_fun, fun_inits, maxcv_inits, maxcv_inits);
        catch
            error("MATLAB:solveOneProblem:merit_fun_error", "Error occurred while calculating the merit values. Please check the merit function.");
        end

        % Create the figure for the summary.
        warning('off');
        if strcmp(problem_type, 'u')
            n_cols = 1;
        else
            n_cols = 3;
        end
        defaultFigurePosition = get(0, 'DefaultFigurePosition');
        default_width = defaultFigurePosition(3);
        default_height = defaultFigurePosition(4);
        fig_summary = figure('Position', [defaultFigurePosition(1:2), n_cols * default_width, 2 * default_height], 'visible', 'off');
        T_summary = tiledlayout(fig_summary, 2, 1, 'Padding', 'compact', 'TileSpacing', 'compact');
        F_title = strrep(profile_options.feature_stamp, '_', '\_');
        P_title = strrep(problem_name, '_', '\_');
        T_title = ['Solving ``', P_title, '" with ``', F_title, '" feature'];
        title_fontsize = min(12, 1.2 * default_width / length(T_title));
        title(T_summary, ['Solving ``', P_title, '" with ``', F_title, '" feature'], 'Interpreter', 'latex', 'FontSize', title_fontsize);
        % Use gobjects to create arrays of handles and axes.
        t_summary = gobjects(2, 1);
        axs_summary = gobjects([2, 1, 1, n_cols]);
        i_axs = 0;
        for i = 1:2
            t_summary(i) = tiledlayout(T_summary, 1, n_cols, 'Padding', 'compact', 'TileSpacing', 'compact');
            t_summary(i).Layout.Tile = i;
            for j = 1:n_cols
                i_axs = i_axs + 1;
                axs_summary(i_axs) = nexttile(t_summary(i));
            end
        end
        ylabel(t_summary(1), "History profiles", 'Interpreter', 'latex', 'FontSize', 14);
        ylabel(t_summary(2), "Cummin history profiles", 'Interpreter', 'latex', 'FontSize', 14);

        if strcmp(problem_type, 'u')
            cell_axs_summary = {axs_summary(1)};
            cell_axs_summary_cum = {axs_summary(2)};
        else
            cell_axs_summary = {axs_summary(1), axs_summary(2), axs_summary(3)};
            cell_axs_summary_cum = {axs_summary(4), axs_summary(5), axs_summary(6)};
        end

        pdf_hist_file_name = [regexprep(regexprep(regexprep(strrep(problem_name,' ','_'),'[^a-zA-Z0-9\-_]',''),'[-_]+','_'),'^[-_]+',''), '.pdf'];
        pdf_summary = fullfile(path_hist_plots, pdf_hist_file_name);
        processed_solver_names = cellfun(@(s) strrep(s, '_', '\_'), solver_names, 'UniformOutput', false);

        drawHist(fun_history, maxcv_history, merit_history, fun_inits, maxcv_inits, merit_inits, processed_solver_names, cell_axs_summary, false, problem_type, problem_dim, n_eval, profile_options, default_height);
        drawHist(fun_history, maxcv_history, merit_history, fun_inits, maxcv_inits, merit_inits, processed_solver_names, cell_axs_summary_cum, true, problem_type, problem_dim, n_eval, profile_options, default_height);

        exportgraphics(fig_summary, pdf_summary, 'ContentType', 'vector');
        warning('on');
        close(fig_summary);
    catch Exception
        if ~profile_options.(ProfileOptionKey.SILENT.value)
            fprintf("INFO: An error occurred while plotting the history plots of the problem %s: %s\n", problem_name, Exception.message);
        end
    end
end