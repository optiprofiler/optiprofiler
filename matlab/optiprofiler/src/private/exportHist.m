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

        warning_state = warning;
        warning('off');
        warning_cleanup = onCleanup(@() warning(warning_state));
        if strcmp(problem_type, 'u')
            n_cols = 1;
        else
            n_cols = 3;
        end
        [~, default_width, default_height] = profileFigurePosition();
        summary_profile_width = default_width + summaryLegendExtraWidth(numel(solver_names), default_width, solver_names);

        pdf_hist_file_name = [regexprep(regexprep(regexprep(strrep(problem_name,' ','_'),'[^a-zA-Z0-9\-_]',''),'[-_]+','_'),'^[-_]+',''), '.pdf'];
        pdf_summary = fullfile(path_hist_plots, pdf_hist_file_name);
        processed_solver_names = cellfun(@escapeLatexText, solver_names, 'UniformOutput', false);

        raw_path = fullfile(path_hist_plots, 'raw');
        cummin_path = fullfile(path_hist_plots, 'cummin');
        if ~exist(raw_path, 'dir')
            mkdir(raw_path);
        end
        if ~exist(cummin_path, 'dir')
            mkdir(cummin_path);
        end

        exportHistoryFigure(pdf_summary, 'combined', problem_name, problem_type, problem_dim, processed_solver_names, ...
            fun_history, maxcv_history, merit_history, fun_inits, maxcv_inits, merit_inits, n_eval, ...
            profile_options, n_cols, default_width, default_height, summary_profile_width);
        exportHistoryFigure(fullfile(raw_path, pdf_hist_file_name), 'raw', problem_name, problem_type, problem_dim, processed_solver_names, ...
            fun_history, maxcv_history, merit_history, fun_inits, maxcv_inits, merit_inits, n_eval, ...
            profile_options, n_cols, default_width, default_height, summary_profile_width);
        exportHistoryFigure(fullfile(cummin_path, pdf_hist_file_name), 'cummin', problem_name, problem_type, problem_dim, processed_solver_names, ...
            fun_history, maxcv_history, merit_history, fun_inits, maxcv_inits, merit_inits, n_eval, ...
            profile_options, n_cols, default_width, default_height, summary_profile_width);
        clear warning_cleanup;
    catch Exception
        if ~profile_options.(ProfileOptionKey.SILENT.value)
            printOptiProfilerMessage('INFO', sprintf('An error occurred while plotting the history plots of the problem %s: %s', problem_name, shortenMessageForLog(Exception.message)));
        end
    end
end

function exportHistoryFigure(pdf_file, mode, problem_name, problem_type, problem_dim, solver_names, ...
    fun_history, maxcv_history, merit_history, fun_inits, maxcv_inits, merit_inits, n_eval, ...
    profile_options, n_cols, default_width, default_height, summary_profile_width)

    switch mode
        case 'combined'
            is_cum_rows = [false, true];
        case 'raw'
            is_cum_rows = false;
        case 'cummin'
            is_cum_rows = true;
        otherwise
            error("MATLAB:exportHist:unknownMode", "Unknown history plot mode: %s.", mode);
    end

    n_rows = numel(is_cum_rows);
    fig_summary = figure('Units', 'pixels', 'Position', profileFigurePosition(n_cols * summary_profile_width, n_rows * default_height), 'visible', 'off');
    fig_cleanup = onCleanup(@() close(fig_summary));
    T_summary = tiledlayout(fig_summary, n_rows, 1, 'Padding', 'compact', 'TileSpacing', 'compact');
    F_title = escapeLatexText(profile_options.(ProfileOptionKey.FEATURE_STAMP.value));
    P_title = escapeLatexText(problem_name);
    T_title = ['Solving "', P_title, '" with "', F_title, '" feature'];
    title_fontsize = min(12, 1.2 * default_width / max(length(T_title), 1));
    title_obj = title(T_summary, T_title, 'Interpreter', 'latex', 'FontSize', title_fontsize);
    set(title_obj, 'Interpreter', 'latex');

    t_summary = gobjects(n_rows, 1);
    axs_summary = gobjects(n_rows * n_cols, 1);
    i_axs = 0;
    for i = 1:n_rows
        t_summary(i) = tiledlayout(T_summary, 1, n_cols, 'Padding', 'compact', 'TileSpacing', 'compact');
        t_summary(i).Layout.Tile = i;
        for j = 1:n_cols
            i_axs = i_axs + 1;
            axs_summary(i_axs) = nexttile(t_summary(i));
        end
    end
    if strcmp(mode, 'combined')
        ylabel(t_summary(1), "History profiles", 'Interpreter', 'latex', 'FontSize', 14);
        ylabel(t_summary(2), "Cummin history profiles", 'Interpreter', 'latex', 'FontSize', 14);
    end

    for i = 1:n_rows
        row_start = (i - 1) * n_cols + 1;
        if strcmp(problem_type, 'u')
            cell_axs_summary = {axs_summary(row_start)};
        else
            cell_axs_summary = {axs_summary(row_start), axs_summary(row_start + 1), axs_summary(row_start + 2)};
        end
        drawHist(fun_history, maxcv_history, merit_history, fun_inits, maxcv_inits, merit_inits, solver_names, cell_axs_summary, is_cum_rows(i), problem_type, problem_dim, n_eval, profile_options, default_height);
    end

    exportgraphics(fig_summary, pdf_file, 'ContentType', 'vector');
    clear fig_cleanup;
end
