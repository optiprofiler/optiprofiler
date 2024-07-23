function [fig_perf_hist, fig_perf_out, fig_data_hist, fig_data_out, fig_log_ratio_hist, fig_log_ratio_out] = drawProfiles(work_hist, work_out, problem_dimensions, labels, tolerance_label, cell_axs_summary, is_perf, is_data, is_log_ratio, profile_options)
%DRAWPROFILES draws the performance, data, and log-ratio profiles with respect to the whole history of function values or the returned values from solvers, respectively.

    n_solvers = size(work_hist, 2);

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Create the individual figures.
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    fig_perf_hist = figure('visible', 'off');
    t_perf_hist = tiledlayout(fig_perf_hist, 1, 1, 'Padding', 'compact', 'TileSpacing', 'compact');
    ax_perf_hist = nexttile(t_perf_hist);
    fig_perf_out = figure('visible', 'off');
    t_perf_out = tiledlayout(fig_perf_out, 1, 1, 'Padding', 'compact', 'TileSpacing', 'compact');
    ax_perf_out = nexttile(t_perf_out);
    fig_data_hist = figure('visible', 'off');
    t_data_hist = tiledlayout(fig_data_hist, 1, 1, 'Padding', 'compact', 'TileSpacing', 'compact');
    ax_data_hist = nexttile(t_data_hist);
    fig_data_out = figure('visible', 'off');
    t_data_out = tiledlayout(fig_data_out, 1, 1, 'Padding', 'compact', 'TileSpacing', 'compact');
    ax_data_out = nexttile(t_data_out);
    if n_solvers > 2
        fig_log_ratio_hist = [];
        ax_log_ratio_hist = [];
        fig_log_ratio_out = [];
        ax_log_ratio_out = [];
    else
        fig_log_ratio_hist = figure('visible', 'off');
        t_log_ratio_hist = tiledlayout(fig_log_ratio_hist, 1, 1, 'Padding', 'compact', 'TileSpacing', 'compact');
        ax_log_ratio_hist = nexttile(t_log_ratio_hist);
        fig_log_ratio_out = figure('visible', 'off');
        t_log_ratio_out = tiledlayout(fig_log_ratio_out, 1, 1, 'Padding', 'compact', 'TileSpacing', 'compact');
        ax_log_ratio_out = nexttile(t_log_ratio_out);
    end
    % Draw the performance and data profiles.
    [x_perf_hist, y_perf_hist, ratio_max_perf_hist, x_data_hist, y_data_hist, ratio_max_data_hist] = getExtendedPerformancesDataProfileAxes(work_hist, problem_dimensions);
    [x_perf_out, y_perf_out, ratio_max_perf_out, x_data_out, y_data_out, ratio_max_data_out] = getExtendedPerformancesDataProfileAxes(work_out, problem_dimensions);
    drawPerformanceDataProfiles(ax_perf_hist, x_perf_hist, y_perf_hist, labels, profile_options);
    drawPerformanceDataProfiles(ax_perf_out, x_perf_out, y_perf_out, labels, profile_options);
    drawPerformanceDataProfiles(ax_data_hist, x_data_hist, y_data_hist, labels, profile_options);
    drawPerformanceDataProfiles(ax_data_out, x_data_out, y_data_out, labels, profile_options);
    % Set x-axis limits.
    set(ax_perf_hist, 'XLim', [0.0, 1.1 * ratio_max_perf_hist]);
    set(ax_perf_out, 'XLim', [0.0, 1.1 * ratio_max_perf_out]);
    set(ax_data_hist, 'XLim', [0.0, 1.1 * ratio_max_data_hist]);
    set(ax_data_out, 'XLim', [0.0, 1.1 * ratio_max_data_out]);
    % Modify x-axis ticks labels of the performance profiles.
    [ticks_perf_hist, tickLabels_perf_hist] = perfTicks(1.1 * ratio_max_perf_hist);
    set(ax_perf_hist, 'XTick', ticks_perf_hist, 'XTickLabel', tickLabels_perf_hist, 'TickLabelInterpreter', 'latex');
    [ticks_perf_out, tickLabels_perf_out] = perfTicks(1.1 * ratio_max_perf_out);
    set(ax_perf_out, 'XTick', ticks_perf_out, 'XTickLabel', tickLabels_perf_out, 'TickLabelInterpreter', 'latex');
    % Modify x-axis ticks labels of the data profiles.
    [ticks_data_hist, tickLabels_data_hist] = dataTicks(1.1 * ratio_max_data_hist);
    set(ax_data_hist, 'XTick', ticks_data_hist, 'XTickLabel', tickLabels_data_hist, 'TickLabelInterpreter', 'latex');
    [ticks_data_out, tickLabels_data_out] = dataTicks(1.1 * ratio_max_data_out);
    set(ax_data_out, 'XTick', ticks_data_out, 'XTickLabel', tickLabels_data_out, 'TickLabelInterpreter', 'latex');
    % Set x-axis labels.
    xlabel(ax_perf_hist, 'Performance ratio', 'Interpreter', 'latex');
    xlabel(ax_perf_out, 'Performance ratio', 'Interpreter', 'latex');
    xlabel(ax_data_hist, 'Number of simplex gradients', 'Interpreter', 'latex');
    xlabel(ax_data_out, 'Number of simplex gradients', 'Interpreter', 'latex');
    % Set y-axis labels.
    ylabel(ax_perf_hist, ['Performance profiles (', tolerance_label, ')'], 'Interpreter', 'latex');
    ylabel(ax_perf_out, ['Performance profiles (', tolerance_label, ')'], 'Interpreter', 'latex');
    ylabel(ax_data_hist, ['Data profiles (', tolerance_label, ')'], 'Interpreter', 'latex');
    ylabel(ax_data_out, ['Data profiles (', tolerance_label, ')'], 'Interpreter', 'latex');
    % Draw the log-ratio profiles.
    if n_solvers <= 2
        copy_work_hist1 = work_hist;
        drawLogRatioProfiles(ax_log_ratio_hist, copy_work_hist1, labels);
        copy_work_out1 = work_out;
        drawLogRatioProfiles(ax_log_ratio_out, copy_work_out1, labels);
        % Set x-axis labels.
        xlabel(ax_log_ratio_hist, 'Problem', 'Interpreter', 'latex');
        xlabel(ax_log_ratio_out, 'Problem', 'Interpreter', 'latex');
        % Set y-axis labels.
        ylabel(ax_log_ratio_hist, ['Log-ratio profile (', tolerance_label, ')'], 'Interpreter', 'latex');
        ylabel(ax_log_ratio_out, ['Log-ratio profile (', tolerance_label, ')'], 'Interpreter', 'latex');
    end

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Create the figures in summary.
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    if is_perf && is_data && is_log_ratio
        % Draw the performance and data profiles.
        drawPerformanceDataProfiles(cell_axs_summary{1}, x_perf_hist, y_perf_hist, labels, profile_options);
        drawPerformanceDataProfiles(cell_axs_summary{2}, x_perf_out, y_perf_out, labels, profile_options);
        drawPerformanceDataProfiles(cell_axs_summary{3}, x_data_hist, y_data_hist, labels, profile_options);
        drawPerformanceDataProfiles(cell_axs_summary{4}, x_data_out, y_data_out, labels, profile_options);
        % Set x-axis limits.
        set(cell_axs_summary{1}, 'XLim', [0.0, 1.1 * ratio_max_perf_hist]);
        set(cell_axs_summary{2}, 'XLim', [0.0, 1.1 * ratio_max_perf_out]);
        set(cell_axs_summary{3}, 'XLim', [0.0, 1.1 * ratio_max_data_hist]);
        set(cell_axs_summary{4}, 'XLim', [0.0, 1.1 * ratio_max_data_out]);
        % Modify x-axis ticks labels of the performance profiles.
        set(cell_axs_summary{1}, 'XTick', ticks_perf_hist, 'XTickLabel', tickLabels_perf_hist, 'TickLabelInterpreter', 'latex');
        set(cell_axs_summary{2}, 'XTick', ticks_perf_out, 'XTickLabel', tickLabels_perf_out, 'TickLabelInterpreter', 'latex');
        % Modify x-axis ticks labels of the data profiles.
        set(cell_axs_summary{3}, 'XTick', ticks_data_hist, 'XTickLabel', tickLabels_data_hist, 'TickLabelInterpreter', 'latex');
        set(cell_axs_summary{4}, 'XTick', ticks_data_out, 'XTickLabel', tickLabels_data_out, 'TickLabelInterpreter', 'latex');
        % Set x-axis labels.
        xlabel(cell_axs_summary{1}, 'Performance ratio', 'Interpreter', 'latex');
        xlabel(cell_axs_summary{2}, 'Performance ratio', 'Interpreter', 'latex');
        xlabel(cell_axs_summary{3}, 'Number of simplex gradients', 'Interpreter', 'latex');
        xlabel(cell_axs_summary{4}, 'Number of simplex gradients', 'Interpreter', 'latex');
        % Set y-axis labels.
        ylabel(cell_axs_summary{1}, ['Performance profiles (', tolerance_label, ')'], 'Interpreter', 'latex');
        ylabel(cell_axs_summary{2}, ['Performance profiles (', tolerance_label, ')'], 'Interpreter', 'latex');
        ylabel(cell_axs_summary{3}, ['Data profiles (', tolerance_label, ')'], 'Interpreter', 'latex');
        ylabel(cell_axs_summary{4}, ['Data profiles (', tolerance_label, ')'], 'Interpreter', 'latex');
        % Draw the log-ratio profiles.
        copy_work_hist2 = work_hist;
        drawLogRatioProfiles(cell_axs_summary{5}, copy_work_hist2, labels);
        copy_work_out2 = work_out;
        drawLogRatioProfiles(cell_axs_summary{6}, copy_work_out2, labels);
        % Set x-axis labels.
        xlabel(cell_axs_summary{5}, 'Problem', 'Interpreter', 'latex');
        xlabel(cell_axs_summary{6}, 'Problem', 'Interpreter', 'latex');
        % Set y-axis labels.
        ylabel(cell_axs_summary{5}, ['Log-ratio profile (', tolerance_label, ')'], 'Interpreter', 'latex');
        ylabel(cell_axs_summary{6}, ['Log-ratio profile (', tolerance_label, ')'], 'Interpreter', 'latex');
    elseif is_perf && is_data
        % Draw the performance and data profiles.
        drawPerformanceDataProfiles(cell_axs_summary{1}, x_perf_hist, y_perf_hist, labels, profile_options);
        drawPerformanceDataProfiles(cell_axs_summary{2}, x_perf_out, y_perf_out, labels, profile_options);
        drawPerformanceDataProfiles(cell_axs_summary{3}, x_data_hist, y_data_hist, labels, profile_options);
        drawPerformanceDataProfiles(cell_axs_summary{4}, x_data_out, y_data_out, labels, profile_options);
        % Set x-axis limits.
        set(cell_axs_summary{1}, 'XLim', [0.0, 1.1 * ratio_max_perf_hist]);
        set(cell_axs_summary{2}, 'XLim', [0.0, 1.1 * ratio_max_perf_out]);
        set(cell_axs_summary{3}, 'XLim', [0.0, 1.1 * ratio_max_data_hist]);
        set(cell_axs_summary{4}, 'XLim', [0.0, 1.1 * ratio_max_data_out]);
        % Modify x-axis ticks labels of the performance profiles.
        set(cell_axs_summary{1}, 'XTick', ticks_perf_hist, 'XTickLabel', tickLabels_perf_hist, 'TickLabelInterpreter', 'latex');
        set(cell_axs_summary{2}, 'XTick', ticks_perf_out, 'XTickLabel', tickLabels_perf_out, 'TickLabelInterpreter', 'latex');
        % Modify x-axis ticks labels of the data profiles.
        set(cell_axs_summary{3}, 'XTick', ticks_data_hist, 'XTickLabel', tickLabels_data_hist, 'TickLabelInterpreter', 'latex');
        set(cell_axs_summary{4}, 'XTick', ticks_data_out, 'XTickLabel', tickLabels_data_out, 'TickLabelInterpreter', 'latex');
        % Set x-axis labels.
        xlabel(cell_axs_summary{1}, 'Performance ratio', 'Interpreter', 'latex');
        xlabel(cell_axs_summary{2}, 'Performance ratio', 'Interpreter', 'latex');
        xlabel(cell_axs_summary{3}, 'Number of simplex gradients', 'Interpreter', 'latex');
        xlabel(cell_axs_summary{4}, 'Number of simplex gradients', 'Interpreter', 'latex');
        % Set y-axis labels.
        ylabel(cell_axs_summary{1}, ['Performance profiles (', tolerance_label, ')'], 'Interpreter', 'latex');
        ylabel(cell_axs_summary{2}, ['Performance profiles (', tolerance_label, ')'], 'Interpreter', 'latex');
        ylabel(cell_axs_summary{3}, ['Data profiles (', tolerance_label, ')'], 'Interpreter', 'latex');
        ylabel(cell_axs_summary{4}, ['Data profiles (', tolerance_label, ')'], 'Interpreter', 'latex');
    elseif is_perf && is_log_ratio
        % Draw the performance profiles.
        drawPerformanceDataProfiles(cell_axs_summary{1}, x_perf_hist, y_perf_hist, labels, profile_options);
        drawPerformanceDataProfiles(cell_axs_summary{2}, x_perf_out, y_perf_out, labels, profile_options);
        % Set x-axis limits.
        set(cell_axs_summary{1}, 'XLim', [0.0, 1.1 * ratio_max_perf_hist]);
        set(cell_axs_summary{2}, 'XLim', [0.0, 1.1 * ratio_max_perf_out]);
        % Modify x-axis ticks labels of the performance profiles.
        set(cell_axs_summary{1}, 'XTick', ticks_perf_hist, 'XTickLabel', tickLabels_perf_hist, 'TickLabelInterpreter', 'latex');
        set(cell_axs_summary{2}, 'XTick', ticks_perf_out, 'XTickLabel', tickLabels_perf_out, 'TickLabelInterpreter', 'latex');
        % Set x-axis labels.
        xlabel(cell_axs_summary{1}, 'Performance ratio', 'Interpreter', 'latex');
        xlabel(cell_axs_summary{2}, 'Performance ratio', 'Interpreter', 'latex');
        % Set y-axis labels.
        ylabel(cell_axs_summary{1}, ['Performance profiles (', tolerance_label, ')'], 'Interpreter', 'latex');
        ylabel(cell_axs_summary{2}, ['Performance profiles (', tolerance_label, ')'], 'Interpreter', 'latex');
        % Draw the log-ratio profiles.
        copy_work_hist2 = work_hist;
        drawLogRatioProfiles(cell_axs_summary{3}, copy_work_hist2, labels);
        copy_work_out2 = work_out;
        drawLogRatioProfiles(cell_axs_summary{4}, copy_work_out2, labels);
        % Set x-axis labels.
        xlabel(cell_axs_summary{3}, 'Problem', 'Interpreter', 'latex');
        xlabel(cell_axs_summary{4}, 'Problem', 'Interpreter', 'latex');
        % Set y-axis labels.
        ylabel(cell_axs_summary{3}, ['Log-ratio profile (', tolerance_label, ')'], 'Interpreter', 'latex');
        ylabel(cell_axs_summary{4}, ['Log-ratio profile (', tolerance_label, ')'], 'Interpreter', 'latex');
    elseif is_data && is_log_ratio
        % Draw the data profiles.
        drawPerformanceDataProfiles(cell_axs_summary{1}, x_data_hist, y_data_hist, labels, profile_options);
        drawPerformanceDataProfiles(cell_axs_summary{2}, x_data_out, y_data_out, labels, profile_options);
        % Set x-axis limits.
        set(cell_axs_summary{1}, 'XLim', [0.0, 1.1 * ratio_max_data_hist]);
        set(cell_axs_summary{2}, 'XLim', [0.0, 1.1 * ratio_max_data_out]);
        % Modify x-axis ticks labels of the data profiles.
        set(cell_axs_summary{1}, 'XTick', ticks_data_hist, 'XTickLabel', tickLabels_data_hist, 'TickLabelInterpreter', 'latex');
        set(cell_axs_summary{2}, 'XTick', ticks_data_out, 'XTickLabel', tickLabels_data_out, 'TickLabelInterpreter', 'latex');
        % Set x-axis labels.
        xlabel(cell_axs_summary{1}, 'Number of simplex gradients', 'Interpreter', 'latex');
        xlabel(cell_axs_summary{2}, 'Number of simplex gradients', 'Interpreter', 'latex');
        % Set y-axis labels.
        ylabel(cell_axs_summary{1}, ['Data profiles (', tolerance_label, ')'], 'Interpreter', 'latex');
        ylabel(cell_axs_summary{2}, ['Data profiles (', tolerance_label, ')'], 'Interpreter', 'latex');
        % Draw the log-ratio profiles.
        copy_work_hist2 = work_hist;
        drawLogRatioProfiles(cell_axs_summary{3}, copy_work_hist2, labels);
        copy_work_out2 = work_out;
        drawLogRatioProfiles(cell_axs_summary{4}, copy_work_out2, labels);
        % Set x-axis labels.
        xlabel(cell_axs_summary{3}, 'Problem', 'Interpreter', 'latex');
        xlabel(cell_axs_summary{4}, 'Problem', 'Interpreter', 'latex');
        % Set y-axis labels.
        ylabel(cell_axs_summary{3}, ['Log-ratio profile (', tolerance_label, ')'], 'Interpreter', 'latex');
        ylabel(cell_axs_summary{4}, ['Log-ratio profile (', tolerance_label, ')'], 'Interpreter', 'latex');
    elseif is_perf
        % Draw the performance profiles.
        drawPerformanceDataProfiles(cell_axs_summary{1}, x_perf_hist, y_perf_hist, labels, profile_options);
        drawPerformanceDataProfiles(cell_axs_summary{2}, x_perf_out, y_perf_out, labels, profile_options);
        % Set x-axis limits.
        set(cell_axs_summary{1}, 'XLim', [0.0, 1.1 * ratio_max_perf_hist]);
        set(cell_axs_summary{2}, 'XLim', [0.0, 1.1 * ratio_max_perf_out]);
        % Modify x-axis ticks labels of the performance profiles.
        set(cell_axs_summary{1}, 'XTick', ticks_perf_hist, 'XTickLabel', tickLabels_perf_hist, 'TickLabelInterpreter', 'latex');
        set(cell_axs_summary{2}, 'XTick', ticks_perf_out, 'XTickLabel', tickLabels_perf_out, 'TickLabelInterpreter', 'latex');
        % Set x-axis labels.
        xlabel(cell_axs_summary{1}, 'Performance ratio', 'Interpreter', 'latex');
        xlabel(cell_axs_summary{2}, 'Performance ratio', 'Interpreter', 'latex');
        % Set y-axis labels.
        ylabel(cell_axs_summary{1}, ['Performance profiles (', tolerance_label, ')'], 'Interpreter', 'latex');
        ylabel(cell_axs_summary{2}, ['Performance profiles (', tolerance_label, ')'], 'Interpreter', 'latex');
    elseif is_data
        % Draw the data profiles.
        drawPerformanceDataProfiles(cell_axs_summary{1}, x_data_hist, y_data_hist, labels, profile_options);
        drawPerformanceDataProfiles(cell_axs_summary{2}, x_data_out, y_data_out, labels, profile_options);
        % Set x-axis limits.
        set(cell_axs_summary{1}, 'XLim', [0.0, 1.1 * ratio_max_data_hist]);
        set(cell_axs_summary{2}, 'XLim', [0.0, 1.1 * ratio_max_data_out]);
        % Modify x-axis ticks labels of the data profiles.
        set(cell_axs_summary{1}, 'XTick', ticks_data_hist, 'XTickLabel', tickLabels_data_hist, 'TickLabelInterpreter', 'latex');
        set(cell_axs_summary{2}, 'XTick', ticks_data_out, 'XTickLabel', tickLabels_data_out, 'TickLabelInterpreter', 'latex');
        % Set x-axis labels.
        xlabel(cell_axs_summary{1}, 'Number of simplex gradients', 'Interpreter', 'latex');
        xlabel(cell_axs_summary{2}, 'Number of simplex gradients', 'Interpreter', 'latex');
        % Set y-axis labels.
        ylabel(cell_axs_summary{1}, ['Data profiles (', tolerance_label, ')'], 'Interpreter', 'latex');
        ylabel(cell_axs_summary{2}, ['Data profiles (', tolerance_label, ')'], 'Interpreter', 'latex');
    elseif is_log_ratio
        % Draw the log-ratio profiles.
        copy_work_hist2 = work_hist;
        drawLogRatioProfiles(cell_axs_summary{1}, copy_work_hist2, labels);
        copy_work_out2 = work_out;
        drawLogRatioProfiles(cell_axs_summary{2}, copy_work_out2, labels);
        % Set x-axis labels.
        xlabel(cell_axs_summary{1}, 'Problem', 'Interpreter', 'latex');
        xlabel(cell_axs_summary{2}, 'Problem', 'Interpreter', 'latex');
        % Set y-axis labels.
        ylabel(cell_axs_summary{1}, ['Log-ratio profile (', tolerance_label, ')'], 'Interpreter', 'latex');
        ylabel(cell_axs_summary{2}, ['Log-ratio profile (', tolerance_label, ')'], 'Interpreter', 'latex');
    end

end