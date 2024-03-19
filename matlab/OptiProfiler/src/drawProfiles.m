function [fig_perf_hist, fig_perf_ret, fig_data_hist, fig_data_ret, fig_log_ratio_hist, fig_log_ratio_ret] = drawProfiles(work_hist, work_ret, problem_dimensions, labels, tolerance_label, cell_axs_summary)
%DRAWPROFILES draws the performance, data, and log-ratio profiles with respect to the whole history of function values or the returned values from solvers, respectively.

    n_solvers = size(work_hist, 2);

    % Create the individual figures.
    fig_perf_hist = figure('visible', 'off');
    t_perf_hist = tiledlayout(fig_perf_hist, 1, 1, 'Padding', 'compact', 'TileSpacing', 'compact');
    ax_perf_hist = nexttile(t_perf_hist);
    fig_perf_ret = figure('visible', 'off');
    t_perf_ret = tiledlayout(fig_perf_ret, 1, 1, 'Padding', 'compact', 'TileSpacing', 'compact');
    ax_perf_ret = nexttile(t_perf_ret);
    fig_data_hist = figure('visible', 'off');
    t_data_hist = tiledlayout(fig_data_hist, 1, 1, 'Padding', 'compact', 'TileSpacing', 'compact');
    ax_data_hist = nexttile(t_data_hist);
    fig_data_ret = figure('visible', 'off');
    t_data_ret = tiledlayout(fig_data_ret, 1, 1, 'Padding', 'compact', 'TileSpacing', 'compact');
    ax_data_ret = nexttile(t_data_ret);
    if n_solvers > 2
        fig_log_ratio_hist = [];
        ax_log_ratio_hist = [];
        fig_log_ratio_ret = [];
        ax_log_ratio_ret = [];
    else
        fig_log_ratio_hist = figure('visible', 'off');
        t_log_ratio_hist = tiledlayout(fig_log_ratio_hist, 1, 1, 'Padding', 'compact', 'TileSpacing', 'compact');
        ax_log_ratio_hist = nexttile(t_log_ratio_hist);
        fig_log_ratio_ret = figure('visible', 'off');
        t_log_ratio_ret = tiledlayout(fig_log_ratio_ret, 1, 1, 'Padding', 'compact', 'TileSpacing', 'compact');
        ax_log_ratio_ret = nexttile(t_log_ratio_ret);
    end

    % Draw the performance and data profiles.
    [x_perf_hist, y_perf_hist, ratio_max_perf_hist, x_data_hist, y_data_hist, ratio_max_data_hist] = getExtendedPerformancesDataProfileAxes(work_hist, problem_dimensions);
    [x_perf_ret, y_perf_ret, ratio_max_perf_ret, x_data_ret, y_data_ret, ratio_max_data_ret] = getExtendedPerformancesDataProfileAxes(work_ret, problem_dimensions);
    drawPerformanceDataProfiles(ax_perf_hist, x_perf_hist, y_perf_hist, labels);
    drawPerformanceDataProfiles(ax_perf_ret, x_perf_ret, y_perf_ret, labels);
    drawPerformanceDataProfiles(ax_data_hist, x_data_hist, y_data_hist, labels);
    drawPerformanceDataProfiles(ax_data_ret, x_data_ret, y_data_ret, labels);
    drawPerformanceDataProfiles(cell_axs_summary{1}, x_perf_hist, y_perf_hist, labels);
    drawPerformanceDataProfiles(cell_axs_summary{2}, x_perf_ret, y_perf_ret, labels);
    drawPerformanceDataProfiles(cell_axs_summary{3}, x_data_hist, y_data_hist, labels);
    drawPerformanceDataProfiles(cell_axs_summary{4}, x_data_ret, y_data_ret, labels);

    % Set x-axis limits.
    set(ax_perf_hist, 'XLim', [0.0, 1.1 * ratio_max_perf_hist]);
    set(ax_perf_ret, 'XLim', [0.0, 1.1 * ratio_max_perf_ret]);
    set(ax_data_hist, 'XLim', [0.0, 1.1 * ratio_max_data_hist]);
    set(ax_data_ret, 'XLim', [0.0, 1.1 * ratio_max_data_ret]);
    set(cell_axs_summary{1}, 'XLim', [0.0, 1.1 * ratio_max_perf_hist]);
    set(cell_axs_summary{2}, 'XLim', [0.0, 1.1 * ratio_max_perf_ret]);
    set(cell_axs_summary{3}, 'XLim', [0.0, 1.1 * ratio_max_data_hist]);
    set(cell_axs_summary{4}, 'XLim', [0.0, 1.1 * ratio_max_data_ret]);

    % Modify x-axis ticks labels of the performance profiles.
    xtl_perf_hist = arrayfun(@(x) ['$' num2str(2 ^ x) '$'], get(ax_perf_hist, 'XTick'), 'UniformOutput', false);
    set(ax_perf_hist, 'XTickLabel', xtl_perf_hist, 'TickLabelInterpreter', 'latex');
    set(cell_axs_summary{1}, 'XTickLabel', xtl_perf_hist, 'TickLabelInterpreter', 'latex');
    xtl_perf_ret = arrayfun(@(x) ['$' num2str(2 ^ x) '$'], get(ax_perf_ret, 'XTick'), 'UniformOutput', false);
    set(ax_perf_ret, 'XTickLabel', xtl_perf_ret, 'TickLabelInterpreter', 'latex');
    set(cell_axs_summary{2}, 'XTickLabel', xtl_perf_ret, 'TickLabelInterpreter', 'latex');

    % Modify x-axis ticks labels of the data profiles.
    xtl_data_hist = arrayfun(@(x) ['$', num2str(2 ^ x - 1), '$'], get(ax_data_hist, 'XTick'), 'UniformOutput', false);
    set(ax_data_hist, 'XTickLabel', xtl_data_hist, 'TickLabelInterpreter', 'latex');
    set(cell_axs_summary{3}, 'XTickLabel', xtl_data_hist, 'TickLabelInterpreter', 'latex');
    xtl_data_ret = arrayfun(@(x) ['$', num2str(2 ^ x - 1), '$'], get(ax_data_ret, 'XTick'), 'UniformOutput', false);
    set(ax_data_ret, 'XTickLabel', xtl_data_ret, 'TickLabelInterpreter', 'latex');
    set(cell_axs_summary{4}, 'XTickLabel', xtl_data_ret, 'TickLabelInterpreter', 'latex');

    % Set x-axis labels.
    xlabel(ax_perf_hist, 'Performance ratio', 'Interpreter', 'latex');
    xlabel(ax_perf_ret, 'Performance ratio', 'Interpreter', 'latex');
    xlabel(ax_data_hist, 'Number of simplex gradients', 'Interpreter', 'latex');
    xlabel(ax_data_ret, 'Number of simplex gradients', 'Interpreter', 'latex');
    xlabel(cell_axs_summary{1}, 'Performance ratio', 'Interpreter', 'latex');
    xlabel(cell_axs_summary{2}, 'Performance ratio', 'Interpreter', 'latex');
    xlabel(cell_axs_summary{3}, 'Number of simplex gradients', 'Interpreter', 'latex');
    xlabel(cell_axs_summary{4}, 'Number of simplex gradients', 'Interpreter', 'latex');

    % Set y-axis labels.
    ylabel(ax_perf_hist, ['Performance profiles (', tolerance_label, ')'], 'Interpreter', 'latex');
    ylabel(ax_perf_ret, ['Performance profiles (', tolerance_label, ')'], 'Interpreter', 'latex');
    ylabel(ax_data_hist, ['Data profiles (', tolerance_label, ')'], 'Interpreter', 'latex');
    ylabel(ax_data_ret, ['Data profiles (', tolerance_label, ')'], 'Interpreter', 'latex');
    ylabel(cell_axs_summary{1}, ['Performance profiles (', tolerance_label, ')'], 'Interpreter', 'latex');
    ylabel(cell_axs_summary{2}, ['Performance profiles (', tolerance_label, ')'], 'Interpreter', 'latex');
    ylabel(cell_axs_summary{3}, ['Data profiles (', tolerance_label, ')'], 'Interpreter', 'latex');
    ylabel(cell_axs_summary{4}, ['Data profiles (', tolerance_label, ')'], 'Interpreter', 'latex');
    
    % Draw the log-ratio profiles.
    if n_solvers <= 2
        copy_work_hist1 = work_hist;
        drawLogRatioProfiles(ax_log_ratio_hist, copy_work_hist1, labels);
        copy_work_hist2 = work_hist;
        drawLogRatioProfiles(cell_axs_summary{5}, copy_work_hist2, labels);
        copy_work_ret1 = work_ret;
        drawLogRatioProfiles(ax_log_ratio_ret, copy_work_ret1, labels);
        copy_work_ret2 = work_ret;
        drawLogRatioProfiles(cell_axs_summary{6}, copy_work_ret2, labels);
        % Set x-axis labels.
        xlabel(ax_log_ratio_hist, 'Problem', 'Interpreter', 'latex');
        xlabel(cell_axs_summary{5}, 'Problem', 'Interpreter', 'latex');
        xlabel(ax_log_ratio_ret, 'Problem', 'Interpreter', 'latex');
        xlabel(cell_axs_summary{6}, 'Problem', 'Interpreter', 'latex');
        % Set y-axis labels.
        ylabel(ax_log_ratio_hist, ['Log-ratio profile (', tolerance_label, ')'], 'Interpreter', 'latex');
        ylabel(cell_axs_summary{5}, ['Log-ratio profile (', tolerance_label, ')'], 'Interpreter', 'latex');
        ylabel(ax_log_ratio_ret, ['Log-ratio profile (', tolerance_label, ')'], 'Interpreter', 'latex');
        ylabel(cell_axs_summary{6}, ['Log-ratio profile (', tolerance_label, ')'], 'Interpreter', 'latex');
    end

end