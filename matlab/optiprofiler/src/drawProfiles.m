function [fig_perf, fig_data, fig_log_ratio, profiles] = drawProfiles(work, problem_dimensions, solver_names, tolerance_label, cell_axs_summary, is_summary, is_perf, is_data, is_log_ratio, profile_options, profiles)
%DRAWPROFILES draws the performance, data, and log-ratio profiles.

    n_solvers = size(work, 2);

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Create the individual figures.
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    fig_perf = figure('visible', 'off');
    t_perf = tiledlayout(fig_perf, 1, 1, 'Padding', 'compact', 'TileSpacing', 'compact');
    ax_perf = nexttile(t_perf);
    fig_data = figure('visible', 'off');
    t_data = tiledlayout(fig_data, 1, 1, 'Padding', 'compact', 'TileSpacing', 'compact');
    ax_data = nexttile(t_data);
    if n_solvers > 2
        fig_log_ratio = [];
        ax_log_ratio = [];
    else
        fig_log_ratio = figure('visible', 'off');
        t_log_ratio = tiledlayout(fig_log_ratio, 1, 1, 'Padding', 'compact', 'TileSpacing', 'compact');
        ax_log_ratio = nexttile(t_log_ratio);
    end

    [x_perf, y_perf, ratio_max_perf, x_data, y_data, ratio_max_data, profiles] = getExtendedPerformancesDataProfileAxes(work, problem_dimensions, profiles);
    if n_solvers == 2
        [x_log_ratio, y_log_ratio, ratio_max_log_ratio, profiles] = getLogRatioProfileAxes(work, profiles);
    end

    if ~profile_options.(ProfileOptionKey.DRAW_PLOTS.value)
        return;
    end

    drawPerformanceDataProfiles(ax_perf, x_perf, y_perf, solver_names, profile_options);
    drawPerformanceDataProfiles(ax_data, x_data, y_data, solver_names, profile_options);
    if n_solvers == 2
        drawLogRatioProfiles(ax_log_ratio, x_log_ratio, y_log_ratio, ratio_max_log_ratio, solver_names);
    end
    % Set x-axis limits.
    set(ax_perf, 'XLim', [0.0, 1.1 * ratio_max_perf]);
    set(ax_data, 'XLim', [0.0, 1.1 * ratio_max_data]);
    % Modify x-axis ticks labels of the performance profiles.
    [ticks_perf, tickLabels_perf] = perfTicks(1.1 * ratio_max_perf);
    set(ax_perf, 'XTick', ticks_perf, 'XTickLabel', tickLabels_perf, 'TickLabelInterpreter', 'latex');
    % Modify x-axis ticks labels of the data profiles.
    [ticks_data, tickLabels_data] = dataTicks(1.1 * ratio_max_data);
    set(ax_data, 'XTick', ticks_data, 'XTickLabel', tickLabels_data, 'TickLabelInterpreter', 'latex');
    % Set x-axis labels.
    xlabel(ax_perf, 'Performance ratio', 'Interpreter', 'latex');
    xlabel(ax_data, 'Number of simplex gradients', 'Interpreter', 'latex');
    % Set y-axis labels.
    ylabel(ax_perf, ['Performance profiles (', tolerance_label, ')'], 'Interpreter', 'latex');
    ylabel(ax_data, ['Data profiles (', tolerance_label, ')'], 'Interpreter', 'latex');
    if n_solvers == 2
        % Set x-axis labels.
        xlabel(ax_log_ratio, 'Problem', 'Interpreter', 'latex');
        % Set y-axis labels.
        ylabel(ax_log_ratio, ['Log-ratio profile (', tolerance_label, ')'], 'Interpreter', 'latex');
    end

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Create the figures in summary.
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    if is_summary
        if is_perf && is_data && is_log_ratio
            % Draw the performance and data profiles.
            drawPerformanceDataProfiles(cell_axs_summary{1}, x_perf, y_perf, solver_names, profile_options);
            drawPerformanceDataProfiles(cell_axs_summary{2}, x_data, y_data, solver_names, profile_options);
            % Set x-axis limits.
            set(cell_axs_summary{1}, 'XLim', [0.0, 1.1 * ratio_max_perf]);
            set(cell_axs_summary{2}, 'XLim', [0.0, 1.1 * ratio_max_data]);
            % Modify x-axis ticks labels of the performance profiles.
            set(cell_axs_summary{1}, 'XTick', ticks_perf, 'XTickLabel', tickLabels_perf, 'TickLabelInterpreter', 'latex');
            % Modify x-axis ticks labels of the data profiles.
            set(cell_axs_summary{2}, 'XTick', ticks_data, 'XTickLabel', tickLabels_data, 'TickLabelInterpreter', 'latex');
            % Set x-axis labels.
            xlabel(cell_axs_summary{1}, 'Performance ratio', 'Interpreter', 'latex');
            xlabel(cell_axs_summary{2}, 'Number of simplex gradients', 'Interpreter', 'latex');
            % Set y-axis labels.
            ylabel(cell_axs_summary{1}, ['Performance profiles (', tolerance_label, ')'], 'Interpreter', 'latex');
            ylabel(cell_axs_summary{2}, ['Data profiles (', tolerance_label, ')'], 'Interpreter', 'latex');
            % Draw the log-ratio profiles.
            drawLogRatioProfiles(cell_axs_summary{3}, x_log_ratio, y_log_ratio, ratio_max_log_ratio, solver_names);
            % Set x-axis labels.
            xlabel(cell_axs_summary{3}, 'Problem', 'Interpreter', 'latex');
            % Set y-axis labels.
            ylabel(cell_axs_summary{3}, ['Log-ratio profile (', tolerance_label, ')'], 'Interpreter', 'latex');
        elseif is_perf && is_data
            % Draw the performance and data profiles.
            drawPerformanceDataProfiles(cell_axs_summary{1}, x_perf, y_perf, solver_names, profile_options);
            drawPerformanceDataProfiles(cell_axs_summary{2}, x_data, y_data, solver_names, profile_options);
            % Set x-axis limits.
            set(cell_axs_summary{1}, 'XLim', [0.0, 1.1 * ratio_max_perf]);
            set(cell_axs_summary{2}, 'XLim', [0.0, 1.1 * ratio_max_data]);
            % Modify x-axis ticks labels of the performance profiles.
            set(cell_axs_summary{1}, 'XTick', ticks_perf, 'XTickLabel', tickLabels_perf, 'TickLabelInterpreter', 'latex');
            % Modify x-axis ticks labels of the data profiles.
            set(cell_axs_summary{2}, 'XTick', ticks_data, 'XTickLabel', tickLabels_data, 'TickLabelInterpreter', 'latex');
            % Set x-axis labels.
            xlabel(cell_axs_summary{1}, 'Performance ratio', 'Interpreter', 'latex');
            xlabel(cell_axs_summary{2}, 'Number of simplex gradients', 'Interpreter', 'latex');
            % Set y-axis labels.
            ylabel(cell_axs_summary{1}, ['Performance profiles (', tolerance_label, ')'], 'Interpreter', 'latex');
            ylabel(cell_axs_summary{2}, ['Data profiles (', tolerance_label, ')'], 'Interpreter', 'latex');
        elseif is_perf && is_log_ratio
            % Draw the performance profiles.
            drawPerformanceDataProfiles(cell_axs_summary{1}, x_perf, y_perf, solver_names, profile_options);
            % Set x-axis limits.
            set(cell_axs_summary{1}, 'XLim', [0.0, 1.1 * ratio_max_perf]);
            % Modify x-axis ticks labels of the performance profiles.
            set(cell_axs_summary{1}, 'XTick', ticks_perf, 'XTickLabel', tickLabels_perf, 'TickLabelInterpreter', 'latex');
            % Set x-axis labels.
            xlabel(cell_axs_summary{1}, 'Performance ratio', 'Interpreter', 'latex');
            % Set y-axis labels.
            ylabel(cell_axs_summary{1}, ['Performance profiles (', tolerance_label, ')'], 'Interpreter', 'latex');
            % Draw the log-ratio profiles.
            drawLogRatioProfiles(cell_axs_summary{2}, x_log_ratio, y_log_ratio, ratio_max_log_ratio, solver_names);
            % Set x-axis labels.
            xlabel(cell_axs_summary{2}, 'Problem', 'Interpreter', 'latex');
            % Set y-axis labels.
            ylabel(cell_axs_summary{2}, ['Log-ratio profile (', tolerance_label, ')'], 'Interpreter', 'latex');
        elseif is_data && is_log_ratio
            % Draw the data profiles.
            drawPerformanceDataProfiles(cell_axs_summary{1}, x_data, y_data, solver_names, profile_options);
            % Set x-axis limits.
            set(cell_axs_summary{1}, 'XLim', [0.0, 1.1 * ratio_max_data]);
            % Modify x-axis ticks labels of the data profiles.
            set(cell_axs_summary{1}, 'XTick', ticks_data, 'XTickLabel', tickLabels_data, 'TickLabelInterpreter', 'latex');
            % Set x-axis labels.
            xlabel(cell_axs_summary{1}, 'Number of simplex gradients', 'Interpreter', 'latex');
            % Set y-axis labels.
            ylabel(cell_axs_summary{1}, ['Data profiles (', tolerance_label, ')'], 'Interpreter', 'latex');
            % Draw the log-ratio profiles.
            drawLogRatioProfiles(cell_axs_summary{2}, x_log_ratio, y_log_ratio, ratio_max_log_ratio, solver_names);
            % Set x-axis labels.
            xlabel(cell_axs_summary{2}, 'Problem', 'Interpreter', 'latex');
            % Set y-axis labels.
            ylabel(cell_axs_summary{2}, ['Log-ratio profile (', tolerance_label, ')'], 'Interpreter', 'latex');
        elseif is_perf
            % Draw the performance profiles.
            drawPerformanceDataProfiles(cell_axs_summary{1}, x_perf, y_perf, solver_names, profile_options);
            % Set x-axis limits.
            set(cell_axs_summary{1}, 'XLim', [0.0, 1.1 * ratio_max_perf]);
            % Modify x-axis ticks labels of the performance profiles.
            set(cell_axs_summary{1}, 'XTick', ticks_perf, 'XTickLabel', tickLabels_perf, 'TickLabelInterpreter', 'latex');
            % Set x-axis labels.
            xlabel(cell_axs_summary{1}, 'Performance ratio', 'Interpreter', 'latex');
            % Set y-axis labels.
            ylabel(cell_axs_summary{1}, ['Performance profiles (', tolerance_label, ')'], 'Interpreter', 'latex');
        elseif is_data
            % Draw the data profiles.
            drawPerformanceDataProfiles(cell_axs_summary{1}, x_data, y_data, solver_names, profile_options);
            % Set x-axis limits.
            set(cell_axs_summary{1}, 'XLim', [0.0, 1.1 * ratio_max_data]);
            % Modify x-axis ticks labels of the data profiles.
            set(cell_axs_summary{1}, 'XTick', ticks_data, 'XTickLabel', tickLabels_data, 'TickLabelInterpreter', 'latex');
            % Set x-axis labels.
            xlabel(cell_axs_summary{1}, 'Number of simplex gradients', 'Interpreter', 'latex');
            % Set y-axis labels.
            ylabel(cell_axs_summary{1}, ['Data profiles (', tolerance_label, ')'], 'Interpreter', 'latex');
        elseif is_log_ratio
            % Draw the log-ratio profiles.
            drawLogRatioProfiles(cell_axs_summary{1}, x_log_ratio, y_log_ratio, ratio_max_log_ratio, solver_names);
            % Set x-axis labels.
            xlabel(cell_axs_summary{1}, 'Problem', 'Interpreter', 'latex');
            % Set y-axis labels.
            ylabel(cell_axs_summary{1}, ['Log-ratio profile (', tolerance_label, ')'], 'Interpreter', 'latex');
        end
    end
end

function [ticks, tickLabels] = perfTicks(ratio_cut_perf)
    % Set the x-axis ticks and tick labels for the performance profiles.

    if ratio_cut_perf >= 5
        max_power = floor(ratio_cut_perf);
        ticks = linspace(0, max_power, 6);
        ticks(2:end-1) = round(ticks(2:end-1));
        ticks = unique(ticks, 'stable');
    elseif ratio_cut_perf >= 1
        max_power = floor(ratio_cut_perf);
        ticks = (0:1:max_power);
    elseif ratio_cut_perf >= 1e-3
        ticks = [0 ratio_cut_perf];
    else
        ticks = [0];
    end

    tickLabels = arrayfun(@(x) num2str(2 ^ x), ticks, 'UniformOutput', false);
end

function [ticks, tickLabels] = dataTicks(ratio_cut_data)
    % Set the x-axis ticks and tick labels for the data profiles.

    if ratio_cut_data >= 5
        max_power = floor(ratio_cut_data);
        ticks = linspace(0, max_power, 6);
        ticks(2:end-1) = round(ticks(2:end-1));
        ticks = unique(ticks, 'stable');
    elseif ratio_cut_data >= 1
        max_power = floor(ratio_cut_data);
        ticks = (0:1:max_power);
    elseif ratio_cut_data >= 1e-1
        ticks = [0 ratio_cut_data];
    else
        ticks = [0];
    end

    tickLabels = arrayfun(@(x) num2str(2 ^ x - 1), ticks, 'UniformOutput', false);
end