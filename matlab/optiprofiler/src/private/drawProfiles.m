function [fig_perf, fig_data, fig_log_ratio, curves] = drawProfiles(work, problem_dimensions, solver_names, tolerance_latex, cell_axs_summary, is_summary, is_perf, is_data, is_log_ratio, profile_options, curves)
%DRAWPROFILES draws the performance, data, and log-ratio profiles.

    solver_names = cellfun(@(s) strrep(s, '_', '\_'), solver_names, 'UniformOutput', false);
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

    [x_perf, y_perf, ratio_max_perf, x_data, y_data, ratio_max_data, curves] = getExtendedPerformancesDataProfileAxes(work, problem_dimensions, profile_options, curves);
    if n_solvers == 2
        [x_log_ratio, y_log_ratio, ratio_max_log_ratio, n_solvers_fail, curves] = getLogRatioProfileAxes(work, curves);
    end

    if profile_options.(ProfileOptionKey.SCORE_ONLY.value)
        return;
    end

    drawPerfDetail(ax_perf, x_perf, y_perf, ratio_max_perf, solver_names, profile_options, tolerance_latex);
    drawDataDetail(ax_data, x_data, y_data, ratio_max_data, solver_names, profile_options, tolerance_latex);
    if n_solvers == 2
        drawLogRatioDetail(ax_log_ratio, x_log_ratio, y_log_ratio, ratio_max_log_ratio, n_solvers_fail, solver_names, profile_options, tolerance_latex);
    end

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Create the figures in summary.
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    if is_summary
        if is_perf && is_data && is_log_ratio
            drawPerfDetail(cell_axs_summary{1}, x_perf, y_perf, ratio_max_perf, solver_names, profile_options, tolerance_latex);
            drawDataDetail(cell_axs_summary{2}, x_data, y_data, ratio_max_data, solver_names, profile_options, tolerance_latex);
            drawLogRatioDetail(cell_axs_summary{3}, x_log_ratio, y_log_ratio, ratio_max_log_ratio, n_solvers_fail, solver_names, profile_options, tolerance_latex);
        elseif is_perf && is_data
            drawPerfDetail(cell_axs_summary{1}, x_perf, y_perf, ratio_max_perf, solver_names, profile_options, tolerance_latex);
            drawDataDetail(cell_axs_summary{2}, x_data, y_data, ratio_max_data, solver_names, profile_options, tolerance_latex);
        elseif is_perf && is_log_ratio
            drawPerfDetail(cell_axs_summary{1}, x_perf, y_perf, ratio_max_perf, solver_names, profile_options, tolerance_latex);
            drawLogRatioDetail(cell_axs_summary{2}, x_log_ratio, y_log_ratio, ratio_max_log_ratio, n_solvers_fail, solver_names, profile_options, tolerance_latex);
        elseif is_data && is_log_ratio
            drawDataDetail(cell_axs_summary{1}, x_data, y_data, ratio_max_data, solver_names, profile_options, tolerance_latex);
            drawLogRatioDetail(cell_axs_summary{2}, x_log_ratio, y_log_ratio, ratio_max_log_ratio, n_solvers_fail, solver_names, profile_options, tolerance_latex);
        elseif is_perf
            drawPerfDetail(cell_axs_summary{1}, x_perf, y_perf, ratio_max_perf, solver_names, profile_options, tolerance_latex);
        elseif is_data
            drawDataDetail(cell_axs_summary{1}, x_data, y_data, ratio_max_data, solver_names, profile_options, tolerance_latex);
        elseif is_log_ratio
            drawLogRatioDetail(cell_axs_summary{1}, x_log_ratio, y_log_ratio, ratio_max_log_ratio, n_solvers_fail, solver_names, profile_options, tolerance_latex);
        end
    end
end

function drawPerfDetail(ax_perf, x_perf, y_perf, ratio_max_perf, solver_names, profile_options, tolerance_latex)
    drawPerformanceDataProfiles(ax_perf, x_perf, y_perf, solver_names, profile_options);
    % Set x-axis limits.
    if profile_options.(ProfileOptionKey.SEMILOGX.value)
        set(ax_perf, 'XLim', [0.0, 1.1 * ratio_max_perf]);
    else
        set(ax_perf, 'XLim', [1.0, 1.1 * ratio_max_perf]);
    end
    % Modify x-axis ticks labels of the performance profiles.
    [ticks_perf, tickLabels_perf] = perfTicks(1.1 * ratio_max_perf, profile_options.(ProfileOptionKey.SEMILOGX.value));
    set(ax_perf, 'XTick', ticks_perf, 'XTickLabel', tickLabels_perf, 'TickLabelInterpreter', 'latex');
    % Set x-axis labels.
    if ~isempty(profile_options.(ProfileOptionKey.XLABEL_PERFORMANCE_PROFILE.value))
        xlabel_str = profile_options.(ProfileOptionKey.XLABEL_PERFORMANCE_PROFILE.value);
        xlabel(ax_perf, xlabel_str, 'Interpreter', 'latex');
    end
    % Set y-axis labels.
    if ~isempty(profile_options.(ProfileOptionKey.YLABEL_PERFORMANCE_PROFILE.value))
        ylabel_str = sprintf(profile_options.(ProfileOptionKey.YLABEL_PERFORMANCE_PROFILE.value), tolerance_latex);
        ylabel(ax_perf, ylabel_str, 'Interpreter', 'latex');
    end
end

function drawDataDetail(ax_data, x_data, y_data, ratio_max_data, solver_names, profile_options, tolerance_latex)
    drawPerformanceDataProfiles(ax_data, x_data, y_data, solver_names, profile_options);
    % Set x-axis limits.
    set(ax_data, 'XLim', [0.0, 1.1 * ratio_max_data]);
    % Modify x-axis ticks labels of the data profiles.
    [ticks_data, tickLabels_data] = dataTicks(1.1 * ratio_max_data, profile_options.(ProfileOptionKey.SEMILOGX.value));
    set(ax_data, 'XTick', ticks_data, 'XTickLabel', tickLabels_data, 'TickLabelInterpreter', 'latex');
    % Set x-axis labels.
    if ~isempty(profile_options.(ProfileOptionKey.XLABEL_DATA_PROFILE.value))
        xlabel_str = profile_options.(ProfileOptionKey.XLABEL_DATA_PROFILE.value);
        xlabel(ax_data, xlabel_str, 'Interpreter', 'latex');
    end
    % Set y-axis labels.
    if ~isempty(profile_options.(ProfileOptionKey.YLABEL_DATA_PROFILE.value))
        ylabel_str = sprintf(profile_options.(ProfileOptionKey.YLABEL_DATA_PROFILE.value), tolerance_latex);
        ylabel(ax_data, ylabel_str, 'Interpreter', 'latex');
    end
end

function drawLogRatioDetail(ax_log_ratio, x_log_ratio, y_log_ratio, ratio_max_log_ratio, n_solvers_fail, solver_names, profile_options, tolerance_latex)
    drawLogRatioProfiles(ax_log_ratio, x_log_ratio, y_log_ratio, ratio_max_log_ratio, n_solvers_fail, solver_names, profile_options);
    % Set x-axis labels.
    if ~isempty(profile_options.(ProfileOptionKey.XLABEL_LOG_RATIO_PROFILE.value)) 
        xlabel_str = profile_options.(ProfileOptionKey.XLABEL_LOG_RATIO_PROFILE.value);
        xlabel(ax_log_ratio, xlabel_str, 'Interpreter', 'latex');
    end
    % Set y-axis labels.
    if ~isempty(profile_options.(ProfileOptionKey.YLABEL_LOG_RATIO_PROFILE.value))
        ylabel_str = sprintf(profile_options.(ProfileOptionKey.YLABEL_LOG_RATIO_PROFILE.value), tolerance_latex);
        ylabel(ax_log_ratio, ylabel_str, 'Interpreter', 'latex');
    end
end

function [ticks, tickLabels] = perfTicks(ratio_cut_perf, is_semilogx)

    if is_semilogx
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
    else
        if ratio_cut_perf >= 5
            max_power = floor(ratio_cut_perf);
            ticks = linspace(1, max_power, 5);
            ticks(2:end-1) = round(ticks(2:end-1));
            ticks = unique(ticks, 'stable');
        elseif ratio_cut_perf >= 2
            max_power = floor(ratio_cut_perf);
            ticks = (1:1:max_power);
        elseif ratio_cut_perf >= 1 + 1e-3
            ticks = [1 ratio_cut_perf];
        else
            ticks = [1];
        end
        tickLabels = arrayfun(@(x) num2str(x), ticks, 'UniformOutput', false);
    end
end

function [ticks, tickLabels] = dataTicks(ratio_cut_data, is_semilogx)

    if is_semilogx
        if ratio_cut_data >= 5
            max_power = floor(ratio_cut_data);
            ticks = linspace(1, max_power, 5);
            ticks = [0 ticks];
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
    else
        if ratio_cut_data >= 5
            max_power = floor(ratio_cut_data);
            ticks = linspace(0, max_power, 5);
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
        tickLabels = arrayfun(@(x) num2str(x), ticks, 'UniformOutput', false);
    end
end