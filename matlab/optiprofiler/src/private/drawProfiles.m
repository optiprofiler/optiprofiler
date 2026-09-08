function [fig_perf, fig_data, fig_log_ratio, curves, graphics_ok, presentation] = drawProfiles(work, problem_dimensions, solver_names, tolerance_latex, cell_axs_summary, is_summary, is_perf, is_data, is_log_ratio, profile_options, curves)
%DRAWPROFILES draws the performance, data, and log-ratio profiles.

    solver_names = cellfun(@escapeLatexText, solver_names, 'UniformOutput', false);
    n_solvers = size(work, 2);

    % Numerical curves are also the scoring input. Compute them independently
    % of graphics so score_only and the no-graphics SVG fallback are identical.
    fig_perf = []; fig_data = []; fig_log_ratio = []; graphics_ok = true;
    [x_perf, y_perf, ratio_max_perf, x_data, y_data, ratio_max_data, curves] = getExtendedPerformancesDataProfileAxes(work, problem_dimensions, profile_options, curves);
    if n_solvers == 2
        if nargout > 5
            [x_log_ratio, y_log_ratio, ratio_max_log_ratio, n_solvers_fail, curves, bar_sources] = getLogRatioProfileAxes(work, curves);
        else
            [x_log_ratio, y_log_ratio, ratio_max_log_ratio, n_solvers_fail, curves] = getLogRatioProfileAxes(work, curves);
        end
    end
    presentation = {};
    if nargout > 5
      try
        % Preparing a report must not allocate a figure. These are precisely
        % the shared arrays consumed by the native renderer below.
        presentation = [profilePresentation(x_perf,y_perf,ratio_max_perf,'performance',profile_options), ...
            profilePresentation(x_data,y_data,ratio_max_data,'data',profile_options)];
        if n_solvers == 2
            % Same record vocabulary as Python prepare_log_ratio_plot_data;
            % MATLAB additionally retains bar identity (problem_mapping).
            opacity = ones(size(y_log_ratio)); opacity(bar_sources.solver1_failed & bar_sources.solver2_failed) = 0.5;
            for field = fieldnames(bar_sources)'
                values = bar_sources.(field{1});
                bar_sources.(field{1}) = num2cell(values(:)');
            end
            n_bars = numel(y_log_ratio);
            n_below = sum(y_log_ratio < 0) - n_solvers_fail;
            n_above = sum(y_log_ratio > 0) - n_solvers_fail;
            % The bar() calls of drawLogRatioProfiles in drawing order.
            groups = {};
            if n_solvers_fail > 0
                groups{end+1} = struct('start_index', 1, 'end_index', n_solvers_fail, 'solver_index', 1, 'opacity', 0.5);
                groups{end+1} = struct('start_index', n_bars - n_solvers_fail + 1, 'end_index', n_bars, 'solver_index', 2, 'opacity', 0.5);
            end
            if n_below > 0
                groups{end+1} = struct('start_index', n_solvers_fail + 1, 'end_index', n_solvers_fail + n_below, 'solver_index', 1, 'opacity', 1);
            end
            if n_above > 0
                groups{end+1} = struct('start_index', n_bars - n_solvers_fail - n_above + 1, 'end_index', n_bars - n_solvers_fail, 'solver_index', 2, 'opacity', 1);
            end
            presentation{end+1} = struct('kind', 'log_ratio', 'status', 'numeric_prepared', ...
                'fidelity', 'exact_rendered_data', 'series', {{struct('geometry', 'bar', 'x', {num2cell(x_log_ratio(:)')}, ...
                    'y', {num2cell(y_log_ratio(:)')}, 'visible', {num2cell(y_log_ratio(:)'~=0)}, ...
                    'opacity', {num2cell(opacity(:)')}, 'solver_index_by_sign', struct('negative', 1, 'positive', 2))}}, ...
                'bar_groups', {groups}, 'ratio_max', ratio_max_log_ratio, ...
                'failure_placeholder', 1.1*ratio_max_log_ratio, ...
                'x_transform', 'sorted_bar_position', 'y_transform', 'log2(work_solver_1/work_solver_2)', ...
                'x_limits', {{0.5, n_bars + 0.5}}, 'y_limits', {{-1.1*ratio_max_log_ratio, 1.1*ratio_max_log_ratio}}, ...
                'both_failed_pairs', n_solvers_fail, 'tie_pairs', sum(y_log_ratio == 0), ...
                'problem_mapping', 'bar_sources', 'bar_sources', bar_sources);
        end
      catch cause
        % A secondary numeric-presentation failure is report evidence, not
        % permission to abort the original scoring/rendering execution.
        presentation = {struct('preparation_error', cause.identifier)};
      end
    end
    if profile_options.(ProfileOptionKey.SCORE_ONLY.value)
        return;
    end

    try
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % Create the individual figures.
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        fig_perf = figure('Units', 'pixels', 'Position', profileFigurePosition(), 'visible', 'off');
        t_perf = tiledlayout(fig_perf, 1, 1, 'Padding', 'compact', 'TileSpacing', 'compact');
        ax_perf = nexttile(t_perf);
        fig_data = figure('Units', 'pixels', 'Position', profileFigurePosition(), 'visible', 'off');
        t_data = tiledlayout(fig_data, 1, 1, 'Padding', 'compact', 'TileSpacing', 'compact');
        ax_data = nexttile(t_data);
        if n_solvers == 2
            fig_log_ratio = figure('Units', 'pixels', 'Position', profileFigurePosition(), 'visible', 'off');
            t_log_ratio = tiledlayout(fig_log_ratio, 1, 1, 'Padding', 'compact', 'TileSpacing', 'compact');
            ax_log_ratio = nexttile(t_log_ratio);
        else
            fig_log_ratio = [];
            ax_log_ratio = [];
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
    catch cause
        % Only graphics are guarded: curve/scoring calculations above must
        % retain their original failure semantics. The caller emits SVG from
        % the completed curves and does not retry or recompute solver data.
        graphics_ok = false;
        for candidate = {fig_perf, fig_data, fig_log_ratio}
            if ~isempty(candidate{1}) && isgraphics(candidate{1}), close(candidate{1}); end
        end
        fig_perf = []; fig_data = []; fig_log_ratio = [];
        printOptiProfilerMessage('WARNING', sprintf('Native profile rendering failed; using SVG fallback: %s',cause.message));
    end
end

function panels = profilePresentation(x,y,maximum,kind,options)
    % Same record vocabulary as the Python _plot_sink data in draw_profiles;
    % x_transform labels and limits are pinned by plot_data.schema.json.
    [xs,means,lower,upper,n_runs] = prepareProfilePlotData(x,y,options);
    transform = 'work/best_work';
    if strcmp(kind,'data'), transform = 'work/(dimension+1)'; end
    if options.semilogx
        transform = 'log2(work/best_work)';
        if strcmp(kind,'data'), transform = 'log2(1+work/(dimension+1))'; end
    end
    % The limits actually given to set(ax, 'XLim', ...) by drawPerfDetail /
    % drawDataDetail (MATLAB tolerates an infinite ratio_max there).
    x_low = 0;
    if strcmp(kind,'performance') && ~options.semilogx, x_low = 1; end
    series = cell(1,numel(xs));
    for solver = 1:numel(xs)
        series{solver} = struct('solver_index', solver, ...
            'x', {num2cell(xs{solver}(:)')}, 'mean', {num2cell(means{solver}(:)')}, ...
            'lower', {num2cell(lower{solver}(:)')}, 'upper', {num2cell(upper{solver}(:)')}, ...
            'geometry', 'step', 'band_visible', n_runs>1);
    end
    panels = {struct('kind', kind, 'status', 'numeric_prepared', ...
        'fidelity', 'exact_rendered_data', 'series', {series}, ...
        'x_transform', transform, 'ratio_max', maximum, 'failure_placeholder', 1.1*maximum, ...
        'x_limits', {{x_low, 1.1*maximum}}, 'y_limits', {{0, 1}}, ...
        'errorbar_type', options.errorbar_type, 'std_ddof', double(n_runs>1), 'n_runs', n_runs)};
end

function drawPerfDetail(ax_perf, x_perf, y_perf, ratio_max_perf, solver_names, profile_options, tolerance_latex)
    drawPerformanceDataProfiles(ax_perf, x_perf, y_perf, solver_names, profile_options);
    set(ax_perf, 'FontSize', 10);
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
        xlabel(ax_perf, xlabel_str, 'Interpreter', 'latex', 'FontSize', 11);
    end
    % Set y-axis labels.
    if ~isempty(profile_options.(ProfileOptionKey.YLABEL_PERFORMANCE_PROFILE.value))
        ylabel_str = sprintf(profile_options.(ProfileOptionKey.YLABEL_PERFORMANCE_PROFILE.value), tolerance_latex);
        ylabel(ax_perf, ylabel_str, 'Interpreter', 'latex', 'FontSize', 11);
    end
    placeSolverLegend(ax_perf, size(x_perf, 2), 'southeast');
end

function drawDataDetail(ax_data, x_data, y_data, ratio_max_data, solver_names, profile_options, tolerance_latex)
    drawPerformanceDataProfiles(ax_data, x_data, y_data, solver_names, profile_options);
    set(ax_data, 'FontSize', 10);
    % Set x-axis limits.
    set(ax_data, 'XLim', [0.0, 1.1 * ratio_max_data]);
    % Modify x-axis ticks labels of the data profiles.
    [ticks_data, tickLabels_data] = dataTicks(1.1 * ratio_max_data, profile_options.(ProfileOptionKey.SEMILOGX.value));
    set(ax_data, 'XTick', ticks_data, 'XTickLabel', tickLabels_data, 'TickLabelInterpreter', 'latex');
    % Set x-axis labels.
    if ~isempty(profile_options.(ProfileOptionKey.XLABEL_DATA_PROFILE.value))
        xlabel_str = profile_options.(ProfileOptionKey.XLABEL_DATA_PROFILE.value);
        xlabel(ax_data, xlabel_str, 'Interpreter', 'latex', 'FontSize', 11);
    end
    % Set y-axis labels.
    if ~isempty(profile_options.(ProfileOptionKey.YLABEL_DATA_PROFILE.value))
        ylabel_str = sprintf(profile_options.(ProfileOptionKey.YLABEL_DATA_PROFILE.value), tolerance_latex);
        ylabel(ax_data, ylabel_str, 'Interpreter', 'latex', 'FontSize', 11);
    end
    placeSolverLegend(ax_data, size(x_data, 2), 'southeast');
end

function drawLogRatioDetail(ax_log_ratio, x_log_ratio, y_log_ratio, ratio_max_log_ratio, n_solvers_fail, solver_names, profile_options, tolerance_latex)
    drawLogRatioProfiles(ax_log_ratio, x_log_ratio, y_log_ratio, ratio_max_log_ratio, n_solvers_fail, solver_names, profile_options);
    set(ax_log_ratio, 'FontSize', 10);
    % Set x-axis labels.
    if ~isempty(profile_options.(ProfileOptionKey.XLABEL_LOG_RATIO_PROFILE.value))
        xlabel_str = profile_options.(ProfileOptionKey.XLABEL_LOG_RATIO_PROFILE.value);
        xlabel(ax_log_ratio, xlabel_str, 'Interpreter', 'latex', 'FontSize', 11);
    end
    % Set y-axis labels.
    if ~isempty(profile_options.(ProfileOptionKey.YLABEL_LOG_RATIO_PROFILE.value))
        ylabel_str = sprintf(profile_options.(ProfileOptionKey.YLABEL_LOG_RATIO_PROFILE.value), tolerance_latex);
        y_label = ylabel(ax_log_ratio, ylabel_str, 'Interpreter', 'latex', 'FontSize', 11);
        set(y_label, 'Units', 'normalized');
        label_position = y_label.Position;
        label_position(1) = -0.05;
        y_label.Position = label_position;
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
        % Data profiles are plotted at z = log2(1 + alpha), where
        % alpha = n_eval / (n + 1). Choose ticks in the original alpha
        % units and then map them to z so that labels remain literal
        % simplex-gradient budgets (0, 1, 2, 4, ...).
        ratio_cut_alpha = max(0, 2 ^ ratio_cut_data - 1);
        if ratio_cut_alpha >= 1
            max_power = floor(log2(ratio_cut_alpha));
            if max_power + 1 <= 5
                powers = 0:1:max_power;
            else
                powers = unique([0, floor(linspace(1, max_power, 4))], 'stable');
            end
            alpha_ticks = [0, 2 .^ powers];
            ticks = log2(1 + alpha_ticks);
            tickLabels = arrayfun(@(x) num2str(x), alpha_ticks, 'UniformOutput', false);
        elseif ratio_cut_alpha >= 1e-1
            ticks = [0 ratio_cut_data];
            tickLabels = {'0', stripTrailingZeros(sprintf('%.2f', ratio_cut_alpha))};
        else
            ticks = [0];
            tickLabels = {'0'};
        end
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

function label = stripTrailingZeros(label)
    label = regexprep(label, '0+$', '');
    label = regexprep(label, '\.$', '');
end
