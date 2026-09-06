function drawFunMaxcvMeritHist(ax, y, solver_names, is_cum, problem_n, y_shift, n_eval, profile_options, show_xlabel)
%DRAWFUNMAXCVMERITHIST draws figures of histories of function values, maximum constraint violation, or merit function values.

    if nargin < 9
        show_xlabel = true;
    end

    line_colors = profile_options.(ProfileOptionKey.LINE_COLORS.value);
    line_styles = profile_options.(ProfileOptionKey.LINE_STYLES.value);
    line_widths = profile_options.(ProfileOptionKey.LINE_WIDTHS.value);

    [x_indices, y_values_m, y_values_l, y_values_u, n_runs] = ...
        prepareHistoryPlotData(y, is_cum, y_shift, n_eval, profile_options);
    n_solvers = size(y, 1);

    xl_lim = 1 / (problem_n + 1);
    xr_lim = 1 / (problem_n + 1);
    is_log_scale = false;

    hold(ax, 'on');
    for i_solver = 1:n_solvers
        % Truncate the histories according to the function evaluations of each solver.
        i_x = x_indices{i_solver};
        i_y_mean = y_values_m{i_solver};
        i_y_lower = y_values_l{i_solver};
        i_y_upper = y_values_u{i_solver};
        i_eval = length(i_x);
        if i_eval == 0
            continue;
        end
        x = i_x / (problem_n + 1);
        xr_lim = max(xr_lim, x(end));

        color = line_colors(mod(i_solver - 1, size(line_colors, 1)) + 1, :);
        line_style = line_styles{mod(i_solver - 1, size(line_styles, 2)) + 1};
        line_width = line_widths(mod(i_solver - 1, length(line_widths)) + 1);

        if n_runs > 1 && i_eval > 1
            fill(ax, [x, fliplr(x)], [i_y_lower, fliplr(i_y_upper)], color, 'FaceAlpha', 0.2, 'EdgeAlpha', 0, 'HandleVisibility', 'off');
        end
        if i_eval == 1
            plot(ax, x, i_y_mean, 'o', 'Color', color, 'DisplayName', solver_names{i_solver});
        elseif i_eval > 1
            plot(ax, x, i_y_mean, line_style, 'Color', color, 'LineWidth', line_width, 'DisplayName', solver_names{i_solver});
        end
        if any(i_y_mean) && any(diff(i_y_mean))
            is_log_scale = true;
        end
    end

    % When the function values are not all zero and there is at least some change in the function values, use log scale for the y-axis.
    if is_log_scale
        set(ax, 'YScale', 'log');
        shown = [y_values_m{:}, y_values_l{:}, y_values_u{:}];
        positive = shown(shown > 0);
        if ~isempty(positive)
            log_low = log10(min(positive));
            log_high = log10(max(positive));
            margin = max(0.05 * (log_high - log_low), 0.05);
            % Explicit representable limits avoid overflow/underflow in
            % automatic logarithmic margins. No floor is applied to data.
            smallest = realmin * eps;
            low = max(smallest, 10 ^ max(log10(smallest), log_low - margin));
            high = 10 ^ min(102, log_high + margin);
            set(ax, 'YLim', [low, high]);
            if log_low < -250 || log_high - log_low > 200
                first = ceil(log10(low));
                last = floor(log10(high));
                step = max(1, ceil((last - first) / 8));
                set(ax, 'YTick', 10 .^ (first:step:last), 'YMinorTick', 'off');
            end
        end
    end

    box(ax, 'on');
    hold(ax, 'off');
    if abs(xl_lim - xr_lim) < eps
        xr_lim = xl_lim + eps;
    end
    set(ax, 'FontSize', 9);
    set(ax, 'XLim', [xl_lim, xr_lim]);
    if show_xlabel
        x_label = xlabel(ax, profile_options.(ProfileOptionKey.XLABEL_DATA_PROFILE.value), 'Interpreter', 'latex', 'FontSize', 10);
        set(x_label, 'Units', 'normalized', 'Position', [0.5, -0.10, 0], 'VerticalAlignment', 'top');
    else
        xlabel(ax, '');
    end
    placeSolverLegend(ax, n_solvers, 'northeast');

end
