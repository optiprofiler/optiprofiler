function drawFunMaxcvMeritHist(ax, y, solver_names, is_cum, problem_n, y_shift, n_eval, profile_options)
%DRAWFUNMAXCVMERITHIST draws figures of histories of function values, maximum constraint violation, or merit function values.

    line_colors = profile_options.(ProfileOptionKey.LINE_COLORS.value);
    line_styles = profile_options.(ProfileOptionKey.LINE_STYLES.value);
    line_widths = profile_options.(ProfileOptionKey.LINE_WIDTHS.value);

    % Shift the y-axis.
    y = y + y_shift;

    n_solvers = size(y, 1);
    n_runs = size(y, 2);
    y_mean = squeeze(mean(y, 2));

    if strcmp(profile_options.(ProfileOptionKey.RANGE_TYPE.value), 'minmax')
        y_lower = squeeze(min(y, [], 2));
        y_upper = squeeze(max(y, [], 2));
    else
        y_std = squeeze(std(y, 0, 2));
        y_lower = y_mean - y_std;
        y_upper = y_mean + y_std;
    end

    if is_cum
        y_mean = cummin(y_mean, 2);
        y_lower = cummin(y_lower, 2);
        y_upper = cummin(y_upper, 2);
    end

    xl_lim = 1 / (problem_n + 1);
    xr_lim = 1 / (problem_n + 1);

    hold(ax, 'on');
    for i_solver = 1:n_solvers
        % Truncate the histories according to the function evaluations of each solver.
        i_eval = max(max(n_eval(i_solver,:)), 1);
        i_eval = min(i_eval, size(y(i_solver,:,:), 3));
        x = (1:i_eval) / (problem_n + 1);
        xr_lim = max(xr_lim, x(end));

        color = line_colors(mod(i_solver - 1, size(line_colors, 1)) + 1, :);
        line_style = line_styles{mod(i_solver - 1, size(line_styles, 2)) + 1};
        line_width = line_widths(mod(i_solver - 1, length(line_widths)) + 1);

        if n_runs > 1 && i_eval > 1
            fill(ax, [x, fliplr(x)], [y_lower(i_solver, 1:i_eval), fliplr(y_upper(i_solver, 1:i_eval))], color, 'FaceAlpha', 0.2, 'EdgeAlpha', 0, 'HandleVisibility', 'off');
        end
        if i_eval == 1
            plot(ax, x, y_mean(i_solver, 1:i_eval), 'o', 'Color', color, 'DisplayName', solver_names{i_solver});
        else
            plot(ax, x, y_mean(i_solver, 1:i_eval), line_style, 'Color', color, 'LineWidth', line_width, 'DisplayName', solver_names{i_solver});
        end
    end

    % When the function values are not all zero and there is at least some change in the function values, , use log scale for the y-axis.
    if any(y_mean(i_solver, :)) && any(diff(y_mean(i_solver, :)))
        set(ax, 'YScale', 'log');
    end

    box(ax, 'on');
    legend(ax, 'Location', 'northeast');
    hold(ax, 'off');
    set(ax, 'XLim', [xl_lim, xr_lim]);
    xlabel(ax, 'Number of simplex gradients', 'Interpreter', 'latex');

end