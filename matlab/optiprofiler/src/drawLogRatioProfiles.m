function drawLogRatioProfiles(ax, x, y, ratio_max, n_solvers_equal, solver_names, profile_options)
%DRAWLOGRATIOPROFILES draws log-ratio profiles.

    n_problems = length(x);
    n_below = sum(y < 0) - n_solvers_equal;
    n_above = sum(y > 0) - n_solvers_equal;

    hold(ax, 'on');
    % Draw the log-ratio profiles.
    if n_solvers_equal > 0
        bar(ax, x(1:n_solvers_equal), y(1:n_solvers_equal), 'FaceColor', profile_options.(ProfileOptionKey.BAR_COLORS.value)(1, :), 'LineWidth', 1.5, 'LineStyle', 'None', 'FaceAlpha', 0.5);
        bar(ax, x(end - n_solvers_equal + 1:end), y(end - n_solvers_equal + 1:end), 'FaceColor', profile_options.(ProfileOptionKey.BAR_COLORS.value)(mod(1, size(profile_options.(ProfileOptionKey.BAR_COLORS.value), 1)) + 1, :), 'LineWidth', 1.5, 'LineStyle', 'None', 'FaceAlpha', 0.5);
    end
    if n_below > 0
        bar(ax, x(n_solvers_equal + 1:n_solvers_equal + n_below), y(n_solvers_equal + 1:n_solvers_equal + n_below), 'FaceColor', profile_options.(ProfileOptionKey.BAR_COLORS.value)(1, :), 'LineWidth', 1.5, 'LineStyle', 'None');
    end
    if n_above > 0
        bar(ax, x(end - n_above + 1:end - n_solvers_equal), y(end - n_above + 1:end - n_solvers_equal), 'FaceColor', profile_options.(ProfileOptionKey.BAR_COLORS.value)(mod(1, size(profile_options.(ProfileOptionKey.BAR_COLORS.value), 1)) + 1, :), 'LineWidth', 1.5, 'LineStyle', 'None');
    end
    text(ax, (n_problems + 1) / 2, -ratio_max, solver_names{1}, 'HorizontalAlignment', 'center', 'VerticalAlignment', 'bottom', 'FontSize', 24);
    text(ax, (n_problems + 1) / 2, ratio_max, solver_names{2}, 'HorizontalAlignment', 'center', 'VerticalAlignment', 'top', 'FontSize', 24);
    xticks(ax, []);
    xlim(ax, [0.5, n_problems + 0.5]);
    yyaxis(ax, 'right');
    set(ax, 'YTickLabel', []);
    yyaxis(ax, 'left');
    ylim(ax, [-1.1 * ratio_max, 1.1 * ratio_max]);
    linkprop([ax.XAxis; ax.YAxis],'color');
    linkprop([ax.YAxis(1), ax.YAxis(2)],{'Limits','TickValues'});
    box(ax,'on');
    hold(ax, 'off');
end