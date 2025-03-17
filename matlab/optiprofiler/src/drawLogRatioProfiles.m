function drawLogRatioProfiles(ax, x, y, ratio_max, solver_names, profile_options)
%DRAWLOGRATIOPROFILES draws log-ratio profiles.

    length_x = length(x);

    % Draw the log-ratio profiles.
    bar(ax, x(y < 0), y(y < 0), 'FaceColor', profile_options.(ProfileOptionKey.BAR_COLORS.value)(1, :), 'LineWidth', 1.5, 'LineStyle', 'None');
    hold(ax, 'on');
    % In case bar_colors only contains one color.
    bar(ax, x(y > 0), y(y > 0), 'FaceColor', profile_options.(ProfileOptionKey.BAR_COLORS.value)(mod(1, size(profile_options.(ProfileOptionKey.BAR_COLORS.value), 1)) + 1, :), 'LineWidth', 1.5, 'LineStyle', 'None');
    text(ax, (length_x + 1) / 2, -ratio_max, solver_names{1}, 'HorizontalAlignment', 'center', 'VerticalAlignment', 'bottom', 'FontSize', 24);
    text(ax, (length_x + 1) / 2, ratio_max, solver_names{2}, 'HorizontalAlignment', 'center', 'VerticalAlignment', 'top', 'FontSize', 24);
    xticks(ax, []);
    xlim(ax, [0.5, length_x + 0.5]);
    yyaxis(ax, 'right');
    set(ax, 'YTickLabel', []);
    yyaxis(ax, 'left');
    ylim(ax, [-1.1 * ratio_max, 1.1 * ratio_max]);
    linkprop([ax.XAxis; ax.YAxis],'color');
    linkprop([ax.YAxis(1), ax.YAxis(2)],{'Limits','TickValues'});
    box(ax,'on');
    hold(ax, 'off');
end