function drawPerformanceDataProfiles(ax, x, y, solver_names, profile_options)
%DRAWPERFORMANCEDATAPROFILES draws performance profiles and data profiles.

    line_color_order = profile_options.line_color_order;
    line_style_order = profile_options.line_style_order;


    n_solvers = size(x, 2);
    n_runs = size(y, 3);
    y_mean = squeeze(mean(y, 3));

    switch profile_options.(ProfileOptionKey.RANGE_TYPE.value)
        case 'minmax'
            y_lower = squeeze(min(y, [], 3));
            y_upper = squeeze(max(y, [], 3));
        case 'meanstd'
            y_std = squeeze(std(y, [], 3));
            y_lower = max(y_mean - y_std, 0);
            y_upper = min(y_mean + y_std, 1);
        otherwise
            error("Unknown range type.");
    end

    hold(ax, 'on');

    for i_solver = 1:n_solvers
        [x_stairs, y_mean_stairs] = stairs(x(:, i_solver), y_mean(:, i_solver));
        [~, y_lower_stairs] = stairs(x(:, i_solver), y_lower(:, i_solver));
        [~, y_upper_stairs] = stairs(x(:, i_solver), y_upper(:, i_solver));

        % Get the color and the line style MATLAB will use for the next plot command in the axes 'ax'.
        color = line_color_order(mod(i_solver - 1, size(line_color_order, 1)) + 1, :);
        line_style = line_style_order{mod(i_solver - 1, size(line_style_order, 2)) + 1};

        % We first plot the shaded area and then the mean line to ensure that the shaded area is behind the mean line in case the installed MATLAB version does not support transparency.
        if n_runs > 1
            fill(ax, [x_stairs; flipud(x_stairs)], [y_lower_stairs; flipud(y_upper_stairs)], color, 'FaceAlpha', 0.2, 'EdgeAlpha', 0, 'HandleVisibility', 'off');
        end
        plot(ax, x_stairs, y_mean_stairs, line_style, 'Color', color, 'DisplayName', solver_names{i_solver});
    end

    set(ax, 'YLim', [0.0, 1.0]);
    set(ax, 'YTick', 0:0.2:1);
    yyaxis(ax, 'right');
    set(ax, 'YMinorTick','on');
    set(get(ax, 'YAxis'), 'MinorTickValues', [0.1 0.3 0.5 0.7 0.9]);
    set(ax, 'YTickLabel', []);
    yyaxis(ax, 'left');
    set(ax, 'YMinorTick','on');
    set(get(ax, 'YAxis'), 'MinorTickValues', [0.1 0.3 0.5 0.7 0.9]);
    linkprop([ax.XAxis; ax.YAxis],'color');
    linkprop([ax.YAxis(1), ax.YAxis(2)],{'Limits','TickValues'});
    box(ax,'on');
    legend(ax, 'Location', 'southeast');
    hold(ax, 'off');
end