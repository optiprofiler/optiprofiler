function drawPerformanceDataProfiles(ax, x, y, labels)
%DRAWPERFORMANCEDATAPROFILES draws performance profiles and data profiles.

    ax.LineStyleCyclingMethod = 'withcolor';
    n_solvers = size(x, 2);
    n_runs = size(y, 3);
    y_mean = squeeze(mean(y, 3));
    y_min = squeeze(min(y, [], 3));
    y_max = squeeze(max(y, [], 3));
    hold(ax, 'on');

    for i_solver = 1:n_solvers
        [x_stairs, y_mean_stairs] = stairs(x(:, i_solver), y_mean(:, i_solver));
        [~, y_min_stairs] = stairs(x(:, i_solver), y_min(:, i_solver));
        [~, y_max_stairs] = stairs(x(:, i_solver), y_max(:, i_solver));

        % Get the color MATLAB will use for the next plot command in the axes 'ax'.
        nextColor = ax.ColorOrder(mod(ax.ColorOrderIndex-1, size(ax.ColorOrder, 1)) + 1, :);

        plot(ax, x_stairs, y_mean_stairs, 'DisplayName', labels{i_solver});
        if n_runs > 1
            fill(ax, [x_stairs; flipud(x_stairs)], [y_min_stairs; flipud(y_max_stairs)], ...
            nextColor, 'FaceAlpha', 0.2, 'EdgeAlpha', 0, 'HandleVisibility', 'off');
        end
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