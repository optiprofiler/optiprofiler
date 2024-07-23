function drawPerformanceDataProfiles(ax, x, y, labels, profile_options)
%DRAWPERFORMANCEDATAPROFILES draws performance profiles and data profiles.

    set(ax, 'DefaultAxesColorOrder', [0, 0.4470, 0.7410; 0.8500, 0.3250, 0.0980; 0.9290, 0.6940, 0.1250; 0.4940, 0.1840, 0.5560; 0.4660, 0.6740, 0.1880; 0.3010, 0.7450, 0.9330; 0.6350, 0.0780, 0.1840]);
    set(ax, 'DefaultAxesLineStyleOrder', {'-', '--', ':', '-.'});

    n_solvers = size(x, 2);
    n_runs = size(y, 3);
    y_mean = squeeze(mean(y, 3));

    switch profile_options.(ProfileOptionKey.RANGE_TYPE.value)
        case 'minmax'
            y_lower = squeeze(min(y, [], 3));
            y_upper = squeeze(max(y, [], 3));
        case 'meanstd'
            y_std = squeeze(std(y, [], 3));
            y_lower = max(y_mean - profile_options.(ProfileOptionKey.STD_FACTOR.value) * y_std, 0);
            y_upper = min(y_mean + profile_options.(ProfileOptionKey.STD_FACTOR.value) * y_std, 1);
        otherwise
            error("Unknown range type.");
    end

    hold(ax, 'on');

    for i_solver = 1:n_solvers
        [x_stairs, y_mean_stairs] = stairs(x(:, i_solver), y_mean(:, i_solver));
        [~, y_lower_stairs] = stairs(x(:, i_solver), y_lower(:, i_solver));
        [~, y_upper_stairs] = stairs(x(:, i_solver), y_upper(:, i_solver));

        % Get the color MATLAB will use for the next plot command in the axes 'ax'.
        nextColor = ax.ColorOrder(mod(ax.ColorOrderIndex-1, size(ax.ColorOrder, 1)) + 1, :);

        plot(ax, x_stairs, y_mean_stairs, 'DisplayName', labels{i_solver});
        if n_runs > 1
            fill(ax, [x_stairs; flipud(x_stairs)], [y_lower_stairs; flipud(y_upper_stairs)], ...
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