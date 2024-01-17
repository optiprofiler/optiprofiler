function [fig, ax] = drawProfile(x, y, perf_or_data, ratio_max, labels, x_label, y_label)
    n_solvers = size(x, 2);
    y_mean = squeeze(mean(y, 3));
    y_min = squeeze(min(y, [], 3));
    y_max = squeeze(max(y, [], 3));
    fig = figure('visible', 'off');
    ax = axes;

    ax.LineStyleCyclingMethod = 'withcolor';

    hold on;

    switch perf_or_data
        case 'perf'
            for i_solver = 1:n_solvers
                [x_stairs, y_mean_stairs] = stairs(x(:, i_solver), y_mean(:, i_solver));
                [~, y_min_stairs] = stairs(x(:, i_solver), y_min(:, i_solver));
                [~, y_max_stairs] = stairs(x(:, i_solver), y_max(:, i_solver));
                x_stairs = log2([x_stairs(1); x_stairs; 2.0 * ratio_max]);
                y_mean_stairs = [0.0; y_mean_stairs; y_mean_stairs(end)];
                y_min_stairs = [0.0; y_min_stairs; y_min_stairs(end)];
                y_max_stairs = [0.0; y_max_stairs; y_max_stairs(end)];

                % Get the color MATLAB will use for the next plot command in the axes 'ax'.
                nextColor = ax.ColorOrder(mod(ax.ColorOrderIndex-1, size(ax.ColorOrder, 1)) + 1, :);
                plot(ax, x_stairs, y_mean_stairs, 'DisplayName', labels{i_solver});
                fill(ax, [x_stairs; flipud(x_stairs)], [y_min_stairs; flipud(y_max_stairs)], ...
                    nextColor, 'FaceAlpha', 0.2, 'EdgeAlpha', 0, 'HandleVisibility', 'off');
            end
            set(ax, 'YTick', 0:0.1:1, 'YTickLabel', {'0', '', '0.2', '', '0.4', '', '0.6', '', '0.8', '', '1'}, ...
                'XLim', [0.0, 1.1 * log2(ratio_max)], 'YLim', [0.0, 1.0], 'TickLabelInterpreter', 'latex');
        case 'data'
            for i_solver = 1:n_solvers
                [x_stairs, y_mean_stairs] = stairs(x(:, i_solver), y_mean(:, i_solver));
                [~, y_min_stairs] = stairs(x(:, i_solver), y_min(:, i_solver));
                [~, y_max_stairs] = stairs(x(:, i_solver), y_max(:, i_solver));
                x_stairs = [x_stairs(1); x_stairs; 2.0 * ratio_max];
                y_mean_stairs = [0.0; y_mean_stairs; y_mean_stairs(end)];
                y_min_stairs = [0.0; y_min_stairs; y_min_stairs(end)];
                y_max_stairs = [0.0; y_max_stairs; y_max_stairs(end)];

                % Get the color MATLAB will use for the next plot command in the axes 'ax'.
                nextColor = ax.ColorOrder(mod(ax.ColorOrderIndex-1, size(ax.ColorOrder, 1)) + 1, :);
                plot(ax, x_stairs, y_mean_stairs, 'DisplayName', labels{i_solver});
                fill(ax, [x_stairs; flipud(x_stairs)], [y_min_stairs; flipud(y_max_stairs)], ...
                    nextColor, 'FaceAlpha', 0.2, 'EdgeAlpha', 0, 'HandleVisibility', 'off');
            end         
            set(ax, 'YTick', 0:0.1:1, 'YTickLabel', {'0', '', '0.2', '', '0.4', '', '0.6', '', '0.8', '', '1'}, ...
                'XLim', [0.0, 1.1 * ratio_max], 'YLim', [0.0, 1.0], 'TickLabelInterpreter', 'latex');
        otherwise
            error("Unknown type of data to plot.");
    end
    
    

    xlabel(ax, x_label, 'Interpreter', 'latex');
    ylabel(ax, y_label, 'Interpreter', 'latex');
    
    legend(ax, 'Location', 'southeast');
    hold off;
end