function [fig, ax] = drawProfile(x, y, perf_or_data, ratio_max, labels, x_label, y_label)
    n_solvers = size(x, 2);
    fig = figure('visible', 'off');
    ax = axes;
    ax.LineStyleCyclingMethod = 'withcolor';

    hold on;

    switch perf_or_data
        case 'perf'
            for i_solver = 1:n_solvers
                [x_stairs, y_stairs] = stairs(x(:, i_solver), y(:, i_solver));
                x_stairs = log2([x_stairs(1); x_stairs; 2.0 * ratio_max]);
                y_stairs = [0.0; y_stairs; y_stairs(end)];
                plot(ax, x_stairs, y_stairs, 'DisplayName', labels{i_solver});
            end
            set(ax, 'YTick', 0:0.1:1, 'YTickLabel', {'0', '', '0.2', '', '0.4', '', '0.6', '', '0.8', '', '1'}, ...
                'XLim', [0.0, 1.1 * log2(ratio_max)], 'YLim', [0.0, 1.0], 'TickLabelInterpreter', 'latex');
        case 'data'
            for i_solver = 1:n_solvers
                [x_stairs, y_stairs] = stairs(x(:, i_solver), y(:, i_solver));
                x_stairs = [x_stairs(1); x_stairs; 2.0 * ratio_max];
                y_stairs = [0.0; y_stairs; y_stairs(end)];
                plot(ax, x_stairs, y_stairs, 'DisplayName', labels{i_solver});
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