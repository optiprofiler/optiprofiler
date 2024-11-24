function drawFunMaxcvMeritHist(ax, y, labels, is_cum, problem_n, y_shift, n_eval)
%DRAWFUNMAXCVMERITHIST draws figures of histories of function values, maximum constraint violation, or merit function values.

    set(ax, 'DefaultAxesColorOrder', [0, 0.4470, 0.7410; 0.8500, 0.3250, 0.0980; 0.9290, 0.6940, 0.1250; 0.4940, 0.1840, 0.5560; 0.4660, 0.6740, 0.1880; 0.3010, 0.7450, 0.9330; 0.6350, 0.0780, 0.1840]);
    set(ax, 'DefaultAxesLineStyleOrder', {'-', '--', ':', '-.'});

    % Shift the y-axis.
    y = y + y_shift;

    n_solvers = size(y, 1);
    n_runs = size(y, 2);
    y_mean = squeeze(mean(y, 2));

    y_lower = squeeze(min(y, [], 2));
    y_upper = squeeze(max(y, [], 2));

    if is_cum
        y_mean = cummin(y_mean, 2);
        y_lower = cummin(y_lower, 2);
        y_upper = cummin(y_upper, 2);
    end

    for i_solver = 1:n_solvers
        % Truncate the histories according to the function evaluations of each solver.
        length = max(n_eval(i_solver,:));
        length = min(length, size(y(i_solver,:,:), 3));
        nextColor = ax.ColorOrder(mod(ax.ColorOrderIndex-1, size(ax.ColorOrder, 1)) + 1, :);
        x = (1:length) / (problem_n + 1);
        plot(ax, x, y_mean(i_solver, 1:length), 'DisplayName', labels{i_solver});
        hold(ax, 'on');
        if n_runs > 1
            fill(ax, [x, fliplr(x)], [y_lower(i_solver, 1:length), fliplr(y_upper(i_solver, 1:length))], nextColor, 'FaceAlpha', 0.2, 'EdgeAlpha', 0, 'HandleVisibility', 'off');
        end
    end

    % When the function values are not all zero, use log scale for the y-axis.
    if any(y_mean(i_solver, :))
        set(ax, 'YScale', 'log');
    end

    box(ax, 'on');
    legend(ax, 'Location', 'northeast');
    hold(ax, 'off');
    set(ax, 'XLim', [min(x), max(x)]);
    xlabel(ax, 'Number of simplex gradients', 'Interpreter', 'latex');

end