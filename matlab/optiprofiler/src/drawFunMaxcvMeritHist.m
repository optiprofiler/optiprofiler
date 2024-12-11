function drawFunMaxcvMeritHist(ax, y, labels, is_cum, problem_n, y_shift, n_eval, profile_options)
%DRAWFUNMAXCVMERITHIST draws figures of histories of function values, maximum constraint violation, or merit function values.

    set(ax, 'DefaultAxesColorOrder', [0, 0.4470, 0.7410; 0.8500, 0.3250, 0.0980; 0.9290, 0.6940, 0.1250; 0.4940, 0.1840, 0.5560; 0.4660, 0.6740, 0.1880; 0.3010, 0.7450, 0.9330; 0.6350, 0.0780, 0.1840]);
    set(ax, 'DefaultAxesLineStyleOrder', {'-', '--', ':', '-.'});

    % Shift the y-axis.
    y = y + y_shift;

    % Set a upper bound to ensure the y-axis is not too large.
    upper_bound = 1e100;

    n_solvers = size(y, 1);
    n_runs = size(y, 2);
    y_mean = squeeze(mean(y, 2));
    % We want to make sure y is not too large (even Inf).
    y_mean = min(y_mean, upper_bound + y_shift);

    if strcmp(profile_options.(ProfileOptionKey.RANGE_TYPE.value), 'minmax')
        y_lower = squeeze(min(y, [], 2));
        y_upper = squeeze(max(y, [], 2));
        % We want to make sure y is not too large (even Inf).
        y_lower = min(y_lower, upper_bound + y_shift);
        y_upper = min(y_upper, upper_bound + y_shift);
    else
        y_std = squeeze(std(y, 0, 2));
        y_lower = y_mean - y_std;
        y_upper = y_mean + y_std;
        % We want to make sure y is not too large (even Inf).
        y_lower = min(y_lower, upper_bound + y_shift);
        y_upper = min(y_upper, upper_bound + y_shift);
    end

    if is_cum
        y_mean = cummin(y_mean, 2);
        y_lower = cummin(y_lower, 2);
        y_upper = cummin(y_upper, 2);
    end

    xl_lim = 1 / (problem_n + 1);
    xr_lim = 1 / (problem_n + 1);
    for i_solver = 1:n_solvers
        % Truncate the histories according to the function evaluations of each solver.
        length = max(n_eval(i_solver,:));
        length = min(length, size(y(i_solver,:,:), 3));
        nextColor = ax.ColorOrder(mod(ax.ColorOrderIndex-1, size(ax.ColorOrder, 1)) + 1, :);
        x = (1:length) / (problem_n + 1);
        xr_lim = max(xr_lim, x(end));
        if length == 1
            plot(ax, x, y_mean(i_solver, 1:length), 'o', 'DisplayName', labels{i_solver});
        else
            plot(ax, x, y_mean(i_solver, 1:length), 'DisplayName', labels{i_solver});
        end
        hold(ax, 'on');
        if n_runs > 1 && length > 1
            fill(ax, [x, fliplr(x)], [y_lower(i_solver, 1:length), fliplr(y_upper(i_solver, 1:length))], nextColor, 'FaceAlpha', 0.2, 'EdgeAlpha', 0, 'HandleVisibility', 'off');
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