function y_shift = drawFunMaxcvMeritHist(ax, y, labels, ismaxcv, profile_options, problem_n)
%DRAWFUNMAXCVMERITHIST draws figures of histories of function values, maximum constraint violation, or merit function values.

    set(ax, 'DefaultAxesColorOrder', [0, 0.4470, 0.7410; 0.8500, 0.3250, 0.0980; 0.9290, 0.6940, 0.1250; 0.4940, 0.1840, 0.5560; 0.4660, 0.6740, 0.1880; 0.3010, 0.7450, 0.9330; 0.6350, 0.0780, 0.1840]);
    set(ax, 'DefaultAxesLineStyleOrder', {'-', '--', ':', '-.'});

    % Find out the smallest and largest value of y.
    y_max = max(y(:));
    y_min = min(y(:));

    % Decide the shift of the y-axis to make the shifted y-axis be positive.
    y_shift = 0;
    if y_min <= 1e-6 * (y_max - y_min)
        y_shift = 1e-6 * (y_max - y_min) - y_min;
    end
    
    % Shift the y-axis.
    y = y + y_shift;

    n_solvers = size(y, 1);
    n_runs = size(y, 2);
    y_mean = squeeze(mean(y, 2));
    x = (1:size(y, 3)) / (problem_n + 1);


    switch profile_options.(ProfileOptionKey.RANGE_TYPE.value)
        case 'minmax'
            y_lower = squeeze(min(y, [], 2));
            y_upper = squeeze(max(y, [], 2));
        case 'meanstd'
            y_std = squeeze(std(y, [], 2));
            y_lower = y_mean - profile_options.(ProfileOptionKey.STD_FACTOR.value) * y_std;
            y_upper = y_mean + profile_options.(ProfileOptionKey.STD_FACTOR.value) * y_std;
            if ismaxcv
                y_lower = max(y_lower, 0);
                y_upper = max(y_upper, 0);
            end
        otherwise
            error("Unknown range type.");
    end

    for i_solver = 1:n_solvers
        nextColor = ax.ColorOrder(mod(ax.ColorOrderIndex-1, size(ax.ColorOrder, 1)) + 1, :);

        % When the function values are all zero, use linear scale.
        if ~any(y_mean(i_solver, :))
            plot(ax, x, y_mean(i_solver, :), 'DisplayName', labels{i_solver});
        else
            semilogy(ax, x, y_mean(i_solver, :), 'DisplayName', labels{i_solver});
        end
        hold(ax, 'on');
        if n_runs > 1
            fill(ax, [x, fliplr(x)], [y_lower(i_solver, :), fliplr(y_upper(i_solver, :))], ...
            nextColor, 'FaceAlpha', 0.2, 'EdgeAlpha', 0, 'HandleVisibility', 'off');
        end
    end

    box(ax, 'on');
    legend(ax, 'Location', 'northeast');
    hold(ax, 'off');
    set(ax, 'XLim', [min(x), max(x)]);
    xlabel(ax, 'Number of simplex gradients', 'Interpreter', 'latex');

end