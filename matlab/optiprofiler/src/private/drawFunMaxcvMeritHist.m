function drawFunMaxcvMeritHist(ax, y, solver_names, is_cum, problem_n, y_shift, n_eval, profile_options)
%DRAWFUNMAXCVMERITHIST draws figures of histories of function values, maximum constraint violation, or merit function values.

    line_colors = profile_options.(ProfileOptionKey.LINE_COLORS.value);
    line_styles = profile_options.(ProfileOptionKey.LINE_STYLES.value);
    line_widths = profile_options.(ProfileOptionKey.LINE_WIDTHS.value);

    % Shift the y-axis.
    y = y + y_shift;

    n_solvers = size(y, 1);
    n_runs = size(y, 2);
    y_mean = squeeze(mean(y, 2));

    if strcmp(profile_options.(ProfileOptionKey.ERRORBAR_TYPE.value), 'minmax')
        y_lower = squeeze(min(y, [], 2));
        y_upper = squeeze(max(y, [], 2));
    else
        y_std = squeeze(std(y, 0, 2));
        y_lower = y_mean - y_std;
        y_upper = y_mean + y_std;
    end

    if is_cum
        y_mean = cummin(y_mean, 2);
        y_lower = cummin(y_lower, 2);
        y_upper = cummin(y_upper, 2);
    end

    % =========================================================================
    % Block Aggregation for Large Evaluation Histories
    %
    % When the number of evaluation points is very large, directly plotting all
    % y-values becomes computationally expensive (both time and memory). To
    % avoid this, we apply block aggregation, which means that we group
    % neighboring points into small segments ("blocks") and keep only one
    % representative value from each block. This greatly reduces the total
    % number of points we need to draw, while still keeping the main shape of
    % the curve.
    %
    % There are three ways (modes) to get the value for each block:
    %
    %   1. 'min'
    %      - Find the smallest y_mean value inside the block.
    %      - Its position (x index) is used as the x coordinate (abs_idx).
    %      - The corresponding y_lower and y_upper at that same position are used.
    %      - This highlights the best (lowest) performance in that part of the curve.
    %
    %   2. 'mean'
    %      - Use the middle x position of the block (the median of start and end indices).
    %      - Compute the average of all y_mean, y_lower, and y_upper values in the block.
    %      - This produces a smooth averaged curve that represents the general trend.
    %
    %   3. 'max'
    %      - Find the largest y_mean value inside the block.
    %      - Its position (x index) is used as the x coordinate (abs_idx).
    %      - The corresponding y_lower and y_upper at that same position are used.
    %      - This shows the worst (highest) value in that part of the curve.
    %
    % The first and last points are always kept. Only the points between them
    % are affected by the block aggregation.
    %
    % After doing this, we only keep around `n_blocks` (1000) points instead of
    % all evaluations. This keeps the figure clear and makes plotting faster.
    % =========================================================================
    max_eval = size(y, 3);
    n_blocks = 1000;

    % We initialize the cell array to store the x indices and y values to be plotted.
    x_indices = repmat({[]}, 1, n_solvers);
    y_values_m = repmat({[]}, 1, n_solvers);
    y_values_l = repmat({[]}, 1, n_solvers);
    y_values_u = repmat({[]}, 1, n_solvers);

    % We will build blocks only for 2:(max_eval-1) so the first and last points are excluded
    % from block aggregation and will be appended unconditionally at the end.
    if max_eval > 2
        inner_len = max_eval - 2;
        n_blocks_eff = min(n_blocks, inner_len);
        q = floor(inner_len / n_blocks_eff);
        r = mod(inner_len, n_blocks_eff);
        blocks = q * ones(1, n_blocks_eff);
        if r > 0
            % Distribute the remaining points (r) evenly across the blocks.
            % We choose roughly r evenly spaced block indices using linspace,
            % then add one extra element to each of these blocks so that
            % the total number of elements sums back to inner_len.
            idxs = round(linspace(1, n_blocks_eff, r));
            blocks(idxs) = blocks(idxs) + 1;
        end

        % We aggregate over inner blocks (indices 2, ..., max_eval-1).
        for i_block = 1:n_blocks_eff
            idx_start = sum(blocks(1:i_block-1)) + 1 + 1; % Add additional 1 to skip the first point.
            idx_end = idx_start + blocks(i_block) - 1;
            idx = idx_start:idx_end;
            for i_solver = 1:n_solvers
                i_eval = max(n_eval(i_solver,:));
                switch profile_options.(ProfileOptionKey.HIST_AGGREGATION.value)
                    case 'min'
                        [y_value_m, rel_idx] = min(y_mean(i_solver, idx), [], 'omitnan');
                        abs_idx = idx(rel_idx);
                        y_value_l = y_lower(i_solver, abs_idx);
                        y_value_u = y_upper(i_solver, abs_idx);
                    case 'mean'
                        abs_idx = floor((idx_start + idx_end) / 2);
                        y_value_m = mean(y_mean(i_solver, idx), 'omitnan');
                        y_value_l = mean(y_lower(i_solver, idx), 'omitnan');
                        y_value_u = mean(y_upper(i_solver, idx), 'omitnan');
                    case 'max'
                        [y_value_m, rel_idx] = max(y_mean(i_solver, idx), [], 'omitnan');
                        abs_idx = idx(rel_idx);
                        y_value_l = y_lower(i_solver, abs_idx);
                        y_value_u = y_upper(i_solver, abs_idx);
                end
                if abs_idx > i_eval
                    continue;
                end
                x_indices{i_solver} = [x_indices{i_solver}, abs_idx];
                y_values_m{i_solver} = [y_values_m{i_solver}, y_value_m];
                y_values_l{i_solver} = [y_values_l{i_solver}, y_value_l];
                y_values_u{i_solver} = [y_values_u{i_solver}, y_value_u];
            end
        end
    end

    % We add the first and last indices unconditionally.
    for i_solver = 1:n_solvers
        i_eval = max(n_eval(i_solver,:));
        if i_eval > 0
            % Add the first point.
            x_indices{i_solver} = [1, x_indices{i_solver}];
            y_values_m{i_solver} = [y_mean(i_solver, 1), y_values_m{i_solver}];
            y_values_l{i_solver} = [y_lower(i_solver, 1), y_values_l{i_solver}];
            y_values_u{i_solver} = [y_upper(i_solver, 1), y_values_u{i_solver}];

            % Add the last point.
            if i_eval ~= 1
                x_indices{i_solver} = [x_indices{i_solver}, i_eval];
                y_values_m{i_solver} = [y_values_m{i_solver}, y_mean(i_solver, i_eval)];
                y_values_l{i_solver} = [y_values_l{i_solver}, y_lower(i_solver, i_eval)];
                y_values_u{i_solver} = [y_values_u{i_solver}, y_upper(i_solver, i_eval)];
            end
        end
    end

    xl_lim = 1 / (problem_n + 1);
    xr_lim = 1 / (problem_n + 1);
    is_log_scale = false;

    hold(ax, 'on');
    for i_solver = 1:n_solvers
        % Truncate the histories according to the function evaluations of each solver.
        i_x = x_indices{i_solver};
        i_y_mean = y_values_m{i_solver};
        i_y_lower = y_values_l{i_solver};
        i_y_upper = y_values_u{i_solver};
        i_eval = length(i_x);
        x = i_x / (problem_n + 1);
        xr_lim = max(xr_lim, x(end));

        color = line_colors(mod(i_solver - 1, size(line_colors, 1)) + 1, :);
        line_style = line_styles{mod(i_solver - 1, size(line_styles, 2)) + 1};
        line_width = line_widths(mod(i_solver - 1, length(line_widths)) + 1);

        if n_runs > 1 && i_eval > 1
            fill(ax, [x, fliplr(x)], [i_y_lower, fliplr(i_y_upper)], color, 'FaceAlpha', 0.2, 'EdgeAlpha', 0, 'HandleVisibility', 'off');
        end
        if i_eval == 1
            plot(ax, x, i_y_mean, 'o', 'Color', color, 'DisplayName', solver_names{i_solver});
        elseif i_eval > 1
            plot(ax, x, i_y_mean, line_style, 'Color', color, 'LineWidth', line_width, 'DisplayName', solver_names{i_solver});
        end
        if any(i_y_mean) && any(diff(i_y_mean))
            is_log_scale = true;
        end
    end

    % When the function values are not all zero and there is at least some change in the function values, use log scale for the y-axis.
    if is_log_scale
        set(ax, 'YScale', 'log');
    end

    box(ax, 'on');
    legend(ax, 'Location', 'northeast');
    hold(ax, 'off');
    if abs(xl_lim - xr_lim) < eps
        xr_lim = xl_lim + eps;
    end
    set(ax, 'XLim', [xl_lim, xr_lim]);
    xlabel(ax, profile_options.(ProfileOptionKey.XLABEL_DATA_PROFILE.value), 'Interpreter', 'latex');

end