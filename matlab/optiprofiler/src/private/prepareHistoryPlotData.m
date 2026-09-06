function [x_indices, y_values_m, y_values_l, y_values_u, n_runs] = prepareHistoryPlotData(y, is_cum, y_shift, n_eval, profile_options)
%PREPAREHISTORYPLOTDATA Shared numeric presentation for native and portable plots.
%   Y is the bounded display copy from processHistYaxes, not saved raw data.
%   Return per-solver evaluation indices and the same mean/errorband curves
%   for both renderers, including cumulative views and block aggregation.

    n_solvers = size(y, 1);
    n_runs = size(y, 2);
    % Keep (solver, evaluation), including one-solver/one-evaluation cases.
    y_mean = reshape(mean(y, 2), n_solvers, []);

    if strcmp(profile_options.(ProfileOptionKey.ERRORBAR_TYPE.value), 'minmax')
        y_lower = reshape(min(y, [], 2), n_solvers, []);
        y_upper = reshape(max(y, [], 2), n_solvers, []);
    else
        % Preserve MATLAB's sample normalization (N-1 for N>1; zero for N=1).
        % Python's established bands use N instead, so unequal-run bands differ.
        y_std = reshape(std(y, 0, 2), n_solvers, []);
        y_lower = y_mean - y_std;
        y_upper = y_mean + y_std;
    end

    % Compute the band before translation, as computeHistoryYShift does. Computing
    % variance after a large shift can reintroduce negative lower bands.
    y_mean = y_mean + y_shift;
    y_lower = y_lower + y_shift;
    y_upper = y_upper + y_shift;

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

end
