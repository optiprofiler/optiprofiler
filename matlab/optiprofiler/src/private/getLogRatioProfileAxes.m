function [x_log_ratio, y_log_ratio, ratio_max, n_solvers_fail, curves, bar_sources] = getLogRatioProfileAxes(work, curves)
%GETLOGRATIOPROFILEAXES computes the axes for the log-ratio profiles.

    [n_problems, n_solvers, n_runs] = size(work);
    work_flat = reshape(permute(work, [1, 3, 2]), n_problems * n_runs, n_solvers);
    y_log_ratio = NaN(n_problems * n_runs, 1);
    log_ratio_finite = isfinite(work_flat(:, 1)) & isfinite(work_flat(:, 2)) & work_flat(:, 1) > 0 & work_flat(:, 2) > 0;
    y_log_ratio(log_ratio_finite) = log2(work_flat(log_ratio_finite, 1)) - log2(work_flat(log_ratio_finite, 2));
    ratio_max = max(max(abs(y_log_ratio(log_ratio_finite)), [], 'all'), eps);
    if isempty(ratio_max) || ~isfinite(ratio_max) || ratio_max <= 0
        ratio_max = eps;
    end
    fail1 = ~isfinite(work_flat(:, 1)) | work_flat(:, 1) <= 0;
    fail2 = ~isfinite(work_flat(:, 2)) | work_flat(:, 2) <= 0;
    y_log_ratio(fail1 & ~fail2) = 1.1 * ratio_max;
    y_log_ratio(~fail1 & fail2) = -1.1 * ratio_max;

    % If both solvers fail to solve one problem (in this case log-ratio is NaN), we let
    % log-ratio of both solvers to be 1.1 * ratio_max or -1.1 * ratio_max.
    n_solvers_fail = sum(fail1 & fail2);
    y_log_ratio(fail1 & fail2) = 1.1 * ratio_max;
    y_log_ratio = [y_log_ratio; -1.1 * ratio_max * ones(n_solvers_fail, 1)];

    if nargout > 5
        [y_log_ratio, order] = sort(y_log_ratio);
    else
        y_log_ratio = sort(y_log_ratio);
    end
    x_log_ratio = (1:length(y_log_ratio))';
    if nargout > 5
        % Retain identity at the actual sort, never guess a problem from an
        % anonymous bar. Both-failed rows intentionally create two bars.
        source = [(1:n_problems*n_runs)'; find(fail1 & fail2)];
        source = source(order);
        bar_sources = struct('problem_index', mod(source-1,n_problems)+1, ...
            'run_index', floor((source-1)/n_problems)+1, ...
            'solver1_failed', fail1(source), 'solver2_failed', fail2(source), ...
            'tie', ~fail1(source) & ~fail2(source) & work_flat(source,1) == work_flat(source,2));
    end

    % Store the curves in the `profiles` struct.
    curves.log_ratio{1} = [x_log_ratio(y_log_ratio < 0)'; y_log_ratio(y_log_ratio < 0)'];
    curves.log_ratio{2} = [x_log_ratio(y_log_ratio > 0)'; y_log_ratio(y_log_ratio > 0)'];
end
