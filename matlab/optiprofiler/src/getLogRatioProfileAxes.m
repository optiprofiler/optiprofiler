function [x_log_ratio, y_log_ratio, ratio_max, profiles] = getLogRatioProfileAxes(work, profiles)
%GETLOGRATIOPROFILEAXES computes the axes for the log-ratio profiles.

    [n_problems, n_solvers, n_runs] = size(work);
    work_flat = reshape(permute(work, [1, 3, 2]), n_problems * n_runs, n_solvers);
    y_log_ratio = NaN(n_problems * n_runs, 1);
    log_ratio_finite = isfinite(work_flat(:, 1)) & isfinite(work_flat(:, 2));
    y_log_ratio(log_ratio_finite) = log2(work_flat(log_ratio_finite, 1)) - log2(work_flat(log_ratio_finite, 2));
    ratio_max = max(max(abs(y_log_ratio(log_ratio_finite)), [], 'all'), eps);
    if isempty(ratio_max)
        ratio_max = eps;
    end
    y_log_ratio(isnan(work_flat(:, 1)) & isfinite(work_flat(:, 2))) = 1.1 * ratio_max;
    y_log_ratio(isfinite(work_flat(:, 1)) & isnan(work_flat(:, 2))) = -1.1 * ratio_max;
    y_log_ratio(isnan(work_flat(:, 1)) & isnan(work_flat(:, 2))) = 0.0;
    y_log_ratio = sort(y_log_ratio);

    x_log_ratio = (1:(n_problems * n_runs))';

    % Store the curves in the `profiles` struct.
    profiles.log_ratio{1} = [x_log_ratio(y_log_ratio < 0)'; y_log_ratio(y_log_ratio < 0)'];
    profiles.log_ratio{2} = [x_log_ratio(y_log_ratio > 0)'; y_log_ratio(y_log_ratio > 0)'];
end