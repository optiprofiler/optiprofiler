function [x_perf, y_perf, ratio_max_perf, x_data, y_data, ratio_max_data] = getExtendedPerformancesDataProfileAxes(work, problem_dimensions)

    [n_problems, n_solvers, n_runs] = size(work);

    denominator_perf = @(i_problem, i_run) min(work(i_problem, :, i_run), [], 'omitnan');
    [x_perf, y_perf, ratio_max_perf] = getPerformanceDataProfileAxes(work, denominator_perf);
    x_perf(isinf(x_perf)) = ratio_max_perf ^ 2;
    x_perf = [ones(1, n_solvers); x_perf];
    y_perf = [zeros(1, n_solvers, n_runs); y_perf];
    if n_problems > 1
        x_perf = [x_perf; ones(1, n_solvers) * (ratio_max_perf ^ 2.0)];
        y_perf = [y_perf; y_perf(end, :, :)];
    end
    % We output the log2(x_perf). This is because the x-axis is in log2 scale.
    x_perf = log2(x_perf);

    denominator_data = @(i_problem, i_run) problem_dimensions(i_problem) + 1;
    [x_data, y_data, ratio_max_data] = getPerformanceDataProfileAxes(work, denominator_data);
    x_data(isinf(x_data)) = 2 * ratio_max_data;
    x_data = [zeros(1, n_solvers); x_data];
    y_data = [zeros(1, n_solvers, n_runs); y_data];
    if n_problems > 1
        x_data = [x_data; ones(1, n_solvers) * (2.0 * ratio_max_data)];
        y_data = [y_data; y_data(end, :, :)];
    end

end