function [x_perf, y_perf, ratio_max_perf, x_data, y_data, ratio_max_data, profiles] = getExtendedPerformancesDataProfileAxes(work, problem_dimensions, profiles)
%GETEXTENDEDPERFORMANCESDATAPROFILEAXES computes the axes for the extended performance
%profiles and data profiles. The word "extended" specifically refers to that the data
%at the startpoint and endpoint of the profiles are specially handled.

    [n_problems, n_solvers, n_runs] = size(work);

    denominator_perf = @(i_problem, i_run) min(work(i_problem, :, i_run), [], 'omitnan');
    [x_perf, y_perf, ratio_max_perf] = getPerformanceDataProfileAxes(work, denominator_perf, 'perf');
    % We output the log2(x_perf) and log2(ratio_max_perf). This is because the x-axis is in log2 scale in the performance profile.
    x_perf(isfinite(x_perf)) = log2(x_perf(isfinite(x_perf)));
    x_perf(isinf(x_perf)) = 1.1 * ratio_max_perf;
    x_perf = [zeros(1, n_solvers); x_perf];
    y_perf = [zeros(1, n_solvers, n_runs); y_perf];
    if n_problems > 0
        x_perf = [x_perf; ones(1, n_solvers) * ratio_max_perf * 1.1];
        y_perf = [y_perf; y_perf(end, :, :)];
    end

    denominator_data = @(i_problem, i_run) problem_dimensions(i_problem) + 1;
    [x_data, y_data, ratio_max_data] = getPerformanceDataProfileAxes(work, denominator_data, 'data');
    x_data(isinf(x_data)) = 1.1 * ratio_max_data;
    x_data = [zeros(1, n_solvers); x_data];
    y_data = [zeros(1, n_solvers, n_runs); y_data];
    if n_problems > 0
        x_data = [x_data; ones(1, n_solvers) * ((1 + ratio_max_data) ^ 1.1 - 1)];
        y_data = [y_data; y_data(end, :, :)];
    end
    % We output the log2(1 + x_data) and log2(1 + ratio_max_data). This is because the x-axis is in log2 scale in the data profile.
    x_data = log2(1 + x_data);
    ratio_max_data = log2(1 + ratio_max_data);

    % Store the curves in the `profiles` struct.
    for i_solver = 1:n_solvers
        for i_run = 1:n_runs
            profiles.perf{i_solver, i_run} = [x_perf(:, i_solver)'; y_perf(:, i_solver, i_run)'];
            profiles.data{i_solver, i_run} = [x_data(:, i_solver)'; y_data(:, i_solver, i_run)'];
        end
        y_mean_perf = squeeze(mean(y_perf(:, i_solver, :), 3));
        y_mean_data = squeeze(mean(y_data(:, i_solver, :), 3));
        profiles.perf{i_solver, n_runs + 1} = [x_perf(:, i_solver)'; y_mean_perf'];
        profiles.data{i_solver, n_runs + 1} = [x_data(:, i_solver)'; y_mean_data'];
    end
end