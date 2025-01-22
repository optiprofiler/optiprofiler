function [x_perf, y_perf, ratio_max_perf, x_data, y_data, ratio_max_data, curves] = getExtendedPerformancesDataProfileAxes(work, problem_dimensions, profile_options, curves)
%GETEXTENDEDPERFORMANCESDATAPROFILEAXES computes the axes for the extended performance
%profiles and data profiles. The word "extended" specifically refers to that the data
%at the startpoint and endpoint of the profiles are specially handled.

    [n_problems, n_solvers, n_runs] = size(work);

    denominator_perf = @(i_problem, i_run) min(work(i_problem, :, i_run), [], 'omitnan');
    [x_perf, y_perf, ratio_max_perf] = getPerformanceDataProfileAxes(work, denominator_perf);
    if profile_options.(ProfileOptionKey.SEMILOGX.value)
        % We output the log2(x_perf) and log2(ratio_max_perf). This is because the x-axis is in log2 scale in the performance profile.
        x_perf(isfinite(x_perf)) = log2(x_perf(isfinite(x_perf)));
        if ratio_max_perf > eps && log2(ratio_max_perf) > 1
            ratio_max_perf = log2(ratio_max_perf);
        end
    end
    x_perf(isinf(x_perf)) = 1.1 * ratio_max_perf;
    if profile_options.(ProfileOptionKey.SEMILOGX.value)
        x_perf = [zeros(1, n_solvers); x_perf];
    else
        x_perf = [ones(1, n_solvers); x_perf];
    end
    y_perf = [zeros(1, n_solvers, n_runs); y_perf];
    if n_problems > 0
        x_perf = [x_perf; ones(1, n_solvers) * ratio_max_perf * 1.1];
        y_perf = [y_perf; y_perf(end, :, :)];
    end

    denominator_data = @(i_problem, i_run) problem_dimensions(i_problem) + 1;
    [x_data, y_data, ratio_max_data] = getPerformanceDataProfileAxes(work, denominator_data);
    if profile_options.(ProfileOptionKey.SEMILOGX.value)
        % We output the log2(1 + x_data) and log2(1 + ratio_max_data). This is because the x-axis is in log2 scale in the data profile.
        x_data(isfinite(x_data)) = log2(1 + x_data(isfinite(x_data)));
        if ratio_max_data > eps
            ratio_max_data = log2(1 + ratio_max_data);
        end
    end
    x_data(isinf(x_data)) = 1.1 * ratio_max_data;
    x_data = [zeros(1, n_solvers); x_data];
    y_data = [zeros(1, n_solvers, n_runs); y_data];
    if n_problems > 0
        x_data = [x_data; ones(1, n_solvers) * ratio_max_data * 1.1];
        y_data = [y_data; y_data(end, :, :)];
    end

    % Store the curves in the `curves` struct.
    for i_solver = 1:n_solvers
        for i_run = 1:n_runs
            curves.perf{i_solver, i_run} = [x_perf(:, i_solver)'; y_perf(:, i_solver, i_run)'];
            curves.data{i_solver, i_run} = [x_data(:, i_solver)'; y_data(:, i_solver, i_run)'];
        end
        y_mean_perf = squeeze(mean(y_perf(:, i_solver, :), 3));
        y_mean_data = squeeze(mean(y_data(:, i_solver, :), 3));
        curves.perf{i_solver, n_runs + 1} = [x_perf(:, i_solver)'; y_mean_perf'];
        curves.data{i_solver, n_runs + 1} = [x_data(:, i_solver)'; y_mean_data'];
    end
end