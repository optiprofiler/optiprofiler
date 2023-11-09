function createProfiles(solvers, labels, problem_names, feature_name)
    %CREATEPROFILES plots performance profiles and data profiles.

    % TODO: the current script version does not support options, only supports "simple" features.

    % Set number of pools to be the same as the number of cores.

    profile_options = struct(ProfileOptionKey.N_JOBS.value, 0);
    % profile_options = struct(ProfileOptionKey.N_JOBS.value, feature('numcores'));

    % Build feature.
    % TODO: now we first use "simple" default feature.
    feature = Feature(feature_name);

    % Solve all the problems.
    max_eval_factor = 500;
    [fun_values, maxcv_values, problem_names, problem_dimensions] = solveAll(problem_names, solvers, feature, max_eval_factor, profile_options);

    % Compute the merit values.
    merit_values = computeMeritValues(fun_values, maxcv_values);

    % Extract the merit values at the initial points.
    merit_init = min(merit_values(:, :, :, 1), [], 2, 'omitnan');

    % Determine the least merit value for each problem.
    merit_min = min(min(merit_values, [], 4, 'omitnan'), [], 2, 'omitnan');
    if ismember(feature_name, {FeatureName.NOISY.value, FeatureName.TOUGH.value, FeatureName.TRUNCATED.value})
        feature_plain = Feature(FeatureName.PLAIN.value);
        [fun_values_plain, maxcv_values_plain, ~, ~] = solveAll(problem_names, solvers, feature_plain, max_eval_factor, profile_options);
        merit_values_plain = computeMeritValues(fun_values_plain, maxcv_values_plain);
        merit_min_plain = min(min(merit_values_plain, [], 4, 'omitnan'), [], 2, 'omitnan');
        merit_min = min(merit_min, merit_min_plain, 'omitnan');
    end

    [n_problems, n_solvers, n_runs, max_eval] = size(merit_values);
    tolerances = 10.^([-1:-1:-10]);
    for i = 1:10
        i_profile = i;
        tolerance = tolerances(i);

        work = NaN(n_problems, n_solvers, n_runs);
        for i_problem = 1:n_problems
            for i_solver = 1:n_solvers
                for i_run = 1:n_runs
                    if isfinite(merit_min(i_problem))
                        threshold = max(tolerance * merit_init(i_problem, i_run) * merit_min(i_problem, i_run), merit_min(i_problem, i_run));
                    else
                        threshold = -Inf;
                    end
                    if min(merit_values(i_problem, i_solver, i_run, :), [], 'omitnan') <= threshold
                        work(i_problem, i_solver, i_run) = find(merit_values(i_problem, i_solver, i_run, :) <= threshold, 1, 'last') + 1;
                    end
                end
            end
        end

        % Calculate the x-axes of the performance profiles.
        x_perf = NaN(n_runs, n_problems, n_solvers);
        for i_run = 1:n_runs
            for i_problem = 1:n_problems
                if ~all(isnan(work(i_problem, :, i_run)))
                    x_perf(i_run, i_problem, :) = work(i_problem, :, i_run) / min(work(i_problem, :, i_run), [], 'omitnan');
                end
            end
        end
        perf_ratio_max = max(max(x_perf, [], 'all', 'omitnan'), 2^eps('double'));
        x_perf(isnan(x_perf)) = 2 * perf_ratio_max;
        x_perf = sort(x_perf, 2);
        x_perf = reshape(x_perf, n_problems * n_runs, n_solvers);
        [~, sort_perf] = sort(x_perf, 1);
        for i = 1:size(x_perf,2)
            x_perf(:,i) = x_perf(sort_perf(:,i),i);
        end

        % Calculate the y-axes of the performance profiles.
        y_perf = zeros(n_problems * n_runs, n_solvers);
        for i_run = 1:n_runs
            for i_solver = 1:n_solvers
                y = NaN(n_problems * n_runs, 1);
                y(i_run * n_problems + 1 : (i_run + 1) * n_problems) = linspace(1 / n_problems, 1, n_problems);
                y = y(sort_perf(:, i_solver));
                for i_problem = 1:n_problems*n_runs
                    if isnan(y(i_problem))
                        if i_problem == 1
                            y(i_problem) = 0;
                        else
                            y(i_problem) = y(i_problem - 1);
                        end
                    end
                end
                y_perf(:, i_solver) = y_perf(:, i_solver) + y;
            end
        end
        y_perf = y_perf / n_runs;

        % Calculate the x-axes of the data profiles.
        x_data = NaN(n_runs, n_problems, n_solvers);
        for i_run = 1:n_runs
            for i_problem = 1:n_problems
                if ~all(isnan(work(i_problem, :, i_run)))
                    x_data(i_run, i_problem, :) = work(i_problem, :, i_run) / (problem_dimensions(i_problem) + 1);
                end
            end
        end
        data_ratio_max = max(max(x_data, [], 'all', 'omitnan'), eps('double'));
        x_data(isnan(x_data)) = 2 * data_ratio_max;
        x_data = sort(x_data, 2);
        x_data = reshape(x_data, n_problems * n_runs, n_solvers);
        [~, sort_data] = sort(x_data, 1);
        for i = 1:size(x_data,2)
            x_data(:,i) = x_data(sort_data(:,i),i);
        end

        % Calculate the y-axes of the data profiles.
        y_data = zeros(n_problems * n_runs, n_solvers);
        for i_runs = 1:n_runs
            for i_solver = 1:n_solvers
                y = NaN(n_problems * n_runs, 1);
                y(i_run * n_problems + 1 : (i_run + 1) * n_problems) = linspace(1 / n_problems, 1, n_problems);
                y = y(sort_data(:, i_solver));
                for i_problem = 1:n_problems*n_runs
                    if isnan(y(i_problem))
                        if i_problem == 1
                            y(i_problem) = 0;
                        else
                            y(i_problem) = y(i_problem - 1);
                        end
                    end
                end
                y_data(:, i_solver) = y_data(:, i_solver) + y;
            end
        end
        y_data = y_data / n_runs;

        % Plot the performance profiles.
        figure
        hold on
        for j = 1:n_solvers
            x = repelem(x_perf(:, j), 2);
            x = [0.0, x(1), x(2:end)', 2.0 * perf_ratio_max];
            y = repelem(y_perf(:, j), 2);
            y = [0.0, 0.0, y(1:end-1)', y(end)];
            semilogx(x, y, 'DisplayName', labels{j})
        end
        xlim([1.0, 1.1 * perf_ratio_max])
        ylim([0.0, 1.0])
        xlabel('Performance ratio')
        ylabel('Performance profiles')
        legend('Location','southwest')
        hold off

        % Plot the data profiles.
        figure
        hold on
        for j = 1:n_solvers
            x = repelem(x_data(:, j), 2);
            x = [0.0, x(1), x(2:end)', 2.0 * data_ratio_max];
            y = repelem(y_data(:, j), 2);
            y = [0.0, 0.0, y(1:end-1)', y(end)];
            plot(x, y, 'DisplayName', labels{j})
        end
        xlim([0.0, 1.1 * data_ratio_max])
        ylim([0.0, 1.0])
        xlabel('Number of simplex gradients')
        ylabel('Data profiles')
        legend('Location','southwest')
        hold off

    end


end


function merit_values = computeMeritValues(fun_values, maxcv_values)
    is_nearly_feasible = maxcv_values <= 1e-12;
    is_very_infeasible = maxcv_values >= 1e-6;
    is_undecided = ~is_nearly_feasible & ~is_very_infeasible;
    merit_values = NaN(size(fun_values));
    merit_values(is_nearly_feasible) = fun_values(is_nearly_feasible);
    merit_values(is_very_infeasible) = Inf;
    merit_values(is_undecided) = fun_values(is_undecided) + 1e8 * maxcv_values(is_undecided);
end