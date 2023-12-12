function profileinprima(solvers, labels, problem_names, feature_name)
    %CREATEPROFILES plots performance profiles and data profiles.

    % TODO: the current script version does not support options, only supports "simple" features.

    % Set number of pools to be the same as the number of cores.

    profile_options = struct(ProfileOptionKey.N_JOBS.value, 1);
    % profile_options = struct(ProfileOptionKey.N_JOBS.value, feature('numcores'));

    % Build feature.
    % TODO: now we first use "simple" default feature.
    feature = Feature(feature_name);

    % Solve all the problems.
    max_eval_factor = 500;
    [fun_values, maxcv_values, ~, ~, ~, problem_names, problem_dimensions] = solveAll(problem_names, solvers, feature, max_eval_factor, profile_options);

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

    options.tau = 1e-1;
    options.natural_stop = false;
    options.solvers = labels;
    % output = perfprof(merit_values, merit_min, options)



    output = perfdata(solvers)

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