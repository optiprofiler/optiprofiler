function compareOneProblem(solvers, labels, problem_name, custom_problem_loader, feature_name, feature_options)

    % Build the feature.
    feature_options.(FeatureOptionKey.N_RUNS.value) = 1;
    feature = Feature(feature_name, feature_options);

    % Build the profile options.
    profile_options = struct();
    profile_options.(ProfileOptionKey.PROJECT_X0.value) = false;
    profile_options.(ProfileOptionKey.MAX_EVAL_FACTOR.value) = 100;

    if ~isempty(custom_problem_loader)
        problem_name = {'CUSTOM_PROBLEM', problem_name};
    end

    % Solve the problem.
    [fun_histories, maxcv_histories, fun_out, maxcv_out, fun_init, maxcv_init, n_eval, problem_name, problem_n, computation_time] = solveOneProblem(problem_name, solvers, labels, feature, custom_problem_loader, profile_options);

    % Display the results.
    for i_solver = 1:length(solvers)
        flattened_array = reshape(fun_histories(i_solver, 1, 1:n_eval(i_solver, 1)), [], 1);
        flattened_array = cummin(flattened_array);
        semilogy(flattened_array, 'DisplayName', labels{i_solver}, 'LineWidth', 1.5);
        hold on;
    end
    
    ylabel('Function values');
    xlabel('Number of funciton evaluations');
    legend('show');

end