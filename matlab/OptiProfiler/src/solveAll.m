function [fun_values, maxcv_values, fun_inits, maxcv_inits, n_evals, problem_names, problem_dimensions] = solveAll(problem_names, problem_options, solvers, labels, feature, max_eval_factor, profile_options)
    %SOLVEALL is defined to solve all problems in the given problem set.

    % Initialize output variables
    fun_values = [];
    maxcv_values = [];
    fun_inits = [];
    maxcv_inits = [];
    n_evals = [];

    % Solve all problems.
    n_problems = length(problem_names);
    results = cell(1, n_problems);

    switch profile_options.(ProfileOptionKey.N_JOBS.value)
        case 1
            % Do not use parallel computing.
            for i_problem = 1:n_problems
                problem_name = problem_names{i_problem};
                [tmp_fun_values, tmp_maxcv_values, tmp_fun_init, tmp_maxcv_init, tmp_n_eval, tmp_problem_name, tmp_problem_n] = solveOne(problem_name, problem_options, solvers, labels, feature, max_eval_factor);
                results{i_problem} = {tmp_fun_values, tmp_maxcv_values, tmp_fun_init, tmp_maxcv_init, tmp_n_eval, tmp_problem_name, tmp_problem_n};
            end
        otherwise
            parpool(profile_options.(ProfileOptionKey.N_JOBS.value));
            parfor i_problem = 1:n_problems
                problem_name = problem_names{i_problem};
                [tmp_fun_values, tmp_maxcv_values, tmp_fun_init, tmp_maxcv_init, tmp_n_eval, tmp_problem_name, tmp_problem_n] = solveOne(problem_name, problem_options, solvers, labels, feature, max_eval_factor);
                results{i_problem} = {tmp_fun_values, tmp_maxcv_values, tmp_fun_init, tmp_maxcv_init, tmp_n_eval, tmp_problem_name, tmp_problem_n};
            end
            delete(gcp);
    end

    
    % TODO: NaN(n_problems,...) may have problems, since some results{i} may be empty and we should delete them.
    % Process results.
    n_runs = feature.options.(FeatureOptionKey.N_RUNS.value);
    problem_dimensions = NaN(n_problems, 1);
    for i_problem = 1:n_problems
        problem_dimensions(i_problem) = results{i_problem}{7};
    end
    max_eval = max_eval_factor * max(problem_dimensions);
    fun_values = NaN(n_problems, length(solvers), n_runs, max_eval);
    maxcv_values = NaN(n_problems, length(solvers), n_runs, max_eval);
    fun_inits = NaN(n_problems, 1);
    maxcv_inits = NaN(n_problems, 1);
    n_evals = NaN(n_problems, length(solvers), n_runs);

    for i_problem = 1:n_problems
        result = results{i_problem};
        fun_value = result{1};
        maxcv_value = result{2};
        n_eval = size(fun_value, 3);
        fun_values(i_problem, :, :, 1:n_eval) = fun_value;
        maxcv_values(i_problem, :, :, 1:n_eval) = maxcv_value;
        fun_inits(i_problem) = result{3};
        maxcv_inits(i_problem) = result{4};
        n_evals(i_problem, :, :) = result{5};
        if n_eval > 0 && n_eval < max_eval
            slice = fun_values(i_problem, :, :, n_eval);
            slice = repmat(slice, [1, 1, 1, max_eval-n_eval]);
            fun_values(i_problem, :, :, n_eval+1:end) = slice;
            slice = maxcv_values(i_problem, :, :, n_eval);
            slice = repmat(slice, [1, 1, 1, max_eval-n_eval]);
            maxcv_values(i_problem, :, :, n_eval+1:end) = slice;
        end
    end
    
end