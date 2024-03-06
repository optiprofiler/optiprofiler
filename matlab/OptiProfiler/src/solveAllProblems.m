function [fun_histories, maxcv_histories, fun_out, maxcv_out, fun_init, maxcv_init, n_eval, problem_names, problem_dimensions] = solveAllProblems(problem_names, problem_options, solvers, labels, feature, max_eval_factor, profile_options)
%SOLVEALLPROBLEMS solves all problems in the problem_names list using solvers in the solvers list and stores the computing results.

    % Solve all problems.
    n_problems = length(problem_names);
    results = cell(1, n_problems);

    switch profile_options.(ProfileOptionKey.N_JOBS.value)
        case 1
            % Do not use parallel computing.
            for i_problem = 1:n_problems
                problem_name = problem_names{i_problem};
                [tmp_fun_histories, tmp_maxcv_histories, tmp_fun_out, tmp_maxcv_out, tmp_fun_init, tmp_maxcv_init, tmp_n_eval, tmp_problem_name, tmp_problem_n] = solveOneProblem(problem_name, solvers, labels, feature, max_eval_factor, problem_options);
                results{i_problem} = {tmp_fun_histories, tmp_maxcv_histories, tmp_fun_out, tmp_maxcv_out, tmp_fun_init, tmp_maxcv_init, tmp_n_eval, tmp_problem_name, tmp_problem_n};
            end
        otherwise
            parpool(profile_options.(ProfileOptionKey.N_JOBS.value));
            parfor i_problem = 1:n_problems
                problem_name = problem_names{i_problem};
                [tmp_fun_histories, tmp_maxcv_histories, tmp_fun_out, tmp_maxcv_out, tmp_fun_init, tmp_maxcv_init, tmp_n_eval, tmp_problem_name, tmp_problem_n] = solveOneProblem(problem_name, solvers, labels, feature, max_eval_factor, problem_options);
                results{i_problem} = {tmp_fun_histories, tmp_maxcv_histories, tmp_fun_out, tmp_maxcv_out, tmp_fun_init, tmp_maxcv_init, tmp_n_eval, tmp_problem_name, tmp_problem_n};
            end
            delete(gcp);
            fprintf("Leaving the parallel section.\n");
    end

    if all(cellfun(@isempty, results))
        fprintf("All problems failed to load.\n");
        fun_histories = [];
        maxcv_histories = [];
        fun_out = [];
        maxcv_out = [];
        fun_init = [];
        maxcv_init = [];
        n_eval = [];
        problem_names = [];
        problem_dimensions = [];
    else
        % Process results.
        n_solvers = length(solvers);
        n_runs = feature.options.(FeatureOptionKey.N_RUNS.value);
        problem_dimensions = NaN(n_problems, 1);
        for i_problem = 1:n_problems
            problem_dimensions(i_problem) = results{i_problem}{9};
        end
        if length(problem_dimensions) > 0
            max_eval = max_eval_factor * max(problem_dimensions);
        else
            max_eval = 1;
        end
        fun_histories = NaN(n_problems, n_solvers, n_runs, max_eval);
        maxcv_histories = NaN(n_problems, n_solvers, n_runs, max_eval);
        fun_out = NaN(n_problems, n_solvers, n_runs);
        maxcv_out = NaN(n_problems, n_solvers, n_runs);
        fun_init = NaN(n_problems, 1);
        maxcv_init = NaN(n_problems, 1);
        n_eval = NaN(n_problems, n_solvers, n_runs);

        for i_problem = 1:n_problems
            max_eval = max_eval_factor * problem_dimensions(i_problem);
            result = results{i_problem};
            fun_hist = result{1};
            maxcv_hist = result{2};
            fun_histories(i_problem, :, :, 1:max_eval) = fun_hist;
            maxcv_histories(i_problem, :, :, 1:max_eval) = maxcv_hist;
            fun_out(i_problem, :, :) = result{3};
            maxcv_out(i_problem, :, :) = result{4};
            fun_init(i_problem) = result{5};
            maxcv_init(i_problem) = result{6};
            n_eval(i_problem, :, :) = result{7};
            if max_eval > 0
                fun_histories(i_problem, :, :, max_eval+1:end) = repmat(fun_histories(i_problem, :, :, max_eval), 1, 1, 1, size(fun_histories, 4) - max_eval);
                maxcv_histories(i_problem, :, :, max_eval+1:end) = repmat(maxcv_histories(i_problem, :, :, max_eval), 1, 1, 1, size(maxcv_histories, 4) - max_eval);
            end
        end
    end
    
end