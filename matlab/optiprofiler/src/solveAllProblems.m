function [fun_histories, maxcv_histories, fun_out, maxcv_out, fun_init, maxcv_init, n_eval, problem_names, problem_dimensions, computation_times, problem_unsolved] = solveAllProblems(cutest_problem_names, custom_problem_loader, custom_problem_names, solvers, labels, feature, profile_options, is_plot, path_hist_plots)
%SOLVEALLPROBLEMS solves all problems in the problem_names list using solvers in the solvers list and stores the computing results.


    extra_problem_names = arrayfun(@(i) {sprintf('EXTRA%d', i), custom_problem_names{i}}, 1:length(custom_problem_names), 'UniformOutput', false);
    extra_problem_names = reshape(extra_problem_names, 1, []);
    cutest_problem_names = reshape(cutest_problem_names, 1, []);
    problem_names = [cutest_problem_names, extra_problem_names];

    % Solve all problems.
    n_problems = length(problem_names);
    len_problem_names = max(cellfun(@length, problem_names));
    max_eval_factor = profile_options.(ProfileOptionKey.MAX_EVAL_FACTOR.value);
    results = cell(1, n_problems);
    if ~profile_options.(ProfileOptionKey.SILENT.value)
        if n_problems > 1
            fprintf("INFO: There are %d problems to solve.\n\n", n_problems);
        else
            fprintf("INFO: There is only %d problem to solve.\n\n", n_problems);
        end
    end

    % Decide whether to delete the pool.
    try
        pool = gcp('nocreate');
        % We will delete the pool only when keep_pool is false and the number of workers is not equal to n_jobs.
        if ~isempty(pool) && (profile_options.n_jobs == 1 || profile_options.n_jobs ~= pool.NumWorkers) && ~profile_options.(ProfileOptionKey.KEEP_POOL.value)
            if ~profile_options.(ProfileOptionKey.SILENT.value)
                delete(pool);
            else
                evalc("delete(pool)");
            end
        end
    catch
        pool = [];
    end

    switch profile_options.(ProfileOptionKey.N_JOBS.value)
        case 1
            % Do not use parallel computing.
            if ~isempty(pool)
                if ~profile_options.(ProfileOptionKey.SILENT.value) && ~profile_options.(ProfileOptionKey.KEEP_POOL.value)
                    delete(pool);
                else
                    evalc("delete(pool)");
                end
            end
            for i_problem = 1:n_problems
                problem_name = problem_names{i_problem};
                [tmp_fun_histories, tmp_maxcv_histories, tmp_fun_out, tmp_maxcv_out, tmp_fun_init, tmp_maxcv_init, tmp_n_eval, tmp_problem_name, tmp_problem_n, tmp_computation_time] = solveOneProblem(problem_name, solvers, labels, feature, len_problem_names, custom_problem_loader, profile_options, is_plot, path_hist_plots);
                results{i_problem} = {tmp_fun_histories, tmp_maxcv_histories, tmp_fun_out, tmp_maxcv_out, tmp_fun_init, tmp_maxcv_init, tmp_n_eval, tmp_problem_name, tmp_problem_n, tmp_computation_time};
            end
        otherwise
            try
                pool = gcp('nocreate');
            catch
                % If the user does not have the Parallel Computing Toolbox, the gcp function will not be available. We will print a error message and exit.
                error("The Parallel Computing Toolbox is not available. Please set the n_jobs option to 1.");
            end
            if isempty(pool)
                if ~profile_options.(ProfileOptionKey.SILENT.value)
                    parpool(profile_options.(ProfileOptionKey.N_JOBS.value));
                else
                    evalc("parpool(profile_options.(ProfileOptionKey.N_JOBS.value))");
                end
            end
            parfor i_problem = 1:n_problems
                problem_name = problem_names{i_problem};
                [tmp_fun_histories, tmp_maxcv_histories, tmp_fun_out, tmp_maxcv_out, tmp_fun_init, tmp_maxcv_init, tmp_n_eval, tmp_problem_name, tmp_problem_n, tmp_computation_time] = solveOneProblem(problem_name, solvers, labels, feature, len_problem_names, custom_problem_loader, profile_options, is_plot, path_hist_plots);
                results{i_problem} = {tmp_fun_histories, tmp_maxcv_histories, tmp_fun_out, tmp_maxcv_out, tmp_fun_init, tmp_maxcv_init, tmp_n_eval, tmp_problem_name, tmp_problem_n, tmp_computation_time};
            end
            if ~profile_options.(ProfileOptionKey.KEEP_POOL.value)
                if ~profile_options.(ProfileOptionKey.SILENT.value)
                    delete(gcp);
                    fprintf("Leaving the parallel section.\n");
                else
                    evalc("delete(gcp)");
                end
            end
    end

    % Select and delete the empty results. Get a list of unsolved problems.
    not_empty_index = cellfun(@(c) all(cellfun(@(x) ~isempty(x), c)), results);
    problem_unsolved = problem_names(~not_empty_index);
    results = results(not_empty_index);

    if all(arrayfun(@(x) all(cellfun(@isempty, results{x}([1:7, 9:10]))), 1:length(results))) % Check if all problems failed to load.
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
        computation_times = [];
    else
        % Process results.
        n_problems = length(results);
        n_solvers = length(solvers);
        n_runs = feature.options.(FeatureOptionKey.N_RUNS.value);
        problem_dimensions = NaN(n_problems, 1);
        computation_times = NaN(n_problems, 1);
        for i_problem = 1:n_problems
            problem_dimensions(i_problem) = results{i_problem}{9};
            computation_times(i_problem) = results{i_problem}{10};
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
        problem_names = [];

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
            problem_names{i_problem} = result{8};
            if max_eval > 0
                fun_histories(i_problem, :, :, max_eval+1:end) = repmat(fun_histories(i_problem, :, :, max_eval), 1, 1, 1, size(fun_histories, 4) - max_eval);
                maxcv_histories(i_problem, :, :, max_eval+1:end) = repmat(maxcv_histories(i_problem, :, :, max_eval), 1, 1, 1, size(maxcv_histories, 4) - max_eval);
            end
        end
    end
    
end