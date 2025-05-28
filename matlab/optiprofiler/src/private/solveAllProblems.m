function results = solveAllProblems(solvers, plib, feature, problem_options, profile_options, is_plot, path_hist_plots)
%SOLVEALLPROBLEMS solves all problems in plib satisfying problem_options using solvers in the solvers and stores the computing results.

    results = struct();

    % Get satisfied problem names.
    option_select = problem_options;
    if isfield(option_select, ProblemOptionKey.PLIBS.value)
        option_select = rmfield(option_select, ProblemOptionKey.PLIBS.value);
    end
    if isfield(option_select, ProblemOptionKey.PROBLEM_NAMES.value)
        option_select = rmfield(option_select, ProblemOptionKey.PROBLEM_NAMES.value);
    end
    selector_name = [plib, '_select'];
    select = str2func(selector_name);
    try
        problem_names = {};
        if isfield(problem_options, ProblemOptionKey.PROBLEM_NAMES.value)
            problem_names = problem_options.(ProblemOptionKey.PROBLEM_NAMES.value);
        end
        % Try to use selector function to select problems.
        selected_problem_names = select(option_select);
        % Make sure selected_problem_names is a cell row vector.
        if size(selected_problem_names, 1) > 1
            selected_problem_names = selected_problem_names';
        end
        problem_names = unique([problem_names, selected_problem_names]);
    catch
    end

    if isempty(problem_names)
        if ~profile_options.(ProfileOptionKey.SILENT.value)
            fprintf('\nINFO: No problem is selected from "%s".\n', plib);
        end
        return;
    end

    % Start solving problems.
    loader_name = [plib, '_load'];
    load = str2func(loader_name);
    n_problems = length(problem_names);
    len_problem_names = max(cellfun(@length, problem_names));
    max_eval_factor = profile_options.(ProfileOptionKey.MAX_EVAL_FACTOR.value);
    tmp_results = cell(1, n_problems);
    if ~profile_options.(ProfileOptionKey.SILENT.value)
        fprintf('\nINFO: There are %d problems from "%s" to test.\n', n_problems, plib);
    end

    % Decide whether to delete the pool.
    try
        pool = gcp('nocreate');
        % We will delete the pool only when keep_pool is false and the number of workers is not equal to n_jobs.
        if (~isempty(pool) && (profile_options.n_jobs == 1 || profile_options.n_jobs ~= pool.NumWorkers)) && ~profile_options.(ProfileOptionKey.KEEP_POOL.value)
            if ~profile_options.(ProfileOptionKey.SILENT.value)
                delete(pool);
            else
                evalc("delete(pool)");
            end
        end
    catch
    end

    if (profile_options.(ProfileOptionKey.N_JOBS.value) == 1) || n_problems  == 1
        for i_problem = 1:n_problems
            problem_name = problem_names{i_problem};
            % Load the problem.
            try
                if ~profile_options.(ProfileOptionKey.SILENT.value)
                    loading_info = sprintf('\\nINFO: Loading problem   %%-%ds from \"%%s\".\\n', len_problem_names);
                    fprintf(loading_info, problem_name, plib);
                end
                problem = load(problem_name);
            catch
                if ~profile_options.(ProfileOptionKey.SILENT.value)
                    fail_load_info = sprintf('\\nINFO: Failed to load    %%-%ds from \"%%s\".\\n', len_problem_names);
                    fprintf(fail_load_info, problem_name, plib);
                end
                tmp_results{i_problem} = struct();
                continue;
            end
            result = solveOneProblem(solvers, problem, feature, problem_name, len_problem_names, profile_options, is_plot, path_hist_plots);
            tmp_results{i_problem} = result;
        end
    else
        try
            pool = gcp('nocreate');
        catch
            % If the user does not have the Parallel Computing Toolbox, the gcp function will not be available. We will print a error message and exit.
            error("The Parallel Computing Toolbox is not available. Please set the option `n_jobs` to 1.");
        end
        if isempty(pool)
            try
                if ~profile_options.(ProfileOptionKey.SILENT.value)
                    parpool(profile_options.(ProfileOptionKey.N_JOBS.value));
                else
                    evalc("parpool(profile_options.(ProfileOptionKey.N_JOBS.value))");
                end
            catch
                if ~profile_options.(ProfileOptionKey.SILENT.value)
                    fprintf("INFO: Failed to create a parallel pool with the specified number of workers. Using the default number of workers.\n");
                end
            end
        end
        parfor i_problem = 1:n_problems
            problem_name = problem_names{i_problem};
            % Load the problem.
            try
                if ~profile_options.(ProfileOptionKey.SILENT.value)
                    loading_info = sprintf('\\nINFO: Loading problem   %%-%ds from \"%%s\".\\n', len_problem_names);
                    fprintf(loading_info, problem_name, plib);
                end
                problem = load(problem_name);
            catch
                if ~profile_options.(ProfileOptionKey.SILENT.value)
                    fail_load_info = sprintf('\\nINFO: Failed to load    %%-%ds from \"%%s\".\\n', len_problem_names);
                    fprintf(fail_load_info, problem_name, plib);
                end
                tmp_results{i_problem} = struct();
                continue;
            end
            result = solveOneProblem(solvers, problem, feature, problem_name, len_problem_names, profile_options, is_plot, path_hist_plots);
            tmp_results{i_problem} = result;
        end
        if ~profile_options.(ProfileOptionKey.KEEP_POOL.value)
            if ~profile_options.(ProfileOptionKey.SILENT.value)
                delete(gcp);
                fprintf("INFO: Leaving the parallel section.\n");
            else
                evalc("delete(gcp)");
            end
        end
    end

    % Delete empty elements (struct) in tmp_results.
    tmp_results = tmp_results(~cellfun(@(x) isempty(fieldnames(x)), tmp_results));
    n_problems = length(tmp_results);

    if all(cellfun(@isempty, tmp_results))
        if ~profile_options.(ProfileOptionKey.SILENT.value)
            fprintf('INFO: All problems from "%s" are not solved successfully.\n', plib);
        end
        return;
    else
        % Process results.
        n_solvers = length(solvers);
        n_runs = feature.options.(FeatureOptionKey.N_RUNS.value);
        problem_types = cell(n_problems, 1);
        problem_dims = NaN(n_problems, 1);
        problem_mbs = NaN(n_problems, 1);
        problem_mlcons = NaN(n_problems, 1);
        problem_mnlcons = NaN(n_problems, 1);
        problem_mcons = NaN(n_problems, 1);
        fun_outs = NaN(n_problems, n_solvers, n_runs);
        maxcv_outs = NaN(n_problems, n_solvers, n_runs);
        fun_inits = NaN(n_problems, 1);
        maxcv_inits = NaN(n_problems, 1);
        n_evals = NaN(n_problems, n_solvers, n_runs);
        problem_names = cell(1, n_problems);
        computation_times = NaN(n_problems, n_solvers, n_runs);
        solvers_successes = NaN(n_problems, n_solvers, n_runs);
        for i_problem = 1:n_problems
            result = tmp_results{i_problem};
            problem_types{i_problem} = result.problem_type;
            problem_dims(i_problem) = result.problem_dim;
            problem_mbs(i_problem) = result.problem_mb;
            problem_mlcons(i_problem) = result.problem_mlcon;
            problem_mnlcons(i_problem) = result.problem_mnlcon;
            problem_mcons(i_problem) = result.problem_mcon;
            fun_outs(i_problem, :, :) = result.fun_out;
            maxcv_outs(i_problem, :, :) = result.maxcv_out;
            fun_inits(i_problem) = result.fun_init;
            maxcv_inits(i_problem) = result.maxcv_init;
            n_evals(i_problem, :, :) = result.n_eval;
            problem_names{i_problem} = result.problem_name;
            computation_times(i_problem, :, :) = result.computation_time;
            solvers_successes(i_problem, :, :) = result.solvers_success;
        end

        if n_problems > 0
            max_max_eval = ceil(max_eval_factor * max(problem_dims));
        else
            max_max_eval = 1;
        end
        fun_histories = NaN(n_problems, n_solvers, n_runs, max_max_eval);
        maxcv_histories = NaN(n_problems, n_solvers, n_runs, max_max_eval);
        for i_problem = 1:n_problems
            max_eval = ceil(max_eval_factor * problem_dims(i_problem));
            result = tmp_results{i_problem};
            fun_histories(i_problem, :, :, 1:max_eval) = result.fun_history;
            maxcv_histories(i_problem, :, :, 1:max_eval) = result.maxcv_history;
            if max_eval > 0
                fun_histories(i_problem, :, :, max_eval+1:end) = repmat(fun_histories(i_problem, :, :, max_eval), 1, 1, 1, size(fun_histories, 4) - max_eval);
                maxcv_histories(i_problem, :, :, max_eval+1:end) = repmat(maxcv_histories(i_problem, :, :, max_eval), 1, 1, 1, size(maxcv_histories, 4) - max_eval);
            end
        end
    end

    results.plib = plib;
    results.solver_names = profile_options.(ProfileOptionKey.SOLVER_NAMES.value);
    results.feature_stamp = profile_options.(ProfileOptionKey.FEATURE_STAMP.value);
    results.fun_histories = fun_histories;
    results.maxcv_histories = maxcv_histories;
    results.fun_outs = fun_outs;
    results.maxcv_outs = maxcv_outs;
    results.fun_inits = fun_inits;
    results.maxcv_inits = maxcv_inits;
    results.n_evals = n_evals;
    results.problem_names = problem_names;
    results.problem_types = problem_types;
    results.problem_dims = problem_dims;
    results.problem_mbs = problem_mbs;
    results.problem_mlcons = problem_mlcons;
    results.problem_mnlcons = problem_mnlcons;
    results.problem_mcons = problem_mcons;
    results.computation_times = computation_times;
    results.solvers_successes = solvers_successes;
end