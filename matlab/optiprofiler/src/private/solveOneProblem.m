function result = solveOneProblem(solvers, problem, feature, problem_name, len_problem_names, profile_options, is_plot, path_hist_plots)
%SOLVEONEPROBLEM solves one problem with all the solvers in solvers list.
%   Note that the input `problem` is a FeaturedProblem object when we call this function in `benchmark`.

    solver_names = profile_options.(ProfileOptionKey.SOLVER_NAMES.value);
    solver_names = cellfun(@(s) strrep(s, '\_', '_'), solver_names, 'UniformOutput', false);    % Remove backslash from the solver names.
    solver_isrand = profile_options.(ProfileOptionKey.SOLVER_ISRAND.value);
    result = struct();

    if ~isa(problem, 'Problem')
        error("MATLAB:solveOneProblem:invalid_problem", "The problem provided is not a valid Problem object.\n");
    end

    % Project the initial point if necessary.
    if profile_options.(ProfileOptionKey.PROJECT_X0.value)
        problem.project_x0;
    end

    % Solve the problem with each solver.
    n_solvers = length(solvers);
    n_runs = feature.options.(FeatureOptionKey.N_RUNS.value);
    max_eval = profile_options.(ProfileOptionKey.MAX_EVAL_FACTOR.value) * problem.n;
    max_eval = ceil(max_eval);  % Ensure max_eval is an integer greater than or equal to 1.
    n_eval = zeros(n_solvers, n_runs);
    fun_history = NaN(n_solvers, n_runs, max_eval);
    fun_out = NaN(n_solvers, n_runs);
    maxcv_history = NaN(n_solvers, n_runs, max_eval);
    maxcv_out = NaN(n_solvers, n_runs);
    fun_inits = NaN(n_runs, 1);
    maxcv_inits = NaN(n_runs, 1);
    computation_time = NaN(n_solvers, n_runs);
    solvers_success = false(n_solvers, n_runs);

    % The number of real runs for each solver, which is determined by feature and solver_isrand.
    real_n_runs = n_runs * ones(n_solvers, 1);
    for i_solver = 1:n_solvers
        % If the solver is deterministic and the feature is deterministic, then the number of runs is 1.
        if ~feature.is_stochastic && ~solver_isrand(i_solver)
            real_n_runs(i_solver) = 1;
        end
    end

    len_solver_names = max(cellfun(@length, solver_names));

    for i_solver = 1:n_solvers
        for i_run = 1:real_n_runs(i_solver)
            if ~profile_options.(ProfileOptionKey.SILENT.value)
                format_info_start = sprintf("\\nINFO: Start  solving    %%-%ds with %%-%ds (run %%2d/%%2d).\\n", len_problem_names, len_solver_names);
                fprintf(format_info_start, problem_name, solver_names{i_solver}, i_run, real_n_runs(i_solver));
            end

            % Construct featured_problem.
            real_seed = mod(23333 * profile_options.(ProfileOptionKey.SEED.value) + 211 * i_run, 2^32);
            featured_problem = FeaturedProblem(problem, feature, max_eval, real_seed);

            % Save the problem information.
            problem_type = featured_problem.ptype;
            problem_dim = featured_problem.n;
            problem_mb = featured_problem.mb;
            problem_mlcon = featured_problem.mlcon;
            problem_mnlcon = featured_problem.mnlcon;
            problem_mcon = featured_problem.mcon;
            fun_inits(i_run) = featured_problem.fun_init;
            maxcv_inits(i_run) = featured_problem.maxcv_init;

            % Solve the problem with the solver.
            warning('off', 'all');
            time_start_solver_run = tic;
            try
                if profile_options.(ProfileOptionKey.SOLVER_VERBOSE.value) == 2
                    switch problem_type
                        case 'u'
                            x = solvers{i_solver}(@(x) featured_problem.fun(x), featured_problem.x0);
                        case 'b'
                            x = solvers{i_solver}(@(x) featured_problem.fun(x), featured_problem.x0, featured_problem.xl, featured_problem.xu);
                        case 'l'
                            x = solvers{i_solver}(@(x) featured_problem.fun(x), featured_problem.x0, featured_problem.xl, featured_problem.xu, featured_problem.aub, featured_problem.bub, featured_problem.aeq, featured_problem.beq);
                        case 'n'
                            x = solvers{i_solver}(@(x) featured_problem.fun(x), featured_problem.x0, featured_problem.xl, featured_problem.xu, featured_problem.aub, featured_problem.bub, featured_problem.aeq, featured_problem.beq, @(x) featured_problem.cub(x), @(x) featured_problem.ceq(x));
                    end
                else
                    switch problem_type
                        case 'u'
                            [~, x] = evalc('solvers{i_solver}(@(x) featured_problem.fun(x), featured_problem.x0)');
                        case 'b'
                            [~, x] = evalc('solvers{i_solver}(@(x) featured_problem.fun(x), featured_problem.x0, featured_problem.xl, featured_problem.xu)');
                        case 'l'
                            [~, x] = evalc('solvers{i_solver}(@(x) featured_problem.fun(x), featured_problem.x0, featured_problem.xl, featured_problem.xu, featured_problem.aub, featured_problem.bub, featured_problem.aeq, featured_problem.beq)');
                        case 'n'
                            [~, x] = evalc('solvers{i_solver}(@(x) featured_problem.fun(x), featured_problem.x0, featured_problem.xl, featured_problem.xu, featured_problem.aub, featured_problem.bub, featured_problem.aeq, featured_problem.beq, @(x) featured_problem.cub(x), @(x) featured_problem.ceq(x))');
                    end
                end

                computation_time(i_solver, i_run) = toc(time_start_solver_run);

                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                % It is very important to transform the solution back to the one related to the original problem. (Note that the problem we solve has the objective function f(A @ x + b). Thus, if x is the output solution, then A @ x + b is the solution of the original problem.)
                [A, b] = featured_problem.feature.modifier_affine(featured_problem.seed, featured_problem.problem);
                x = A * x + b;
                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                
                % Use problem.fun and problem.maxcv to evaluate the solution since it is possible that featured_problem.fun and featured_problem.maxcv are modified.
                fun_out(i_solver, i_run) = problem.fun(x);
                maxcv_out(i_solver, i_run) = problem.maxcv(x);
                % Calculate the minimum function value and the minimum constraint violation, omitting the NaN values.
                fun_min = min(featured_problem.fun_hist, [], 'omitnan');
                maxcv_min = min(featured_problem.maxcv_hist, [], 'omitnan');
                if ~profile_options.(ProfileOptionKey.SILENT.value)
                    format_info_end = sprintf("INFO: Finish solving    %%-%ds with %%-%ds (run %%2d/%%2d) in %%.2f seconds.\\n", len_problem_names, len_solver_names);
                    fprintf(format_info_end, problem_name, solver_names{i_solver}, i_run, real_n_runs(i_solver), computation_time(i_solver, i_run));
                    switch problem_type
                        case 'u'
                            format_info_output = sprintf("INFO: Output result for %%-%ds with %%-%ds (run %%2d/%%2d): f = %%10.4e.\\n", len_problem_names, len_solver_names);
                            fprintf(format_info_output, problem_name, solver_names{i_solver}, i_run, real_n_runs(i_solver), fun_out(i_solver, i_run));
                            format_info_best = sprintf("INFO: Best   result for %%-%ds with %%-%ds (run %%2d/%%2d): f = %%10.4e.\\n", len_problem_names, len_solver_names);
                            fprintf(format_info_best, problem_name, solver_names{i_solver}, i_run, real_n_runs(i_solver), fun_min);
                        otherwise
                            format_info_output = sprintf("INFO: Output result for %%-%ds with %%-%ds (run %%2d/%%2d): f = %%10.4e, maxcv = %%10.4e.\\n", len_problem_names, len_solver_names);
                            fprintf(format_info_output, problem_name, solver_names{i_solver}, i_run, real_n_runs(i_solver), fun_out(i_solver, i_run), maxcv_out(i_solver, i_run));
                            format_info_best = sprintf("INFO: Best   result for %%-%ds with %%-%ds (run %%2d/%%2d): f = %%10.4e, maxcv = %%10.4e.\\n", len_problem_names, len_solver_names);
                            fprintf(format_info_best, problem_name, solver_names{i_solver}, i_run, real_n_runs(i_solver), fun_min, maxcv_min);
                    end
                end
            catch Exception
                if profile_options.(ProfileOptionKey.SOLVER_VERBOSE.value) ~= 0
                    fprintf("INFO: An error occurred while solving %s with %s (run %d/%d): %s\n", problem_name, solver_names{i_solver}, i_run, real_n_runs(i_solver), Exception.message);
                end
            end
            warning('on', 'all');
            n_eval(i_solver, i_run) = featured_problem.n_eval_fun;
            fun_history(i_solver, i_run, 1:n_eval(i_solver, i_run)) = featured_problem.fun_hist(1:n_eval(i_solver, i_run));
            maxcv_history(i_solver, i_run, 1:n_eval(i_solver, i_run)) = featured_problem.maxcv_hist(1:n_eval(i_solver, i_run));
            if n_eval(i_solver, i_run) > 0
                fun_history(i_solver, i_run, n_eval(i_solver, i_run)+1:end) = fun_history(i_solver, i_run, n_eval(i_solver, i_run));
                maxcv_history(i_solver, i_run, n_eval(i_solver, i_run)+1:end) = maxcv_history(i_solver, i_run, n_eval(i_solver, i_run));
                solvers_success(i_solver, i_run) = true;
            end
        end
        % If real_n_runs(i_solver) == 1 ~= n_runs, then we need to copy the result to the other runs.
        if real_n_runs(i_solver) == 1 && n_runs > 1
            for j = 2:n_runs
                n_eval(i_solver, j) = n_eval(i_solver, 1);
                fun_history(i_solver, j, :) = fun_history(i_solver, 1, :);
                maxcv_history(i_solver, j, :) = maxcv_history(i_solver, 1, :);
                fun_out(i_solver, j) = fun_out(i_solver, 1);
                maxcv_out(i_solver, j) = maxcv_out(i_solver, 1);
            end
        end
    end
    fprintf("\n");

    % Complete fun_inits and maxcv_inits if real_n_runs(i_solver) < n_runs for all solvers.
    % In this case, real_n_runs(i_solver) should be 1 for all solvers.
    if max(real_n_runs) < n_runs
        for i_run = 2:n_runs
            fun_inits(i_run) = fun_inits(1);
            maxcv_inits(i_run) = maxcv_inits(1);
        end
    end

    % Return the result.
    result.fun_history = fun_history;
    result.maxcv_history = maxcv_history;
    result.fun_out = fun_out;
    result.maxcv_out = maxcv_out;
    result.fun_inits = fun_inits;
    result.maxcv_inits = maxcv_inits;
    result.n_eval = n_eval;
    result.problem_name = problem_name;
    result.problem_type = problem_type;
    result.problem_dim = problem_dim;
    result.problem_mb = problem_mb;
    result.problem_mlcon = problem_mlcon;
    result.problem_mnlcon = problem_mnlcon;
    result.problem_mcon = problem_mcon;
    result.computation_time = computation_time;
    result.solvers_success = solvers_success;

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%% History plots of the computation %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    if is_plot && strcmp(profile_options.(ProfileOptionKey.DRAW_HIST_PLOTS.value), 'parallel')
        exportHist(problem_name, problem_type, problem_dim, solver_names, solvers_success, fun_history, maxcv_history, fun_inits, maxcv_inits, n_eval, profile_options, path_hist_plots)
    end
end