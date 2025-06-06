function result = solveOneProblem(solvers, problem, feature, problem_name, len_problem_names, profile_options, is_plot, path_hist_plots)
%SOLVEONEPROBLEM solves one problem with all the solvers in solvers list.

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

    % Evaluate the functions at the initial point.
    fun_init = problem.fun(problem.x0);
    maxcv_init = problem.maxcv(problem.x0);

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
            time_start_solver_run = tic;

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

            % Solve the problem with the solver.
            warning('off', 'all');
            try
                switch problem.ptype
                    case 'u'
                        if profile_options.(ProfileOptionKey.SOLVER_VERBOSE.value) == 2
                            x = solvers{i_solver}(@(x) featured_problem.fun(x), featured_problem.x0);
                        else
                            [~, x] = evalc('solvers{i_solver}(@(x) featured_problem.fun(x), featured_problem.x0)');
                        end
                    case 'b'
                        if profile_options.(ProfileOptionKey.SOLVER_VERBOSE.value) == 2
                            x = solvers{i_solver}(@(x) featured_problem.fun(x), featured_problem.x0, featured_problem.xl, featured_problem.xu);
                        else
                            [~, x] = evalc('solvers{i_solver}(@(x) featured_problem.fun(x), featured_problem.x0, featured_problem.xl, featured_problem.xu)');
                        end
                    case 'l'
                        if profile_options.(ProfileOptionKey.SOLVER_VERBOSE.value) == 2
                            x = solvers{i_solver}(@(x) featured_problem.fun(x), featured_problem.x0, featured_problem.xl, featured_problem.xu, featured_problem.aub, featured_problem.bub, featured_problem.aeq, featured_problem.beq);
                        else
                            [~, x] = evalc('solvers{i_solver}(@(x) featured_problem.fun(x), featured_problem.x0, featured_problem.xl, featured_problem.xu, featured_problem.aub, featured_problem.bub, featured_problem.aeq, featured_problem.beq)');
                        end
                    case 'n'
                        if profile_options.(ProfileOptionKey.SOLVER_VERBOSE.value) == 2
                            x = solvers{i_solver}(@(x) featured_problem.fun(x), featured_problem.x0, featured_problem.xl, featured_problem.xu, featured_problem.aub, featured_problem.bub, featured_problem.aeq, featured_problem.beq, @(x) featured_problem.cub(x), @(x) featured_problem.ceq(x));
                        else
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
                    format_info_end = sprintf("INFO: Finish solving    %%-%ds with %%-%ds (run %%2d/%%2d) (in %%.2f seconds).\\n", len_problem_names, len_solver_names);
                    fprintf(format_info_end, problem_name, solver_names{i_solver}, i_run, real_n_runs(i_solver), computation_time(i_solver, i_run));
                    switch problem.ptype
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
            n_eval(i_solver, i_run) = featured_problem.n_eval;
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

    % Return the result.
    result.fun_history = fun_history;
    result.maxcv_history = maxcv_history;
    result.fun_out = fun_out;
    result.maxcv_out = maxcv_out;
    result.fun_init = fun_init;
    result.maxcv_init = maxcv_init;
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

    if ~is_plot || profile_options.(ProfileOptionKey.SCORE_ONLY.value) || all(~solvers_success(:))
        return;
    end

    try
        merit_fun = profile_options.(ProfileOptionKey.MERIT_FUN.value);
        try
            merit_history = meritFunCompute(merit_fun, fun_history, maxcv_history, maxcv_init);
            merit_init = meritFunCompute(merit_fun, fun_init, maxcv_init, maxcv_init);
        catch
            error("MATLAB:solveOneProblem:merit_fun_error", "Error occurred while calculating the merit values. Please check the merit function.");
        end

        % Create the figure for the summary.
        warning('off');
        if strcmp(problem.ptype, 'u')
            n_cols = 1;
        else
            n_cols = 3;
        end
        defaultFigurePosition = get(0, 'DefaultFigurePosition');
        default_width = defaultFigurePosition(3);
        default_height = defaultFigurePosition(4);
        fig_summary = figure('Position', [defaultFigurePosition(1:2), n_cols * default_width, 2 * default_height], 'visible', 'off');
        T_summary = tiledlayout(fig_summary, 2, 1, 'Padding', 'compact', 'TileSpacing', 'compact');
        F_title = strrep(profile_options.feature_stamp, '_', '\_');
        P_title = strrep(problem_name, '_', '\_');
        T_title = ['Solving ``', P_title, '" with ``', F_title, '" feature'];
        title_fontsize = min(12, 1.2 * default_width / length(T_title));
        title(T_summary, ['Solving ``', P_title, '" with ``', F_title, '" feature'], 'Interpreter', 'latex', 'FontSize', title_fontsize);
        % Use gobjects to create arrays of handles and axes.
        t_summary = gobjects(2, 1);
        axs_summary = gobjects([2, 1, 1, n_cols]);
        i_axs = 0;
        for i = 1:2
            t_summary(i) = tiledlayout(T_summary, 1, n_cols, 'Padding', 'compact', 'TileSpacing', 'compact');
            t_summary(i).Layout.Tile = i;
            for j = 1:n_cols
                i_axs = i_axs + 1;
                axs_summary(i_axs) = nexttile(t_summary(i));
            end
        end
        ylabel(t_summary(1), "History profiles", 'Interpreter', 'latex', 'FontSize', 14);
        ylabel(t_summary(2), "Cummin history profiles", 'Interpreter', 'latex', 'FontSize', 14);

        if strcmp(problem.ptype, 'u')
            cell_axs_summary = {axs_summary(1)};
            cell_axs_summary_cum = {axs_summary(2)};
        else
            cell_axs_summary = {axs_summary(1), axs_summary(2), axs_summary(3)};
            cell_axs_summary_cum = {axs_summary(4), axs_summary(5), axs_summary(6)};
        end

        pdf_hist_file_name = [regexprep(regexprep(regexprep(strrep(problem_name,' ','_'),'[^a-zA-Z0-9\-_]',''),'[-_]+','_'),'^[-_]+',''), '.pdf'];
        pdf_summary = fullfile(path_hist_plots, pdf_hist_file_name);
        processed_solver_names = cellfun(@(s) strrep(s, '_', '\_'), solver_names, 'UniformOutput', false);

        drawHist(fun_history, maxcv_history, merit_history, fun_init, maxcv_init, merit_init, processed_solver_names, cell_axs_summary, false, problem.ptype, problem_dim, n_eval, profile_options, default_height);
        drawHist(fun_history, maxcv_history, merit_history, fun_init, maxcv_init, merit_init, processed_solver_names, cell_axs_summary_cum, true, problem.ptype, problem_dim, n_eval, profile_options, default_height);

        exportgraphics(fig_summary, pdf_summary, 'ContentType', 'vector');
        warning('on');
        close(fig_summary);
    catch Exception
        if ~profile_options.(ProfileOptionKey.SILENT.value)
            fprintf("INFO: An error occurred while plotting the history plots of the problem %s: %s\n", problem_name, Exception.message);
        end
    end
end