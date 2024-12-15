function [fun_histories, maxcv_histories, fun_out, maxcv_out, fun_init, maxcv_init, n_eval, problem_name, problem_n, computation_time, solvers_success] = solveOneProblem(problem_name, solvers, solver_names, solver_isrand, feature, len_problem_names, custom_problem_loader, profile_options, is_plot, path_hist_plots)
%SOLVEONEPROBLEM solves one problem with all the solvers in solvers list.

    fun_histories = [];
    maxcv_histories = [];
    fun_out = [];
    maxcv_out = [];
    fun_init = [];
    maxcv_init = [];
    n_eval = [];
    problem_n = [];
    computation_time = [];
    solvers_success = false;

    if length(problem_name) == 2
        problem = custom_problem_loader(problem_name{2});

        % Verify whether problem is a Problem object.
        if ~isa(problem, 'Problem')
            fprintf("INFO: Custom problem %s cannot be loaded by custom_problem_loader.\n", problem_name{1});
            return;
        end
        problem_name = sprintf('%s %s', problem_name{1}, problem_name{2});
    else
        try
            problem = s_load(problem_name);
        catch
            return;
        end
    end

    problem_n = problem.n;

    % Project the initial point if necessary.
    if profile_options.(ProfileOptionKey.PROJECT_X0.value)
        problem.project_x0;
    end

    % Evaluate the functions at the initial point.
    fun_init = problem.fun(problem.x0);
    maxcv_init = problem.maxcv(problem.x0);

    % Solve the problem with each solver.
    time_start = tic;
    n_solvers = length(solvers);
    n_runs = feature.options.(FeatureOptionKey.N_RUNS.value);
    max_eval = profile_options.(ProfileOptionKey.MAX_EVAL_FACTOR.value) * problem.n;
    n_eval = zeros(n_solvers, n_runs);
    fun_histories = NaN(n_solvers, n_runs, max_eval);
    fun_out = NaN(n_solvers, n_runs);
    maxcv_histories = NaN(n_solvers, n_runs, max_eval);
    maxcv_out = NaN(n_solvers, n_runs);

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
                format_info_start = sprintf("INFO: Solving %%-%ds with %%-%ds (run %%2d/%%2d).\\n\\n", len_problem_names, len_solver_names);
                fprintf(format_info_start, problem_name, solver_names{i_solver}, i_run, real_n_runs(i_solver));
            end
            time_start_solver_run = tic;
            % Construct featured_problem.
            real_seed = mod(23333 * profile_options.(ProfileOptionKey.SEED.value) + 211 * i_run, 2^32);
            featured_problem = FeaturedProblem(problem, feature, max_eval, real_seed);
            warning('off', 'all');
            try
                switch problem.p_type
                    case 'unconstrained'
                        if profile_options.(ProfileOptionKey.SOLVER_VERBOSE.value) == 2
                            if solver_isrand(i_solver)
                                x = solvers{i_solver}(@(x) featured_problem.fun(x), featured_problem.x0, real_seed);
                            else
                                x = solvers{i_solver}(@(x) featured_problem.fun(x), featured_problem.x0);
                            end
                        else
                            if solver_isrand(i_solver)
                                [~, x] = evalc('solvers{i_solver}(@(x) featured_problem.fun(x), featured_problem.x0, real_seed)');
                            else
                                [~, x] = evalc('solvers{i_solver}(@(x) featured_problem.fun(x), featured_problem.x0)');
                            end
                        end
                    case 'bound-constrained'
                        if profile_options.(ProfileOptionKey.SOLVER_VERBOSE.value) == 2
                            if solver_isrand(i_solver)
                                x = solvers{i_solver}(@(x) featured_problem.fun(x), featured_problem.x0, featured_problem.xl, featured_problem.xu, real_seed);
                            else
                                x = solvers{i_solver}(@(x) featured_problem.fun(x), featured_problem.x0, featured_problem.xl, featured_problem.xu);
                            end
                        else
                            if solver_isrand(i_solver)
                                [~, x] = evalc('solvers{i_solver}(@(x) featured_problem.fun(x), featured_problem.x0, featured_problem.xl, featured_problem.xu, real_seed)');
                            else
                                [~, x] = evalc('solvers{i_solver}(@(x) featured_problem.fun(x), featured_problem.x0, featured_problem.xl, featured_problem.xu)');
                            end
                        end
                    case 'linearly constrained'
                        if profile_options.(ProfileOptionKey.SOLVER_VERBOSE.value) == 2
                            if solver_isrand(i_solver)
                                x = solvers{i_solver}(@(x) featured_problem.fun(x), featured_problem.x0, featured_problem.xl, featured_problem.xu, featured_problem.aub, featured_problem.bub, featured_problem.aeq, featured_problem.beq, real_seed);
                            else
                                x = solvers{i_solver}(@(x) featured_problem.fun(x), featured_problem.x0, featured_problem.xl, featured_problem.xu, featured_problem.aub, featured_problem.bub, featured_problem.aeq, featured_problem.beq);
                            end
                        else
                            if solver_isrand(i_solver)
                                [~, x] = evalc('solvers{i_solver}(@(x) featured_problem.fun(x), featured_problem.x0, featured_problem.xl, featured_problem.xu, featured_problem.aub, featured_problem.bub, featured_problem.aeq, featured_problem.beq, real_seed)');
                            else
                                [~, x] = evalc('solvers{i_solver}(@(x) featured_problem.fun(x), featured_problem.x0, featured_problem.xl, featured_problem.xu, featured_problem.aub, featured_problem.bub, featured_problem.aeq, featured_problem.beq)');
                            end
                        end
                    case 'nonlinearly constrained'
                        if profile_options.(ProfileOptionKey.SOLVER_VERBOSE.value) == 2
                            if solver_isrand(i_solver)
                                x = solvers{i_solver}(@(x) featured_problem.fun(x), featured_problem.x0, featured_problem.xl, featured_problem.xu, featured_problem.aub, featured_problem.bub, featured_problem.aeq, featured_problem.beq, @(x) featured_problem.cub(x), @(x) featured_problem.ceq(x), real_seed);
                            else
                                x = solvers{i_solver}(@(x) featured_problem.fun(x), featured_problem.x0, featured_problem.xl, featured_problem.xu, featured_problem.aub, featured_problem.bub, featured_problem.aeq, featured_problem.beq, @(x) featured_problem.cub(x), @(x) featured_problem.ceq(x));
                            end
                        else
                            if solver_isrand(i_solver)
                                [~, x] = evalc('solvers{i_solver}(@(x) featured_problem.fun(x), featured_problem.x0, featured_problem.xl, featured_problem.xu, featured_problem.aub, featured_problem.bub, featured_problem.aeq, featured_problem.beq, @(x) featured_problem.cub(x), @(x) featured_problem.ceq(x), real_seed)');
                            else
                                [~, x] = evalc('solvers{i_solver}(@(x) featured_problem.fun(x), featured_problem.x0, featured_problem.xl, featured_problem.xu, featured_problem.aub, featured_problem.bub, featured_problem.aeq, featured_problem.beq, @(x) featured_problem.cub(x), @(x) featured_problem.ceq(x))');
                            end
                        end
                end

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
                    format_info_end = sprintf("INFO: Finish solving     %%-%ds with %%-%ds (run %%2d/%%2d) (in %%.2f seconds).\\n", len_problem_names, len_solver_names);
                    fprintf(format_info_end, problem_name, solver_names{i_solver}, i_run, real_n_runs(i_solver), toc(time_start_solver_run));
                    switch problem.p_type
                        case 'unconstrained'
                            format_info_output = sprintf("INFO: Output results for %%-%ds with %%-%ds (run %%2d/%%2d): f = %%10.4e.\\n", len_problem_names, len_solver_names);
                            fprintf(format_info_output, problem_name, solver_names{i_solver}, i_run, real_n_runs(i_solver), fun_out(i_solver, i_run));
                            format_info_best = sprintf("INFO: Best   results for %%-%ds with %%-%ds (run %%2d/%%2d): f = %%10.4e.\\n", len_problem_names, len_solver_names);
                            fprintf(format_info_best, problem_name, solver_names{i_solver}, i_run, real_n_runs(i_solver), fun_min);
                        otherwise
                            format_info_output = sprintf("INFO: Output results for %%-%ds with %%-%ds (run %%2d/%%2d): f = %%10.4e, maxcv = %%10.4e.\\n", len_problem_names, len_solver_names);
                            fprintf(format_info_output, problem_name, solver_names{i_solver}, i_run, real_n_runs(i_solver), fun_out(i_solver, i_run), maxcv_out(i_solver, i_run));
                            format_info_best = sprintf("INFO: Best   results for %%-%ds with %%-%ds (run %%2d/%%2d): f = %%10.4e, maxcv = %%10.4e.\\n", len_problem_names, len_solver_names);
                            fprintf(format_info_best, problem_name, solver_names{i_solver}, i_run, real_n_runs(i_solver), fun_min, maxcv_min);
                    end
                end
                % If one solver succeeds once, then we say this problem is solved.
                solvers_success = true;
            catch Exception
                if profile_options.(ProfileOptionKey.SOLVER_VERBOSE.value) ~= 0
                    fprintf("INFO: An error occurred while solving %s with %s (run %d/%d): %s\n", problem_name, solver_names{i_solver}, i_run, real_n_runs(i_solver), Exception.message);
                end
            end
            warning('on', 'all');
            n_eval(i_solver, i_run) = featured_problem.n_eval;
            fun_histories(i_solver, i_run, 1:n_eval(i_solver, i_run)) = featured_problem.fun_hist(1:n_eval(i_solver, i_run));
            maxcv_histories(i_solver, i_run, 1:n_eval(i_solver, i_run)) = featured_problem.maxcv_hist(1:n_eval(i_solver, i_run));
            if n_eval(i_solver, i_run) > 0
                fun_histories(i_solver, i_run, n_eval(i_solver, i_run)+1:end) = fun_histories(i_solver, i_run, n_eval(i_solver, i_run));
                maxcv_histories(i_solver, i_run, n_eval(i_solver, i_run)+1:end) = maxcv_histories(i_solver, i_run, n_eval(i_solver, i_run));
            end
        end
        % If real_n_runs(i_solver) == 1 ~= n_runs, then we need to copy the results to the other runs.
        if real_n_runs(i_solver) == 1 && n_runs > 1
            n_eval(i_solver, 2:n_runs) = n_eval(i_solver, 1);
            fun_histories(i_solver, 2:n_runs, :) = fun_histories(i_solver, 1, :);
            maxcv_histories(i_solver, 2:n_runs, :) = maxcv_histories(i_solver, 1, :);
            fun_out(i_solver, 2:n_runs) = fun_out(i_solver, 1);
            maxcv_out(i_solver, 2:n_runs) = maxcv_out(i_solver, 1);
        end
    end
    computation_time = toc(time_start);

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%% History plots of the computation %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    if ~is_plot || ~solvers_success
        return;
    end

    try
        merit_histories = computeMeritValues(fun_histories, maxcv_histories, maxcv_init);
        merit_init = computeMeritValues(fun_init, maxcv_init, maxcv_init);

        % Create the figure for the summary.
        warning('off');
        if strcmp(problem.p_type, 'unconstrained')
            n_cols = 1;
        else
            n_cols = 3;
        end
        defaultFigurePosition = get(0, 'DefaultFigurePosition');
        default_width = defaultFigurePosition(3);
        default_height = defaultFigurePosition(4);
        fig_summary = figure('Position', [defaultFigurePosition(1:2), n_cols * default_width, 2 * default_height], ...
        'visible', 'off');
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

        if strcmp(problem.p_type, 'unconstrained')
            cell_axs_summary = {axs_summary(1)};
            cell_axs_summary_cum = {axs_summary(2)};
        else
            cell_axs_summary = {axs_summary(1), axs_summary(2), axs_summary(3)};
            cell_axs_summary_cum = {axs_summary(4), axs_summary(5), axs_summary(6)};
        end

        hist_file_name = [regexprep(regexprep(regexprep(strrep(problem_name,' ','_'),'[^a-zA-Z0-9\-_]',''),'[-_]+','_'),'^[-_]+',''), '.pdf'];
        pdf_summary = fullfile(path_hist_plots, hist_file_name);
        processed_solver_names = cellfun(@(s) strrep(s, '_', '\_'), solver_names, 'UniformOutput', false);

        drawHist(fun_histories, maxcv_histories, merit_histories, fun_init, maxcv_init, merit_init, processed_solver_names, cell_axs_summary, false, problem.p_type, problem_n, n_eval, profile_options, default_height);
        drawHist(fun_histories, maxcv_histories, merit_histories, fun_init, maxcv_init, merit_init, processed_solver_names, cell_axs_summary_cum, true, problem.p_type, problem_n, n_eval, profile_options, default_height);

        exportgraphics(fig_summary, pdf_summary, 'ContentType', 'vector');
        warning('on');
        close(fig_summary);
    catch Exception
        fprintf("INFO: An error occurred while plotting the history plots of the problem %s: %s\n", problem_name, Exception.message);
    end

end