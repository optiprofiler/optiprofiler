function [fun_histories, maxcv_histories, fun_out, maxcv_out, fun_init, maxcv_init, n_eval, problem_name, problem_n, computation_time] = solveOneProblem(problem_name, solvers, labels, feature, custom_problem_loader, profile_options)
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

    if length(problem_name) == 2
        problem = custom_problem_loader(problem_name{2});
        problem_name = sprintf('%s (%s)', problem_name{1}, problem_name{2});
    else
        try
            problem = loader(problem_name);
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

    for i_solver = 1:n_solvers
        for i_run = 1:n_runs
            fprintf("Solving %s with %s (run %d/%d).\n", problem_name, labels{i_solver}, i_run, n_runs);
            time_start_solver_run = tic;
            % Construct featured_problem.
            featured_problem = FeaturedProblem(problem, feature, max_eval, i_run);
            if ~ismember(nargin(solvers{i_solver}), [2, 4, 8, 10])
                error("Solver %s has unknown signature.", labels{i_solver});
            end
            warning('off', 'all');
            try
                if nargin(solvers{i_solver}) == 2
                    [~, x] = evalc('solvers{i_solver}(@featured_problem.fun, featured_problem.x0)');
                elseif nargin(solvers{i_solver}) == 4
                    [~, x] = evalc('solvers{i_solver}(@featured_problem.fun, featured_problem.x0, featured_problem.xl, featured_problem.xu)');
                elseif nargin(solvers{i_solver}) == 8
                    [~, x] = evalc('solvers{i_solver}(@featured_problem.fun, featured_problem.x0, featured_problem.xl, featured_problem.xu, featured_problem.aub, featured_problem.bub, featured_problem.aeq, featured_problem.beq)');
                elseif nargin(solvers{i_solver}) == 10
                    [~, x] = evalc('solvers{i_solver}(@featured_problem.fun, featured_problem.x0, featured_problem.xl, featured_problem.xu, featured_problem.aub, featured_problem.bub, featured_problem.aeq, featured_problem.beq, @featured_problem.cub, @featured_problem.ceq)');
                end
                if strcmp(feature.name, FeatureName.PERMUTED.value)
                    [~, reverse_permutation] = sort(featured_problem.permutation);
                    x = x(reverse_permutation);
                end
                fun_out(i_solver, i_run) = problem.fun(x);
                maxcv_out(i_solver, i_run) = problem.maxcv(x);
                fprintf("Results for %s with %s (run %d/%d): f = %.4e, maxcv = %.4e (%.2f seconds).\n", problem_name, labels{i_solver}, i_run, n_runs, fun_out(i_solver, i_run), maxcv_out(i_solver, i_run), toc(time_start_solver_run));
            catch Exception
                fprintf("An error occurred while solving %s with %s: %s\n", problem_name, labels{i_solver}, Exception.message);
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
    end
    computation_time = toc(time_start);

end