function [fun_histories, maxcv_histories, fun_out, maxcv_out, fun_init, maxcv_init, n_eval, problem_name, problem_n, computation_time] = solveOneProblem(problem_name, solvers, labels, feature, max_eval_factor, problem_options)
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
    load_success = false;

    % Try to load the test problem from CUTEst.
    try
        problem = loadCutest(problem_name, problem_options);
        problem_n = problem.n;
        load_success = true;
    catch
    end

    % Try to load the test problem.
    try
        problem = loadTestProblem(problem_name);
        problem_n = problem.n;
        load_success = true;
    catch
    end

    if ~load_success
        return;
    end
    

    % Evaluate the functions at the initial point.
    fun_init = problem.fun(problem.x0);
    maxcv_init = problem.maxcv(problem.x0);

    % Solve the problem with each solver.
    time_start = tic;
    n_solvers = length(solvers);
    n_runs = feature.options.(FeatureOptionKey.N_RUNS.value);
    max_eval = max_eval_factor * problem.n;
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
            warning('off', 'all');
            try
                [~, x] = evalc('solvers{i_solver}(@(x) featured_problem.fun(x), featured_problem.x0, featured_problem.xl, featured_problem.xu, featured_problem.aub, featured_problem.bub, featured_problem.aeq, featured_problem.beq, @featured_problem.cub, @featured_problem.ceq, max_eval)');
                fun_out(i_solver, i_run) = problem.fun(x);
                maxcv_out(i_solver, i_run) = problem.maxcv(x);
                fprintf("Results for %s with %s (run %d/%d): f = %.4e, maxcv = %.4e (%.2f seconds).\n", problem_name, labels{i_solver}, i_run, n_runs, fun_out(i_solver, i_run), maxcv_out(i_solver, i_run), toc(time_start_solver_run));
            catch Exception
                fprintf("An error occurred while solving %s with %s: %s.\n", problem_name, labels{i_solver}, Exception.message);
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