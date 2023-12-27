function [fun_values, maxcv_values, fun_init, maxcv_init, n_eval, problem_name, problem_n] = solveOne(problem_name, problem_options, solvers, labels, feature, max_eval_factor)
    %SOLVEONE is defined to solve one test problem with a set of solvers.
    
    fun_values = [];
    maxcv_values = [];
    fun_init = [];
    maxcv_init = [];
    n_eval = [];
    problem_n = [];
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
    n_solvers = length(solvers);
    n_runs = feature.options.(FeatureOptionKey.N_RUNS.value);
    max_eval = max_eval_factor * problem.n;
    n_eval = zeros(n_solvers, n_runs);
    fun_values = NaN(n_solvers, n_runs, max_eval);
    maxcv_values = NaN(n_solvers, n_runs, max_eval);

    for i_solver = 1:n_solvers
        for i_run = 1:n_runs
            % clear featured_problem
            featured_problem = FeaturedProblem(problem, feature);
            warning('off', 'all');
            try
                [~] = solvers{i_solver}(@(x) featured_problem.fun(x, i_run), featured_problem.x0, featured_problem.xl, featured_problem.xu, featured_problem.aub, featured_problem.bub, featured_problem.aeq, featured_problem.beq, @featured_problem.cub, @featured_problem.ceq, max_eval);
            catch
                % Ignore all the errors.
            end
            warning('on', 'all');
            n_eval(i_solver, i_run) = min(featured_problem.n_eval, max_eval);
            fun_values(i_solver, i_run, 1:n_eval(i_solver, i_run)) = featured_problem.fun_values(1:n_eval(i_solver, i_run));
            maxcv_values(i_solver, i_run, 1:n_eval(i_solver, i_run)) = featured_problem.maxcv_values(1:n_eval(i_solver, i_run));
            if n_eval(i_solver, i_run) > 0
                fun_values(i_solver, i_run, n_eval(i_solver, i_run)+1:end) = fun_values(i_solver, i_run, n_eval(i_solver, i_run));
                maxcv_values(i_solver, i_run, n_eval(i_solver, i_run)+1:end) = maxcv_values(i_solver, i_run, n_eval(i_solver, i_run));
            end
        end
    end

end