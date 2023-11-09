function [fun_values, maxcv_values, problem_name, problem_n] = solveOne(problem_name, solvers, feature, max_eval_factor)
    %SOLVEONE is defined to solve one test problem with a set of solvers.
    
    % Try to load the test problem.
    try
        problem = loadTestProblem(problem_name);
    catch ME
        fun_values = [];
        maxcv_values = [];
        problem_n = [];
       return; 
    end

    problem_n = problem.n;

    % Solve the problem with each solver.
    n_solvers = length(solvers);
    n_runs = feature.options.(FeatureOptionKey.N_RUNS.value);
    max_eval = max_eval_factor * problem.n;
    fun_values = NaN(n_solvers, n_runs, max_eval);
    maxcv_values = NaN(n_solvers, n_runs, max_eval);

    for i_solver = 1:n_solvers
        for i_run = 1:n_runs
            % clear featured_problem
            featured_problem = FeaturedProblem(problem, feature);
            warning('off', 'all');
            try
                [~] = solvers{i_solver}(@(x) featured_problem.fun(x, i_run), featured_problem.x0, featured_problem.xl, featured_problem.xu, featured_problem.aub, featured_problem.bub, featured_problem.aeq, featured_problem.beq, @featured_problem.cub, @featured_problem.ceq, max_eval);
            catch ME
                % Ignore all the errors.
            end
            warning('on', 'all');
            n_eval = min(featured_problem.n_eval, max_eval);
            fun_values(i_solver, i_run, 1:n_eval) = featured_problem.fun_values(1:n_eval);
            maxcv_values(i_solver, i_run, 1:n_eval) = featured_problem.maxcv_values(1:n_eval);
            if n_eval > 0
                fun_values(i_solver, i_run, n_eval+1:end) = fun_values(i_solver, i_run, n_eval);
                maxcv_values(i_solver, i_run, n_eval+1:end) = maxcv_values(i_solver, i_run, n_eval);
            end
        end
    end

end