function self_comparison()

    clc

    % This is a simple example to show how to test one solver with one
    % changeable hyperparameter.

    solvers = cell(1, 5);
    solver_names = cell(1, 5);
    for i_solver = 1:5
        solvers{i_solver} = @(fun, x0) fminsearch_test(fun, x0, 100 * i_solver);
        solver_names{i_solver} = sprintf('simplex_%d', 100 * i_solver);
    end
    options.solver_names = solver_names;
    benchmark(solvers, options)
end

function x = fminsearch_test(fun, x0, max_fun_evals)

    options = optimset('MaxFunEvals', max_fun_evals);
    x = fminsearch(fun, x0, options);
end