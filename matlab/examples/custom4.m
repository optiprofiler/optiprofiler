function custom4()

    clc

    % This is a simple example to show how to test one solver with one
    % changeable hyperparameter.

    n_solvers = 5;
    solvers = cell(1, n_solvers);
    solver_names = cell(1, n_solvers);
    for i_solver = 1:n_solvers
        solvers{i_solver} = @(fun, x0) test_solver(fun, x0, exp(-0.5 * i_solver));
        solver_names{i_solver} = sprintf('solver_%d', i_solver);
    end
    options.solver_names = solver_names;
    benchmark(solvers, options)
end

function x = test_solver(fun, x0, theta)

    x = x0;
    f = fun(x);
    alpha = 1;
    randstream = RandStream('mt19937ar', 'Seed', 2);
    for i = 1:500 * length(x0)
        d = randstream.randn(size(x0));
        d = d / norm(d);
        f_new = fun(x + alpha * d);
        if f - f_new >= 1e-3 * alpha^2
            x = x + alpha * d;
            f = f_new;
        else
            f_new = fun(x - alpha * d);
            if f - f_new >= 1e-3 * alpha^2
                x = x - alpha * d;
                f = f_new;
            else
                alpha = alpha * theta;
            end
        end
    end
end