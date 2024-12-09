function custom1()

    clc
    
    % Define a custom feature that combines "perturbed_x0" and
    % "linearly_transformed".
    
    solvers = {@fminsearch_test1, @fminsearch_test2};
    options.feature_name = 'custom';
    options.n_runs = 5;
    options.mod_x0 = @mod_x0;
    options.mod_affine = @mod_affine;

    benchmark(solvers, options)
end

function x0 = mod_x0(rand_stream, problem)

    [Q, R] = qr(rand_stream.randn(problem.n));
    Q(:, diag(R) < 0) = -Q(:, diag(R) < 0);
    x0 = Q * problem.x0;
    x0 = x0 + 1e-3 * max(1, norm(x0)) * rand_stream.randn(problem.n, 1) / norm(rand_stream.randn(problem.n, 1));
end

function [A, b, inv] = mod_affine(rand_stream, problem)

    [Q, R] = qr(rand_stream.randn(problem.n));
    Q(:, diag(R) < 0) = -Q(:, diag(R) < 0);
    A = Q';
    b = zeros(problem.n, 1);
    inv = Q;
end

function x = fminsearch_test1(fun, x0)

    options = optimset('MaxFunEvals', 200);
    x = fminsearch(fun, x0, options);
end

function x = fminsearch_test2(fun, x0)

    options = optimset('MaxFunEvals', 500);
    x = fminsearch(fun, x0, options);
end