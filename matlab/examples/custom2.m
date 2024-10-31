function custom2()
    % Define a custom feature that combines "noisy" and "linearly_transformed".
    % We only modify mod_fun since we are dealing with unconstrained problems.
    
    solvers = {@fminsearch_test, @fminunc_test};
    options.feature_name = 'custom';
    options.n_runs = 3;
    options.mod_fun = @mod_fun;
    options.mod_affine = @mod_affine;

    benchmark(solvers, options)
end

function f = mod_fun(x, f, rand_stream, problem)

    f = f + max(1, abs(f)) * 1e-3 * rand_stream.randn(1);
end

function [A, b, inv] = mod_affine(rand_stream, problem)

    [Q, R] = qr(rand_stream.randn(problem.n));
    Q(:, diag(R) < 0) = -Q(:, diag(R) < 0);
    A = Q';
    b = zeros(problem.n, 1);
    inv = Q;
end

function x = fminsearch_test(fun, x0)

    x = fminsearch(fun, x0);
    
end

function x = fminunc_test(fun, x0)

    x = fminunc(fun, x0);

end