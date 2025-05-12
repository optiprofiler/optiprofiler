function custom_feature()

    % Define a custom feature that combines "perturbed_x0", "noisy", and "linearly_transformed".

    solvers = {@fminsearch, @fminunc};
    options.feature_name = 'custom';    % Set the feature name to "custom" if you want to use a custom feature.
    options.n_runs = 5;
    % We need mod_x0 to make sure that the linearly transformed problem is mathematically equivalent to the original problem.
    options.mod_x0 = @mod_x0;
    % We use mod_fun to add noise to the function value.
    options.mod_fun = @mod_fun;
    % We use mod_affine to define the linear transformation.
    options.mod_affine = @mod_affine;
    % We can define the feature stamp to represent the custom feature (optional).
    options.feature_stamp = 'perturbed_x0_noisy_linearly_transformed';

    scores = benchmark(solvers, options)
end

function x0 = mod_x0(rand_stream, problem)

    [Q, R] = qr(rand_stream.randn(problem.n));
    Q(:, diag(R) < 0) = -Q(:, diag(R) < 0);
    x0 = Q * problem.x0;
    x0 = x0 + 1e-3 * max(1, norm(x0)) * rand_stream.randn(problem.n, 1) / norm(rand_stream.randn(problem.n, 1));
end

function f = mod_fun(x, rand_stream, problem)

    f = problem.fun(x);
    f = f + max(1, abs(f)) * 1e-3 * rand_stream.randn(1);
end

function [A, b, inv] = mod_affine(rand_stream, problem)

    [Q, R] = qr(rand_stream.randn(problem.n));
    Q(:, diag(R) < 0) = -Q(:, diag(R) < 0);
    A = Q';
    b = zeros(problem.n, 1);
    inv = Q;
end