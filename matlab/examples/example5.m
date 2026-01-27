function example5()
%EXAMPLE5 corresponds to "Example 5: customizing the test suite" in the "Usage for MATLAB" part of
% our official website (www.optprof.com).
%
% This example shows how to customize your own test suite.

    % Print the information about this example.
    fprintf('\nThis is an example to benchmark two toy solvers on two problem libraries (one is custom) with a custom feature.\n');
    pause(1.5);
    fprintf('\nStart Example 5...\n\n');

    % Start example 5.
    % You can custom your own problem library by following the steps in `optiprofiler/problem_libs/README.txt`.
    % After you finish the construction of your problem library, you can use it via the `plibs` option.
    % Here we use a custom problem library named "custom" as an example.
    options.plibs = {'custom', 's2mpj'};    % Here we use two problem libraries: "custom" and "s2mpj". You can use one or more problem libraries.
    % You can also select problems as in the previous examples.
    options.ptype = 'u';    % Select unconstrained optimization problems.
    options.mindim = 2;     % Select problems with dimension at least 2.

    % You can customize your own feature. More details can be found by typing `help benchmark` in the MATLAB command window.
    % In this example, we define a custom feature that combines "perturbed_x0", "noisy", and "linearly_transformed".
    options.feature_name = 'custom';    % Set the feature name to "custom" if you want to use a custom feature.
    options.mod_x0 = @mod_x0;           % Use mod_x0 to modify the initial guess.
    options.mod_fun = @mod_fun;         % Use mod_fun to add noise to the objective function value.
    options.mod_affine = @mod_affine;   % Use mod_affine to define the linear transformation.
    options.feature_stamp = 'perturbed_x0_noisy_linearly_transformed'; % Define the feature stamp to represent the custom feature (optional).
    options.n_runs = 3;                 % Set the number of runs for each solver on each problem.

    % `solver1` and `solver2` are two toy solvers for unconstrained optimization problems defined in
    % the `examples` folder only for demonstration.
    scores = benchmark({@solver1, @solver2}, options)
end

% Custom feature functions for the example5
function x0 = mod_x0(rand_stream, problem)

    [Q, R] = qr(rand_stream.randn(problem.n));
    Q(:, diag(R) < 0) = -Q(:, diag(R) < 0);
    x0 = Q * problem.x0;
    x0 = x0 + 1e-5 * max(1, norm(x0)) * rand_stream.randn(problem.n, 1) / norm(rand_stream.randn(problem.n, 1));
end

function f = mod_fun(x, rand_stream, problem)

    f = problem.fun(x);
    f = f + max(1, abs(f)) * 1e-6 * rand_stream.randn(1);
end

function [A, b, inv] = mod_affine(rand_stream, problem)

    [Q, R] = qr(rand_stream.randn(problem.n));
    Q(:, diag(R) < 0) = -Q(:, diag(R) < 0);
    A = Q';
    b = zeros(problem.n, 1);
    inv = Q;
end