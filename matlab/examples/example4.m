function example4()
%EXAMPLE4 corresponds to "Example 4: testing parametrized solvers" in the "Usage for MATLAB" part of
% our official website (www.optprof.com).
%
% This example shows how to benchmark with a parameterized solver.

    % Print the information about this example.
    fprintf('\nThis is an example to benchmark a parameterized toy solver with different parameters on the default problem set.\n');
    pause(1.5);
    fprintf('\nStart Example 4...\n\n');

    % Start example 4.
    % We will test the toy solver `solver` with parameters `1`, `2`, and `3`.
    solvers = cell(1, 3);
    solver_names = cell(1, 3);
    for i_solver = 1:3
        solvers{i_solver} = @(fun, x0) solver(fun, x0, i_solver);   % Create a function handle for the solver with the specific parameter.
        solver_names{i_solver} = sprintf('solver_%d', i_solver);    % Name the solver.
    end
    options.solver_names = solver_names;
    scores = benchmark(solvers, options)
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% A simple parameterized toy solver for demonstration.

function x = solver(fun, x0, para)

    % solver will randomly sample 50 * n points in the search space and return the one with the lowest
    % function value, where n is the dimension of the problem.
    % `para` is a parameter that controls the level of the random perturbation.
    n = length(x0);
    rand_stream = RandStream('mt19937ar', 'Seed', 1);
    xhist = NaN(n, 50 * n);
    fhist = NaN(50 * n, 1);
    for i = 1:50 * n
        xhist(:, i) = x0 + rand_stream.randn(n, 1) * para * 0.1;
        fhist(i) = fun(xhist(:, i));
    end
    [~, idx] = min(fhist);  % Find the index of the minimum function value.
    x = xhist(:, idx);  % Return the corresponding point.
end