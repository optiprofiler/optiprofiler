function example2()
%EXAMPLE2 corresponds to "Example 2: one step further by adding options" in the "Usage for MATLAB"
% part of our official website (www.optprof.com).
%
% This example shows how to use OptiProfiler with options. More options can be found by typing
% `help benchmark` in the MATLAB command window.

    options.ptype = 'u';    % Select unconstrained optimization problems. You can change it to any combination of 'u' (unconstrained), 'b' (box-constrained), 'l' (linearly constrained), and 'n' (nonlinearly constrained).
    options.mindim = 2;     % Select problems with dimension at least 2.
    options.maxdim = 5;     % Select problems with dimension at most 5.
    options.feature_name = 'noisy';    % Choose the 'noisy' feature.
    options.n_runs = 3;     % Repeat each solver 3 (default is 5) times on each problem.
    options.max_eval_factor = 10;    % Set the maximum number of function evaluations to 10 (default is 500) times the dimension of the problem.

    % `solver1` and `solver2` are two toy solvers for unconstrained optimization problems defined in
    % the `examples` folder only for demonstration.
    scores = benchmark({@solver1, @solver2, @fminsearch}, options)
end