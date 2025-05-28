function example1()
%EXAMPLE1 corresponds to "Example 1: first example to try out" in the "Usage for MATLAB" part of our
% official website (www.optprof.com).
%
% This example shows the simplest way to use OptiProfiler.

    % Call OptiProfiler via `benchmark` function.
    % `solver1` and `solver2` are two toy solvers for unconstrained optimization problems defined in
    % the `examples` folder only for demonstration.
    % They can be replaced with function handles accepting a function and an initial point,
    % e.g., `@fminunc`, `@fminsearch`, etc.
    % More possible signatures can be found by typing `help benchmark` in the MATLAB command window.

    scores = benchmark({@solver1, @solver2})
end