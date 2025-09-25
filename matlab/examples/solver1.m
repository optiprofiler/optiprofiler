function x = solver1(fun, x0)
%SOLVER1 is a TOY SOLVER for demonstration. It will be used in the examples under the `examples` folder.
%
% Note: This is a toy solver and does not guarantee convergence or optimality.
% It is intended for demonstration purposes only.
% In practice, you would replace this with a more sophisticated optimization algorithm.

    % solver1 will randomly sample 50 * n points in the search space and return the one with the lowest
    % function value, where n is the dimension of the problem.
    n = length(x0);
    n_eval_factor = 50;
    n_eval = n_eval_factor * n;
    rand_stream = RandStream('mt19937ar', 'Seed', 1);
    xhist = NaN(n, n_eval);
    fhist = NaN(n_eval, 1);
    for i = 1:n_eval
        xhist(:, i) = x0 + rand_stream.randn(n, 1);
        fhist(i) = fun(xhist(:, i));
    end
    [~, idx] = min(fhist);  % Find the index of the minimum function value.
    x = xhist(:, idx);  % Return the corresponding point.
end