function merit_values = meritFunCompute(merit_fun, fun_values, maxcv_values, maxcv_inits)
%MERITFUNCOMPUTE computes the merit function values.
%   `merit_fun` is a function handle given by the user.
%   We assume that `merit_fun` only accepts scalar inputs and returns scalar outputs.
%   This function is used to `vectorize` the merit function computation.
%   Note that it will not effect the results if `merit_fun` is already vectorized.
%
%   Note that `fun_values` and `maxcv_values` should have the same size.
%   There are four possible sizes of `fun_values` and `maxcv_inits`:
%   1. fun_values has size (n_problems, n_solvers, n_runs, n_evals), and maxcv_inits has size (n_problems, n_runs)
%   2. fun_values has size (n_problems, n_solvers, n_runs), and maxcv_inits has size (n_problems, n_runs)
%   3. fun_values has size (n_solvers, n_runs, n_evals), and maxcv_inits has size (n_runs, 1)
%   4. fun_values has size (n_runs, 1), and maxcv_inits has size (n_runs, 1)
%

    % Scale `maxcv_init` to match the dimensions of `fun_values`.
    size_fun = size(fun_values);
    size_inits = size(maxcv_inits);
    if numel(size_fun) == 4 && isequal(size_inits, [size_fun(1), size_fun(3)])
        maxcv_inits = reshape(maxcv_inits, [size_fun(1), 1, size_fun(3), 1]);
        maxcv_inits = repmat(maxcv_inits, [1, size_fun(2), 1, size_fun(4)]);
    elseif numel(size_fun) == 3 && isequal(size_inits, [size_fun(1), size_fun(3)])
        maxcv_inits = reshape(maxcv_inits, [size_fun(1), 1, size_fun(3)]);
        maxcv_inits = repmat(maxcv_inits, [1, size_fun(2), 1]);
    elseif numel(size_fun) == 3 && (isequal(size_inits, [size_fun(2), 1]) || isequal(size_inits, [1, size_fun(2)]))
        maxcv_inits = reshape(maxcv_inits, [1, size_fun(2), 1]);
        maxcv_inits = repmat(maxcv_inits, [size_fun(1), 1, size_fun(3)]);
    elseif numel(size_fun) == 2 && isequal(size_inits, size_fun)
        % Already the target size, do nothing.
    elseif isequal(size_inits, size_fun)
        % Already the target size, do nothing.
    else
        error("MATLAB:meritFunCompute:size_mismatch", "The size of maxcv_inits does not match the size of fun_values.");
    end

    % Compute the merit function values by arrayfun.
    merit_values = arrayfun(@(f, cv, cv_init) merit_fun(f, cv, cv_init), fun_values, maxcv_values, maxcv_inits, 'UniformOutput', false);

    % Convert the cell array to a numeric array.
    merit_values = cell2mat(merit_values);
end