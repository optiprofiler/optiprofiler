function merit_values = meritFunCompute(merit_fun, fun_values, maxcv_values, maxcv_init)
%MERITFUNCOMPUTE computes the merit function values.
%   `merit_fun` is a function handle given by the user.
%   We assume that `merit_fun` only accepts scalar inputs and returns scalar outputs.
%   This function is used to `vectorize` the merit function computation.
%   Note that it will not effect the results if `merit_fun` is already vectorized.

    % Scale `maxcv_init` to match the dimensions of `fun_values`.
    copied_dim = size(fun_values);
    if numel(maxcv_init) == 1
        maxcv_init = repmat(maxcv_init, copied_dim);
    else
        maxcv_init = repmat(maxcv_init, [1, copied_dim(2:end)]);
    end

    % Compute the merit function values by arrayfun.
    merit_values = arrayfun(@(f, cv, cv_init) merit_fun(f, cv, cv_init), fun_values, maxcv_values, maxcv_init, 'UniformOutput', false);

    % Convert the cell array to a numeric array.
    merit_values = cell2mat(merit_values);
end