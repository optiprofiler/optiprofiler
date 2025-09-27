function merit_values = meritFunCompute(merit_fun, fun_values, maxcv_values, maxcv_inits, type)
%MERITFUNCOMPUTE computes the merit function values.
%   `merit_fun` is a function handle given by the user.
%   We assume that `merit_fun` only accepts scalar inputs and returns scalar outputs.
%   This function is used to `vectorize` the merit function computation.
%   Note that it will not effect the results if `merit_fun` is already vectorized.
%
%   Note that `fun_values` and `maxcv_values` should have the same size.
%   There are following possible sizes of `fun_values` and `maxcv_inits`:
%
%       fun_values has the same size as maxcv_inits
%
%       fun_values has size (n_problems, n_solvers, n_runs, n_evals),                and maxcv_inits has size (n_problems, n_runs)      case 1
%                        or (n_problems, n_solvers, n_runs) (if n_evals = 1),        and maxcv_inits has size (n_problems, n_runs)      case 2
%                        or (n_problems, n_solvers) (if n_runs = 1 and n_evals = 1), and maxcv_inits has size (n_problems, 1)           case 3
%                        or (n_solvers, n_runs, n_evals),                            and maxcv_inits has size (n_runs, 1)               case 4
%                        or (n_solvers, n_runs) (if n_evals = 1),                    and maxcv_inits has size (n_runs, 1)               case 5
%
%   To distinguish the above cases, we can use the argument `type`:
%       type = 'multiple' for case 1, 2, or 3
%       type = 'single' for case 4 and 5
%   We do not need to distinguish the case where fun_values has the same size as maxcv_inits, since we will first check if their sizes are the same.

    % If `fun_values`, `maxcv_values`, or `maxcv_inits` is a row vector, convert it to a column vector.
    if isrealrow(fun_values)
        fun_values = fun_values(:);
    end
    if isrealrow(maxcv_values)
        maxcv_values = maxcv_values(:);
    end
    if isrealrow(maxcv_inits)
        maxcv_inits = maxcv_inits(:);
    end

    % Scale `maxcv_init` to match the dimensions of `fun_values`.
    size_fun = size(fun_values);
    size_inits = size(maxcv_inits);

    if isequal(size_inits, size_fun)
        % Already the target size, do nothing.
    elseif nargin < 5
        error("MATLAB:meritFunCompute:size_mismatch", "The size of maxcv_inits cannot match the size of fun_values.");
    elseif strcmp(type, 'multiple')
        if numel(size_fun) == 4 && isequal(size_inits, [size_fun(1), size_fun(3)])  % Case 1
            maxcv_inits = reshape(maxcv_inits, [size_fun(1), 1, size_fun(3), 1]);
            maxcv_inits = repmat(maxcv_inits, [1, size_fun(2), 1, size_fun(4)]);
        elseif numel(size_fun) == 3 && isequal(size_inits, [size_fun(1), size_fun(3)])  % Case 2
            maxcv_inits = reshape(maxcv_inits, [size_fun(1), 1, size_fun(3)]);
            maxcv_inits = repmat(maxcv_inits, [1, size_fun(2), 1]);
        elseif numel(size_fun) == 2 && isequal(size_inits, [size_fun(1), 1])  % Case 3
            maxcv_inits = reshape(maxcv_inits, [size_fun(1), 1]);
            maxcv_inits = repmat(maxcv_inits, [1, size_fun(2)]);
        else
            error("MATLAB:meritFunCompute:size_mismatch", "The size of maxcv_inits cannot match the size of fun_values.");
        end
    elseif strcmp(type, 'single')
        if numel(size_fun) == 3 && isequal(size_inits, [size_fun(2), 1])  % Case 4
            maxcv_inits = reshape(maxcv_inits, [1, size_fun(2), 1]);
            maxcv_inits = repmat(maxcv_inits, [size_fun(1), 1, size_fun(3)]);
        elseif numel(size_fun) == 2 && isequal(size_inits, [size_fun(2), 1])  % Case 5
            maxcv_inits = reshape(maxcv_inits, [1, size_fun(2)]);
            maxcv_inits = repmat(maxcv_inits, [size_fun(1), 1]);
        else
            error("MATLAB:meritFunCompute:size_mismatch", "The size of maxcv_inits cannot match the size of fun_values.");
        end
    else
        error("MATLAB:meritFunCompute:size_mismatch", "The size of maxcv_inits cannot match the size of fun_values.");
    end

    % Compute the merit function values by arrayfun.
    merit_values = arrayfun(@(f, cv, cv_init) merit_fun(f, cv, cv_init), fun_values, maxcv_values, maxcv_inits, 'UniformOutput', false);

    % Convert the cell array to a numeric array.
    merit_values = cell2mat(merit_values);
end