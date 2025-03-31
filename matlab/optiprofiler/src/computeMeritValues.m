function merit_values = computeMeritValues(fun_values, maxcv_values, maxcv_init)
    % Special case 1: computeMeritValues(NaN, maxcv_values, maxcv_init) = Inf
    % Sepcial case 2: computeMeritValues(fun_values, NaN, maxcv_init) = Inf

    copied_dim = size(fun_values);
    maxcv_init = repmat(maxcv_init, [1, copied_dim(2:end)]);
    infeasibility_thresholds = max(1e-5, maxcv_init);   % Equal to 1e-5 if `maxcv_init` is NaN.
    is_infeasible = maxcv_values > infeasibility_thresholds;    % Equal to False if `maxcv_values` is NaN.
    is_almost_feasible = (1e-10 < maxcv_values) & (maxcv_values <= infeasibility_thresholds);   % Equal to False if `maxcv_values` is NaN.
    merit_values = fun_values;
    merit_values(is_infeasible | isnan(merit_values)) = Inf;    % Covers the special case 1.
    merit_values(is_almost_feasible) = merit_values(is_almost_feasible) + 1e5 * maxcv_values(is_almost_feasible);
    merit_values(isnan(maxcv_values)) = Inf;    % Covers the special case 2.
end