function merit_values = computeMeritValues(fun_values, maxcv_values, maxcv_init)

    copied_dim = size(fun_values);
    maxcv_init = repmat(maxcv_init, [1, copied_dim(2:end)]);
    infeasibility_thresholds = max(1e-5, maxcv_init);
    is_infeasible = maxcv_values > infeasibility_thresholds;
    is_almost_feasible = (1e-10 < maxcv_values) & (maxcv_values <= infeasibility_thresholds);
    merit_values = fun_values;
    merit_values(is_infeasible | isnan(merit_values)) = Inf;
    merit_values(is_almost_feasible) = merit_values(is_almost_feasible) + 1e5 * maxcv_values(is_almost_feasible);
end