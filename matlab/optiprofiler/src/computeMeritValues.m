function merit_values = computeMeritValues(fun_values, maxcv_values, maxcv_init)
% The merit function varphi(x) is defined by the objective function f(x) and the maximum constraint violation v(x):
%   varphi(x) = f(x)                        if v(x) <= v1
%   varphi(x) = f(x) + 1e5 * (v(x) - v1)    if v1 < v(x) <= v2
%   varphi(x) = Inf                         if v(x) > v2
% where v1 = max(1e-5, v0) and v2 = min(0.01, 1e-10 * max(1, v0)), and v0 is the initial maximum constraint violation.
%
% Special case 1: computeMeritValues(NaN, maxcv_values, maxcv_init) = Inf
% Sepcial case 2: computeMeritValues(fun_values, NaN, maxcv_init) = Inf

    copied_dim = size(fun_values);
    merit_values = fun_values;
    maxcv_init = repmat(maxcv_init, [1, copied_dim(2:end)]);
    almost_feasible_upper_bound = max(1e-5, maxcv_init);   % Equal to 1e-5 if `maxcv_init` is NaN.
    almost_feasible_lower_bound = min(0.01, 1e-10 * max(1, maxcv_init));   % Equal to 1e-10 if `maxcv_init` is NaN.
    is_severely_infeasible = (maxcv_values > almost_feasible_upper_bound) | isnan(merit_values);    % Equal to True if `maxcv_values` is NaN.
    is_almost_feasible = (almost_feasible_lower_bound < maxcv_values) & (maxcv_values <= almost_feasible_upper_bound);   % Equal to False if `maxcv_values` is NaN.
    merit_values(is_severely_infeasible) = Inf;    % Covers the special case 1.
    merit_values(isnan(maxcv_values)) = Inf;    % Covers the special case 2.
    maxcv_values = max(maxcv_values - almost_feasible_lower_bound, 0);   % Shift the constraint violation.
    merit_values(is_almost_feasible) = merit_values(is_almost_feasible) + 1e5 * (maxcv_values(is_almost_feasible));
    
end