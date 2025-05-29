function merit_value = defaultMerit(fun_value, maxcv_value, maxcv_init)
% The merit function varphi(x) is defined by the objective function f(x) and the maximum constraint violation v(x):
%   varphi(x) = f(x)                        if v(x) <= v1
%   varphi(x) = f(x) + 1e5 * (v(x) - v1)    if v1 < v(x) <= v2
%   varphi(x) = Inf                         if v(x) > v2
% where v1 = max(1e-5, v0) and v2 = min(0.01, 1e-10 * max(1, v0)), and v0 is the initial maximum constraint violation.
%
% Special case 1: defaultMerit(NaN, maxcv_value, maxcv_init) = Inf
% Sepcial case 2: defaultMerit(fun_value, NaN, maxcv_init) = Inf

    if isnan(fun_value) || isnan(maxcv_value)
        merit_value = Inf;   % Special case 1 and 2.
        return;
    end

    merit_value = fun_value;   % Initialize the merit value with the objective function value.
    v1 = max(1e-5, maxcv_init);   % Equal to 1e-5 if `maxcv_init` is NaN.
    v2 = min(0.01, 1e-10 * max(1, maxcv_init));   % Equal to 1e-10 if `maxcv_init` is NaN.
    
    if maxcv_value > v2
        merit_value = Inf;   % Case 3: v(x) > v2.
    elseif maxcv_value > v1
        merit_value = fun_value + 1e5 * (maxcv_value - v1);   % Case 2: v1 < v(x) <= v2.
    end
end