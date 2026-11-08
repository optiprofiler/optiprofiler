function integral = integrate(curve, profile_type, profile_options)
    % Compute the integral of the curve from different types of profiles.

    if isempty(curve)
        integral = 0;
        return;
    end

    kernel = profile_options.(ProfileOptionKey.SCORE_WEIGHT_FUN.value);
    integral = 0;
    switch profile_type
        case {'perf', 'data'}
            % The curve is a right-continuous step function.
            contribution = diff(curve(1, :)) .* curve(2, 1:end-1) .* kernel(curve(1, 1:end-1));
            integral = integral + sum(contribution(isfinite(contribution)));
        case 'log_ratio'
            % We do not modify the integral of log_ratio even a score_weight_fun is provided.
            contribution = abs(curve(2, :));
            integral = integral + sum(contribution(isfinite(contribution)));
    end
end
