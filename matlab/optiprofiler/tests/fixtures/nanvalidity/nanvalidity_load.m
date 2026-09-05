function problem = nanvalidity_load(name, varargin)
    problem = Problem(struct('fun', @(x) objective(name, x), 'x0', 0, ...
        'cub', @(x) constraint(name, x), 'name', name));
end

function value = objective(name, x)
    if (strcmp(name, 'NAN_INIT_FUN') && x(1) == 0) || ...
            (strcmp(name, 'NAN_LATER_FUN') && x(1) ~= 0)
        value = NaN;
    elseif strcmp(name, 'VALID_INF') || ...
            (strcmp(name, 'INF_INIT_NAN_LATER') && x(1) == 0)
        value = Inf;
    else
        value = 7;
    end
end

function value = constraint(name, x)
    if (strcmp(name, 'NAN_INIT_CV') && x(1) == 0) || ...
            (ismember(name, {'NAN_LATER_CV', 'INF_INIT_NAN_LATER'}) && x(1) ~= 0)
        value = NaN;
    else
        value = -1;
    end
end
