function problem = composed_fixture_load(name)
% Problems with zero coordinates at the initial point (the payload case that
% collapsed the legacy product mixer) and a small shifted quadratic.
    switch name
        case 'zeros3'
            problem = Problem(struct('fun', @(x) sum(x.^2), 'x0', [0; 0; 0], 'name', name));
        case 'shift2'
            problem = Problem(struct('fun', @(x) sum((x - 1).^2), 'x0', [0; 2], 'name', name));
        otherwise
            error('Review:UnknownFixtureProblem', 'Unknown fixture problem %s.', name);
    end
end
