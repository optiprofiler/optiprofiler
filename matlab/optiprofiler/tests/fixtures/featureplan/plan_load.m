function problem = plan_load(name)
% Real Problem objects; no replacement Feature or solver-driver shim.
    switch name
        case 'PLAN2', dimension = 2;
        case 'PLAN3', dimension = 3;
        otherwise, error('FeaturePlanFixture:UnknownProblem', 'Unknown problem: %s.', name);
    end
    problem = Problem(struct('fun', @(x) sum(x.^2), ...
        'x0', (1:dimension)', 'name', name));
end
