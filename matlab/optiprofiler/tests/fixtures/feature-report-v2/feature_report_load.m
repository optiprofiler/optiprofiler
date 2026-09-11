function problem = feature_report_load(name)
% Hand-checkable full-library fixture, with no external provider installation.
    if strcmp(getenv('OP_FEATURE_D_FORBID_PROVIDER'), '1')
        error('OptiProfiler:ForbiddenProvider', 'Loading saved results must not execute a provider.');
    end
    if strcmp(name, 'SMALL'), n = 2;
    elseif strcmp(name, 'WIDE'), n = 3;
    else, error('OptiProfiler:UnknownFixtureProblem', 'Unknown fixture problem.');
    end
    problem = Problem(struct('name', name, 'fun', @(x) sum(x.^2), 'x0', (1:n)'));
end
