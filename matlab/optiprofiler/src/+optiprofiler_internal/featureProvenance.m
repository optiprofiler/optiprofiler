function payload = featureProvenance(feature, plan, context)
%FEATUREPROVENANCE Describe explicit specification and one role's plan.
% This is native metadata, not an executable replay specification. The caller
% uses EvalReport.encodeMetadata for JSON; callback objects are never invoked.
% Invocation route is distinct from the object's original construction route.
    if nargin < 3, context = struct(); end
    unknown = struct('eval_report_null', true);
    declaration = feature.declared;
    declared = unknown;
    declaration_route = unknown;
    declared_name = unknown;
    if ~isempty(declaration)
        declared = declaration.entries;
        declaration_route = declaration.route;
        declared_name = feature.declared_name;
    end
    stages = feature.stages;
    records = cell(1, numel(stages));
    for i = 1:numel(stages)
        stage = stages{i};
        % Shared feature-pipeline positions and same-kind occurrences are
        % zero-based, independent of MATLAB's one-based numerical run indexes.
        records{i} = struct('position', i - 1, 'name', stage.name, ...
            'code', stage.code, 'occurrence', stage.occurrence, ...
            'identity', stage.identity, 'options', stage.options);
    end
    if numel(stages) > 1
        seed_policy = 'matlab-stage-horner32-v1';
    else
        seed_policy = 'legacy-run-seed';
    end
    feature_record = struct('route', optional(context, 'route', unknown), ...
        'declared_name', declared_name, 'declared', {declared}, ...
        'declaration_route', declaration_route, 'effective_name', feature.name, ...
        'seed_policy', seed_policy, 'feature_stamp', optional(context, 'feature_stamp', unknown), ...
        'full_feature_stamp', optional(context, 'full_feature_stamp', unknown), ...
        'feature_stamp_origin', optional(context, 'feature_stamp_origin', unknown), ...
        'stages', {records});
    experiment = unknown;
    if ~isempty(plan)
        % These facts are an intention/policy, not observed solver calls.
        experiment = struct('role', plan.role, 'n_runs', plan.n_runs, ...
            'origin', plan.origin, 'run_policy', plan.run_policy, ...
            'execution_strategy', plan.execution_strategy);
    end
    payload = struct('schema', 'feature_pipeline-v3', ...
        'feature', feature_record, 'experiment', experiment);
end

function value = optional(source, name, fallback)
    value = fallback;
    if isfield(source, name), value = source.(name); end
end
