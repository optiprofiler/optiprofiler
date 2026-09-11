function options = refinedFeatureOptions(feature, plan, problem_options, profile_options, context)
%REFINEDFEATUREOPTIONS Versioned native settings for an explicit fresh replay.
% The exported effective stages retain native callback handles. Their JSON
% report descriptions are deliberately not accepted as executable substitutes.
    options = profile_options;
    names = fieldnames(problem_options);
    for i = 1:numel(names), options.(names{i}) = problem_options.(names{i}); end
    options.schema = 'options_refined-v2';
    % Preserve the canonical native state through Feature.saveobj/loadobj.
    % Reconstructing from its display spelling would lose unknown declaration
    % provenance and could reapply changed defaults to already-resolved stages.
    options.feature = feature;
    options.feature_route = context.route;
    options.feature_name = feature.name;
    stages = feature.stages;
    specification = cell(1, numel(stages));
    for i = 1:numel(stages)
        specification{i} = struct('name', stages{i}.name, 'options', stages{i}.options);
    end
    % Empty effective stages are valid identity but not a public empty input.
    if isempty(specification), specification = {struct('name', 'plain', 'options', struct())}; end
    % This duplicate is inspection-only. The replay loader rejects disagreement
    % instead of silently selecting one of two user-edited representations.
    options.feature_specification = specification;
    options.n_runs = plan.n_runs;
    options.feature_provenance = optiprofiler_internal.featureProvenance(feature, plan, context);
end
