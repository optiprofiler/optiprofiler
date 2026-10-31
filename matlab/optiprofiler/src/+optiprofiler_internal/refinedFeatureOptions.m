function options = refinedFeatureOptions(feature, plan, problem_options, profile_options, context, direct_problem)
%REFINEDFEATUREOPTIONS Versioned native settings for an explicit fresh replay.
% The exported effective stages retain native callback handles. Their JSON
% report descriptions are deliberately not accepted as executable substitutes.
% A directly supplied Problem is not stored: problem_route records that a
% replay needs it, and the unused library selectors are omitted.
    if nargin < 6, direct_problem = false; end
    options = profile_options;
    if direct_problem
        options.problem_route = 'direct_problem';
    else
        names = fieldnames(problem_options);
        for i = 1:numel(names), options.(names{i}) = problem_options.(names{i}); end
        options.problem_route = 'library_selection';
    end
    options.schema = 'options_refined-v2';
    % Preserve the canonical native state through Feature.saveobj/loadobj.
    % Reconstructing from its display spelling would lose unknown declaration
    % provenance and could reapply changed defaults to already-resolved stages.
    options.feature = feature;
    options.feature_route = [];
    if isfield(context, 'route'), options.feature_route = context.route; end
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
