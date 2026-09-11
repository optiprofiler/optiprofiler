function plan = resolveFeatureExperiment(feature, profile_options, role)
%RESOLVEFEATUREEXPERIMENT Resolve one role without mutating a feature spec.
% Inputs are already validated by benchmark. In particular, an explicitly empty
% n_runs is invalid and must not arrive here as a synonym for omission. The
% literal hints are not equivalent to is_stochastic (custom: 1; affine: 5).
% planned_runs records intended calls, never observed execution/success facts.
    if nargin < 3, role = 'primary'; end
    stages = feature.stages;
    if isempty(stages)
        strategy = 'identity';
    elseif numel(stages) == 1
        strategy = 'legacy-single';
    else
        strategy = 'composed-views';
    end
    % The generic strategy and its language-specific implementation version
    % are distinct provenance facts. Match the actual FeaturedProblem policy
    % without constructing a runtime or changing any legacy numerical stream.
    runtime_policy = 'matlab-legacy-single-v1';
    if strcmp(strategy, 'composed-views'), runtime_policy = 'matlab-composed-views-v1'; end
    requested = isfield(profile_options, 'n_runs');
    count = [];
    if requested, count = profile_options.n_runs; end
    solver_isrand = profile_options.solver_isrand(:);
    if strcmp(role, 'plain_reference')
        % This role is intentionally independent of the primary request, even
        % for randomized solvers. Sharing its count would change the baseline.
        requested = false;
        count = [];
        n_runs = 1;
        origin = 'reference_policy';
    elseif ~strcmp(role, 'primary')
        error('OptiProfiler:UnknownExperimentRole', 'Unknown experiment role: %s.', role);
    elseif requested
        n_runs = count;
        origin = 'explicit';
    elseif any(solver_isrand)
        n_runs = 5;
        origin = 'randomized_solvers';
    else
        n_runs = 1;
        for i = 1:numel(stages)
            definition = optiprofiler_internal.featureDefinitions(stages{i}.name, stages{i}.options);
            n_runs = max(n_runs, definition.replicate_hint);
        end
        origin = 'stage_hints';
    end
    planned_runs = ones(numel(solver_isrand), 1);
    if strcmp(role, 'primary')
        planned_runs(feature.is_stochastic | solver_isrand) = n_runs;
    end
    plan = struct('role', role, 'n_runs', n_runs, 'origin', origin, ...
        'run_policy', 'legacy-hints-v1', 'execution_strategy', strategy, ...
        'runtime_policy', runtime_policy, ...
        'request_present', requested, 'requested_n_runs', count, ...
        'planned_runs', planned_runs);
end
