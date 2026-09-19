function record = propagateProblemReference(record, name, options)
%PROPAGATEPROBLEMREFERENCE Reference fact of a problem after one feature stage.
%   RECORD = propagateProblemReference(RECORD, NAME, OPTIONS) returns the
%   validated reference RECORD of the stage's predecessor unchanged if the
%   stage NAME with its normalized stage OPTIONS retains it, and [] (unknown)
%   otherwise. An unknown input stays unknown.
%
%   The record is retained unchanged or dropped. It holds no point, so nothing
%   is transported, and nothing is derived from it. The rule is the same for
%   every kind (lower_bound, optimum, best_known, target): each is a claim
%   over the feasible points of the problem, so it survives exactly the stages
%   that leave the feasible points and the objective values over them
%   unchanged up to a change of variables.
%
%   - plain, noisy, truncated, random_nan, nonquantifiable_constraints,
%     unrelaxable_constraints: retained. They change only what a solver
%     observes; the truth recorded in the histories is the original problem's.
%   - perturbed_x0: retained. Only the initial point moves.
%   - permuted, linearly_transformed: retained. They are invertible changes of
%     variables built by the framework.
%   - quantized: retained with ground_truth=false, where the mesh is an
%     observation. With ground_truth=true (the default) the truth becomes the
%     mesh problem, whose feasible values are not those of the continuous
%     problem; no proof that the same feasible reference holds is available,
%     so the record becomes unknown.
%   - custom: retained only if the option names are a subset of
%     {mod_x0, mod_affine}. mod_x0 moves the initial point and mod_affine is a
%     change of variables whose inverse is checked when the problem is built
%     (A * inv == I), so it is a valid affine coordinate change. This is a
%     whitelist on purpose: mod_fun, mod_cub, mod_ceq, mod_bounds,
%     mod_linear_ub, mod_linear_eq and any option added later are arbitrary
%     user code that may change values, constraints or bounds, and the
%     framework cannot prove that the feasible reference survives them.
%   - any other name: unknown (fail closed).
%
%   A composition applies this rule stage by stage. Once a stage reports
%   unknown no later stage can restore the record, so a composition is safe
%   exactly when every stage is. Only the stage name and its option names are
%   read; no callback is called, so applying the rule evaluates nothing.

    if isempty(record)
        record = [];
        return
    end
    switch char(name)
        case {'plain', 'noisy', 'truncated', 'random_nan', 'nonquantifiable_constraints', ...
                'unrelaxable_constraints', 'perturbed_x0', 'permuted', 'linearly_transformed'}
            retained = true;
        case 'quantized'
            % Retained only for an explicit logical false; a missing or
            % non-logical value counts as ground truth.
            retained = isstruct(options) && isscalar(options) && isfield(options, 'ground_truth') ...
                && islogical(options.ground_truth) && isscalar(options.ground_truth) && ~options.ground_truth;
        case 'custom'
            retained = isstruct(options) && isscalar(options) ...
                && all(ismember(fieldnames(options), {'mod_x0', 'mod_affine'}));
        otherwise
            retained = false;
    end
    if ~retained
        record = [];
    end
end
