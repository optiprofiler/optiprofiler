function record = propagateProblemReference(record, name, options)
%PROPAGATEPROBLEMREFERENCE Reference fact of a problem after one feature stage.
%   RECORD = propagateProblemReference(RECORD, NAME, OPTIONS) returns the
%   validated reference RECORD of the stage's predecessor unchanged if the
%   stage NAME with its normalized stage OPTIONS retains it, and [] (unknown)
%   otherwise. An unknown input stays unknown.
%
%   The record is retained unchanged or dropped. It holds no point, so nothing
%   is transported, and nothing is derived from it. The rule is the same for
%   every kind (lower_bound, optimum, best_known, target) at zero violation.
%   A nonzero maxcv in best_known/target also requires the residual magnitude
%   to be invariant. General affine changes cannot guarantee this, so those
%   records become unknown even when feasibility is preserved.
%
%   - plain, noisy, truncated, random_nan, nonquantifiable_constraints,
%     unrelaxable_constraints: retained. They change only what a solver
%     observes; the truth recorded in the histories is the original problem's.
%   - perturbed_x0: retained. Only the initial point moves.
%   - permuted: retained. It preserves component residual magnitudes.
%   - linearly_transformed: retained only at maxcv=0. Bounds may become linear
%     rows with a different residual scale under a general affine map.
%   - quantized: retained with ground_truth=false, where the mesh is an
%     observation. With ground_truth=true (the default) the truth becomes the
%     mesh problem, whose feasible values are not those of the continuous
%     problem; no proof that the same feasible reference holds is available,
%     so the record becomes unknown.
%   - custom: retained only if the option names are a subset of
%     {mod_x0, mod_affine}. mod_x0 moves the initial point and mod_affine is a
%     change of variables subject to numerical matrix and transport checks
%     when the problem is built, not an exact inverse or uniform accuracy
%     certificate. The whitelist is deliberate: mod_fun, mod_cub, mod_ceq, mod_bounds,
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
    % Feasibility survives invertible maps, but a nonzero residual magnitude
    % may change when bounds become linear rows. No point is stored, so it
    % cannot safely be recomputed here.
    if retained && isfield(record, 'maxcv') && record.maxcv ~= 0
        if strcmp(char(name), 'linearly_transformed') ...
                || (strcmp(char(name), 'custom') && isfield(options, 'mod_affine'))
            retained = false;
        end
    end
    if ~retained
        record = [];
    end
end
