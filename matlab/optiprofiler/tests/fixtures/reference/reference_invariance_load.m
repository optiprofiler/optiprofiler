function problem = reference_invariance_load(name, varargin)
%REFERENCE_INVARIANCE_LOAD Problems of the reference-invariance fixture library.
%   The environment variable OPTIPROFILER_TEST_STATE_REFERENCES selects whether
%   the provider states reference facts ('1') or not (anything else). Nothing
%   else differs between the two variants: the library name, its root, the
%   problem names, the callbacks and therefore every seed are the same, so a
%   benchmark of one variant can be compared bit for bit with the other.
%
%   STEEP    minimize 1e7 * x subject to -x <= 0. The feasible optimum is 0 at
%            x = 0; an infeasible point has a negative objective.
%   SPHERE   minimize sum((x - 1).^2), unconstrained; -1 is a lower bound.
    with_reference = strcmp(getenv('OPTIPROFILER_TEST_STATE_REFERENCES'), '1');
    mapping = 'feasible_objective/1';
    switch name
        case 'STEEP'
            s = struct('name', name, 'fun', @(x) 1e7 * x(1), 'x0', 1, 'cub', @(x) -x(1));
            reference = struct('merit', 0, 'kind', 'optimum', 'source', 'analytic', 'mapping', mapping);
        case 'SPHERE'
            s = struct('name', name, 'fun', @(x) sum((x(:) - 1).^2), 'x0', [0; 2]);
            reference = struct('merit', -1, 'kind', 'lower_bound', 'source', 'analytic', 'mapping', mapping);
        otherwise
            error('ReferenceInvariance:UnknownProblem', 'Unknown problem %s.', name);
    end
    if with_reference
        s.reference = reference;
    end
    problem = Problem(s);
end
