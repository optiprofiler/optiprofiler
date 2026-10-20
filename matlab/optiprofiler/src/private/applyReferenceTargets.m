function targets = applyReferenceTargets(result, targets, merit_fun)
%APPLYREFERENCETARGETS Convert raw facts in the current merit space per run.
%   Missing records (V1 archives) or empty records keep the historical target.
%   Keep raw facts: load can change merit_fun, and maxcv_init varies by run.
    if ~isfield(result, 'reference_records')
        return
    end
    records = result.reference_records;
    if ~iscell(records) || ~isequal(size(records), size(targets))
        error('MATLAB:applyReferenceTargets:shape', 'Reference records must match the problem/run rows.');
    end
    for p = 1:size(targets, 1)
        for r = 1:size(targets, 2)
            if isempty(records{p, r})
                continue
            end
            reference = normalizeProblemReference(records{p, r});
            name = '';
            if isfield(result, 'problem_names') && numel(result.problem_names) >= p
                name = result.problem_names{p};
            end
            targets(p, r) = referenceTarget(reference, merit_fun, result.maxcv_inits(p, r), name, r);
        end
    end
end
