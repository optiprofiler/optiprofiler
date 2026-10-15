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
            value = merit_fun(reference.fun, reference.maxcv, result.maxcv_inits(p, r));
            if ~(isnumeric(value) && isreal(value) && isscalar(value) && isfinite(value))
                error('MATLAB:applyReferenceTargets:invalidMerit', ...
                    'The merit function must produce a finite real scalar for a supplied reference.');
            end
            targets(p, r) = double(value);
        end
    end
end
