function names = feature_report_select(options)
% Small real provider boundary: one whole-library selection, never one-by-one.
    names = {'SMALL', 'WIDE'};
    dims = [2, 3];
    keep = true(size(dims));
    if isfield(options, 'mindim'), keep = keep & dims >= options.mindim; end
    if isfield(options, 'maxdim'), keep = keep & dims <= options.maxdim; end
    if isfield(options, 'ptype'), keep = keep & contains(options.ptype, 'u'); end
    names = names(keep);
end
