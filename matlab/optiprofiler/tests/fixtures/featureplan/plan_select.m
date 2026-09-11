function names = plan_select(options)
% Two real unconstrained quadratics, selected through the public registry.
    names = {};
    if ~contains(options.ptype, 'u'), return; end
    for key = {'minb', 'minlcon', 'minnlcon', 'mincon'}
        if options.(key{1}) > 0, return; end
    end
    for dimension = [2, 3]
        if dimension >= options.mindim && dimension <= options.maxdim
            names{end + 1} = sprintf('PLAN%d', dimension); %#ok<AGROW>
        end
    end
end
