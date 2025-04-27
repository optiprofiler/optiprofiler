function cutest_options = getDefaultCutestOptions(cutest_options, other_options)
    
    % If the user does not provide the cutest_options (empty struct) and provides cutest_problem_names or custom_problem_names with custom_problem_loader, we will not set the default options.
    if numel(fieldnames(cutest_options)) == 0 && (~isempty(other_options.(OtherOptionKey.CUTEST_PROBLEM_NAMES.value)) || (~isempty(other_options.(OtherOptionKey.CUSTOM_PROBLEM_NAMES.value)) && ~isempty(other_options.(OtherOptionKey.CUSTOM_PROBLEM_LOADER.value))))
        return;
    end

    if ~isfield(cutest_options, CutestOptionKey.EXCLUDELIST.value)
        cutest_options.(CutestOptionKey.EXCLUDELIST.value) = {};
    end
    
    if ~isfield(cutest_options, CutestOptionKey.PTYPE.value)
        cutest_options.(CutestOptionKey.PTYPE.value) = 'u';
    end
    if ~isfield(cutest_options, CutestOptionKey.MINDIM.value)
        cutest_options.(CutestOptionKey.MINDIM.value) = 1;
    end
    if ~isfield(cutest_options, CutestOptionKey.MAXDIM.value)
        cutest_options.(CutestOptionKey.MAXDIM.value) = cutest_options.(CutestOptionKey.MINDIM.value) + 1;
    end
    if ~isfield(cutest_options, CutestOptionKey.MINB.value)
        cutest_options.(CutestOptionKey.MINB.value) = 0;
    end
    if ~isfield(cutest_options, CutestOptionKey.MAXB.value)
        cutest_options.(CutestOptionKey.MAXB.value) = cutest_options.(CutestOptionKey.MINB.value) + 10;
    end
    if ~isfield(cutest_options, CutestOptionKey.MINCON.value)
        cutest_options.(CutestOptionKey.MINCON.value) = 0;
    end
    if ~isfield(cutest_options, CutestOptionKey.MAXCON.value)
        cutest_options.(CutestOptionKey.MAXCON.value) = cutest_options.(CutestOptionKey.MINCON.value) + 10;
    end
end