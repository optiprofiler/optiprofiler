function other_options = getDefaultOtherOptions(solvers, other_options)

    if ~isfield(other_options, OtherOptionKey.SOLVER_NAMES.value)
        other_options.(OtherOptionKey.SOLVER_NAMES.value) = cellfun(@(s) func2str(s), solvers, 'UniformOutput', false);
    end
    if ~isfield(other_options, OtherOptionKey.SOLVER_ISRAND.value)
        other_options.(OtherOptionKey.SOLVER_ISRAND.value) = false(1, numel(solvers));
    end
    if isfield(other_options, OtherOptionKey.PROBLEM.value)
        other_options.(OtherOptionKey.CUTEST_PROBLEM_NAMES.value) = {};
        other_options.(OtherOptionKey.CUSTOM_PROBLEM_NAMES.value) = {other_options.(OtherOptionKey.PROBLEM.value).name};
        other_options.(OtherOptionKey.CUSTOM_PROBLEM_LOADER.value) = @(x) other_options.(OtherOptionKey.PROBLEM.value);
    end
    if ~isfield(other_options, OtherOptionKey.CUTEST_PROBLEM_NAMES.value)
        other_options.(OtherOptionKey.CUTEST_PROBLEM_NAMES.value) = {};
    end
    if ~isfield(other_options, OtherOptionKey.CUSTOM_PROBLEM_LOADER.value)
        other_options.(OtherOptionKey.CUSTOM_PROBLEM_LOADER.value) = {};
    end
    if ~isfield(other_options, OtherOptionKey.CUSTOM_PROBLEM_NAMES.value)
        other_options.(OtherOptionKey.CUSTOM_PROBLEM_NAMES.value) = {};
    end
end