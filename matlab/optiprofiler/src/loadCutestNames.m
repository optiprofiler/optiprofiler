function cutest_problem_names = loadCutestNames(cutest_options, other_options)
%LOADCUTEST Load the CUTEst problems according to the options in cutest_options

    cutest_problem_names = other_options.(OtherOptionKey.CUTEST_PROBLEM_NAMES.value);
    custom_problem_names = other_options.(OtherOptionKey.CUSTOM_PROBLEM_NAMES.value);
    
    % Select the problems based on cutest_options.
    if numel(fieldnames(cutest_options)) == 0
        cutest_problem_names_options = {};
    else
        cutest_problem_names_options = s_select(cutest_options);
    end

    % Merge the problem names selected by cutest_options and given by the user.
    % N.B.: Duplicate names are and MUST BE removed.
    if isempty(cutest_problem_names) && isempty(custom_problem_names)
        cutest_problem_names = cutest_problem_names_options;
    elseif isempty(cutest_problem_names) && ~isempty(custom_problem_names)
        % Do nothing.
    elseif numel(fieldnames(cutest_options)) == 0
        cutest_problem_names = unique(cutest_problem_names);
    else
        cutest_problem_names = unique([cutest_problem_names, cutest_problem_names_options]);
    end

    % Add some problems to excludelist. These problems have bugs or are not suitable for benchmarking.
    % 'DANWOODLS' and 'MISRA1CLS' have bugs at this moment (bugs in SIF files where bound constraints are not defined).
    % 'ROSSIMP1', 'ROSSIMP2', and 'ROSSIMP3' are duplicated with 'ROSENBR'.
    if isfield(cutest_options, CutestOptionKey.EXCLUDELIST.value)
        cutest_options.(CutestOptionKey.EXCLUDELIST.value) = unique([cutest_options.(CutestOptionKey.EXCLUDELIST.value), {'DANWOODLS', 'MISRA1CLS', 'ROSSIMP1', 'ROSSIMP2', 'ROSSIMP3'}]);
    else
        cutest_options.(CutestOptionKey.EXCLUDELIST.value) = {'DANWOODLS', 'MISRA1CLS', 'ROSSIMP1', 'ROSSIMP2', 'ROSSIMP3'};
    end

    % Exclude the problems specified by cutest_options.excludelist.
    if isfield(cutest_options, CutestOptionKey.EXCLUDELIST.value) && ~isempty(cutest_options.(CutestOptionKey.EXCLUDELIST.value))
        % Convert to a cell row vector of chars.
        cutest_options.(CutestOptionKey.EXCLUDELIST.value) = cellfun(@char, cutest_options.(CutestOptionKey.EXCLUDELIST.value), 'UniformOutput', false);
        cutest_problem_names = setdiff(cutest_problem_names, cutest_options.(CutestOptionKey.EXCLUDELIST.value));
    end

    % Check whether the total number of problems is zero.
    if isempty(cutest_problem_names) && isempty(custom_problem_names)
        error("MATLAB:benchmark:AtLeastOneProblem", "At least one problem must be given.");
    end

end