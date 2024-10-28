function cutest_problem_names = loadCutestNames(cutest_options, cutest_problem_names, custom_problem_loader, custom_problem_names)
%LOADCUTEST Load the CUTEst problems according to the options in cutest_options

    % Select the problems based on cutest_options.
    if isempty(cutest_options) || ~isempty(cutest_problem_names) || (~isempty(custom_problem_loader) || ~isempty(custom_problem_names))
        cutest_problem_names_options = {};
    else
        cutest_problem_names_options = selector(cutest_options);
    end

    % Preprocess the CUTEst problem names given by the user.
    if ~isempty(cutest_problem_names)
        if ~ischarstr(cutest_problem_names) && ~(iscell(cutest_problem_names) && all(cellfun(@ischarstr, cutest_problem_names)))
            error("MATLAB:benchmark:cutest_problem_namesNotValid", "The CUTEst problem names must be a charstr or a cell array of charstr.");
        end
        if ischarstr(cutest_problem_names)
            cutest_problem_names = {cutest_problem_names};
        end
        % Convert to a cell row vector of chars.
        cutest_problem_names = cellfun(@char, cutest_problem_names, 'UniformOutput', false);
        cutest_problem_names = cutest_problem_names(:)';
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

    % Exclude the problems specified by cutest_options.excludelist.
    if isfield(cutest_options, CutestOptionKey.EXCLUDELIST.value) && ~isempty(cutest_options.(CutestOptionKey.EXCLUDELIST.value))
        cutest_problem_names = setdiff(cutest_problem_names, cutest_options.(CutestOptionKey.EXCLUDELIST.value));
    end

    % Check whether the total number of problems is zero.
    if isempty(cutest_problem_names) && isempty(custom_problem_names)
        error("MATLAB:benchmark:AtLeastOneProblem", "At least one problem must be given.");
    end

end