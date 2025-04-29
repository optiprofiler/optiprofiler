function other_options = checkValidityOtherOptions(solvers, other_options)
%CHECKVALIDITYOTHEROPTIONS Check the validity of the options in other_options

    % Judge whether other_options.problem is a Problem object.
    if isfield(other_options, OtherOptionKey.PROBLEM.value)
        if ~isa(other_options.(OtherOptionKey.PROBLEM.value), 'Problem')
            error("MATLAB:checkValidityOtherOptions:problemNotProblem", "The field `problem` of `options` for `benchmark` must be a Problem object.");
        end
    end
    % Judge whether other_options.cutest_problem_names is a cell of chars or strings.
    if isfield(other_options, OtherOptionKey.CUTEST_PROBLEM_NAMES.value)
        if ~isempty(other_options.(OtherOptionKey.CUTEST_PROBLEM_NAMES.value))
            if ~ischarstr(other_options.(OtherOptionKey.CUTEST_PROBLEM_NAMES.value)) && ~(iscell(other_options.(OtherOptionKey.CUTEST_PROBLEM_NAMES.value)) && all(cellfun(@ischarstr, other_options.(OtherOptionKey.CUTEST_PROBLEM_NAMES.value))))
                error("MATLAB:checkValidityOtherOptions:cutest_problem_namesNotValid", "The CUTEst problem names must be a charstr or a cell array of charstr.");
            end
            if ischarstr(other_options.(OtherOptionKey.CUTEST_PROBLEM_NAMES.value))
                other_options.(OtherOptionKey.CUTEST_PROBLEM_NAMES.value) = {other_options.(OtherOptionKey.CUTEST_PROBLEM_NAMES.value)};
            end
            % Convert to a cell row vector of chars.
            other_options.(OtherOptionKey.CUTEST_PROBLEM_NAMES.value) = cellfun(@char, other_options.(OtherOptionKey.CUTEST_PROBLEM_NAMES.value), 'UniformOutput', false);
            other_options.(OtherOptionKey.CUTEST_PROBLEM_NAMES.value) = other_options.(OtherOptionKey.CUTEST_PROBLEM_NAMES.value)(:)';
        end
    end
    % Judge whether other_options.custom_problem_loader is a function handle.
    if isfield(other_options, OtherOptionKey.CUSTOM_PROBLEM_LOADER.value) && ~isempty(other_options.(OtherOptionKey.CUSTOM_PROBLEM_LOADER.value))
        if ~isa(other_options.(OtherOptionKey.CUSTOM_PROBLEM_LOADER.value), 'function_handle')
            error("MATLAB:checkValidityOtherOptions:customloaderNotFunctionHandle", "The field `custom_problem_loader` of `options` for `benchmark` must be a function handle.");
        end
    end
    % Judge whether other_options.custom_problem_names is a cell of chars or strings.
    if isfield(other_options, OtherOptionKey.CUSTOM_PROBLEM_NAMES.value)
        if ~ischarstr(other_options.(OtherOptionKey.CUSTOM_PROBLEM_NAMES.value)) && ~(iscell(other_options.(OtherOptionKey.CUSTOM_PROBLEM_NAMES.value)) && all(cellfun(@ischarstr, other_options.(OtherOptionKey.CUSTOM_PROBLEM_NAMES.value))))
            error("MATLAB:checkValidityOtherOptions:customnamesNotcharstrOrCellOfcharstr", "The field `custom_problem_names` of `options` for `benchmark` must be a cell array of chars or strings.");
        end
        if ischarstr(other_options.(OtherOptionKey.CUSTOM_PROBLEM_NAMES.value))
            other_options.(OtherOptionKey.CUSTOM_PROBLEM_NAMES.value) = {other_options.(OtherOptionKey.CUSTOM_PROBLEM_NAMES.value)};
        end
    end
    % Judge whether the custom_problem_loader can load the first custom problem.
    if isfield(other_options, OtherOptionKey.CUSTOM_PROBLEM_LOADER.value) && ~isempty(other_options.(OtherOptionKey.CUSTOM_PROBLEM_LOADER.value)) && isfield(other_options, OtherOptionKey.CUSTOM_PROBLEM_NAMES.value) && ~isempty(other_options.(OtherOptionKey.CUSTOM_PROBLEM_NAMES.value))
        custom_problem_loader = other_options.(OtherOptionKey.CUSTOM_PROBLEM_LOADER.value);
        custom_problem_names = other_options.(OtherOptionKey.CUSTOM_PROBLEM_NAMES.value);
        try
            [~, p] = evalc('custom_problem_loader(custom_problem_names{1})');
        catch
            p = [];
        end
        if isempty(p) || ~isa(p, 'Problem')
            error("MATLAB:checkValidityOtherOptions:customloaderNotAcceptcustomnames", "The `custom_problem_loader` failed to load the first custom problem %s.", custom_problem_names{1});
        end
    end

end