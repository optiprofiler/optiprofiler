function other_options = checkValidityOtherOptions(solvers, other_options)
%CHECKVALIDITYOTHEROPTIONS Check the validity of the options in other_options

    % Judge whether other_options.solver_names is a cell of chars or strings.
    if isfield(other_options, OtherOptionKey.SOLVER_NAMES.value)
        if ~iscell(other_options.(OtherOptionKey.SOLVER_NAMES.value)) || ~all(cellfun(@(l) ischarstr(l), other_options.(OtherOptionKey.SOLVER_NAMES.value)))
            error("MATLAB:checkValidityOtherOptions:solver_namesNotCellOfcharstr", "The field `solver_names` of `options` for `benchmark` must be a cell of chars or strings.");
        end
        if numel(other_options.(OtherOptionKey.SOLVER_NAMES.value)) ~= 0 && numel(other_options.(OtherOptionKey.SOLVER_NAMES.value)) ~= numel(solvers)
            error("MATLAB:checkValidityOtherOptions:solver_namesAndsolversLengthNotSame", "The number of the field `solver_names` of `options` for `benchmark` must equal the number of solvers.");
        end
        if numel(other_options.(OtherOptionKey.SOLVER_NAMES.value)) == 0
            other_options.(OtherOptionKey.SOLVER_NAMES.value) = cellfun(@(s) func2str(s), solvers, 'UniformOutput', false);
        end
    end
    % Judge whether other_options.solver_isrand is a logical array of the same length as the number of solvers.
    if isfield(other_options, OtherOptionKey.SOLVER_ISRAND.value)
        if ~islogical(other_options.(OtherOptionKey.SOLVER_ISRAND.value)) || numel(other_options.(OtherOptionKey.SOLVER_ISRAND.value)) ~= numel(solvers)
            error("MATLAB:checkValidityOtherOptions:solver_israndNotLogical", "The field `solver_isrand` of `options` for `benchmark` must be a logical array of the same length as the number of solvers.");
        end
    end
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
    if isfield(other_options, OtherOptionKey.CUSTOM_PROBLEM_LOADER.value)
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
    if isfield(other_options, OtherOptionKey.CUSTOM_PROBLEM_LOADER.value) && isfield(other_options, OtherOptionKey.CUSTOM_PROBLEM_NAMES.value)
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