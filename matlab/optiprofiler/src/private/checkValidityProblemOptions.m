function problem_options = checkValidityProblemOptions(problem_options, profile_options)
%CHECKVALIDITYCUTESTOPTIONS Check the validity of the options in cuteset_options

    % Validate problem_options.plibs through the explicit MATLAB registry.
    if isfield(problem_options, ProblemOptionKey.PLIBS.value)
        if ~iscell(problem_options.(ProblemOptionKey.PLIBS.value))
            problem_options.(ProblemOptionKey.PLIBS.value) = {problem_options.(ProblemOptionKey.PLIBS.value)};
        end
        valid_name_types = cellfun( ...
            @(name) ischarstr(name) && isscalar(string(name)), ...
            problem_options.(ProblemOptionKey.PLIBS.value));
        if isempty(problem_options.(ProblemOptionKey.PLIBS.value)) || ...
                ~all(valid_name_types)
            error("MATLAB:checkValidityProblemOptions:plibsNotValid", ...
                "The option `plibs` should be a nonempty cell array of strings or chars.");
        end
        % Convert to cell array of chars.
        problem_options.(ProblemOptionKey.PLIBS.value) = cellfun( ...
            @(name) lower(char(name)), ...
            problem_options.(ProblemOptionKey.PLIBS.value), ...
            'UniformOutput', false);
        % Convert to a row vector if not.
        if size(problem_options.(ProblemOptionKey.PLIBS.value), 1) > 1
            problem_options.(ProblemOptionKey.PLIBS.value) = problem_options.(ProblemOptionKey.PLIBS.value)';
        end

        is_loading = isfield(profile_options, ProfileOptionKey.LOAD.value) && ...
            ~isempty(profile_options.(ProfileOptionKey.LOAD.value));
        % MatCUTEst can label previously saved results on every OS, but live
        % problem resolution remains Linux-only.
        if (~isunix || ismac) && ...
                any(ismember(problem_options.(ProblemOptionKey.PLIBS.value), 'matcutest')) && ...
                ~is_loading
            error("MATLAB:checkValidityProblemOptions:plibsNotLinux", ...
                "The option `plibs` cannot contain 'matcutest' on non-Linux OS.");
        end

        for i_plib = 1:numel(problem_options.(ProblemOptionKey.PLIBS.value))
            plib = problem_options.(ProblemOptionKey.PLIBS.value){i_plib};
            try
                if is_loading
                    if ~isKnownProblemLibrary(plib)
                        error("MATLAB:resolveProblemLibrary:notRegistered", ...
                            "Problem library '%s' is not registered.", plib);
                    end
                else
                    resolveProblemLibrary(plib);
                end
            catch ME
                error("MATLAB:checkValidityProblemOptions:plibsNotValid", ...
                    "Problem library '%s' is unavailable or invalid: %s", ...
                    plib, ME.message);
            end
        end
    end


    % Judge whether problem_options.ptype is among all possible problem types.
    if isfield(problem_options, ProblemOptionKey.PTYPE.value)
        if ~ischarstr(problem_options.(ProblemOptionKey.PTYPE.value))
            error("MATLAB:checkValidityProblemOptions:ptypeNotcharstr", "The option `ptype` should be a char or a string.");
        else
            % Convert to lower case CHAR.
            problem_options.(ProblemOptionKey.PTYPE.value) = lower(char(problem_options.(ProblemOptionKey.PTYPE.value)));
            % Check whether the problem type belongs to four types ('u', 'b', 'l', 'n') and their combinations.
            if ~all(ismember(problem_options.(ProblemOptionKey.PTYPE.value), 'ubln'))
                error("MATLAB:checkValidityProblemOptions:ptypeNotubln", "The option `ptype` should be a string containing only 'u', 'b', 'l', 'n'.");
            end
        end
    end


    % Judge whether problem_options.mindim is a integer greater or equal to 1.
    if isfield(problem_options, ProblemOptionKey.MINDIM.value)
        if ~isintegerscalar(problem_options.(ProblemOptionKey.MINDIM.value)) || problem_options.(ProblemOptionKey.MINDIM.value) < 1
            error("MATLAB:checkValidityProblemOptions:mindimNotValid", "The option of `mindim` should be a integer greater or equal to 1.");
        end
    end
    % Judge whether problem_options.maxdim is a integer greater or equal to 1, or equal to Inf.
    if isfield(problem_options, ProblemOptionKey.MAXDIM.value)
        if (~isintegerscalar(problem_options.(ProblemOptionKey.MAXDIM.value)) || problem_options.(ProblemOptionKey.MAXDIM.value) < 1) && problem_options.(ProblemOptionKey.MAXDIM.value) ~= Inf
            error("MATLAB:checkValidityProblemOptions:maxdimNotValid", "The option `maxdim` should be a integer greater or equal to 1, or equal to Inf.");
        end
    end
    % Judge whether problem_options.mindim is smaller or equal to problem_options.maxdim.
    if isfield(problem_options, ProblemOptionKey.MINDIM.value) && isfield(problem_options, ProblemOptionKey.MAXDIM.value)
        if problem_options.(ProblemOptionKey.MINDIM.value) > problem_options.(ProblemOptionKey.MAXDIM.value)
            error("MATLAB:checkValidityProblemOptions:maxdimSmallerThanmindim", "The option `mindim` should be smaller or equal to `maxdim`.");
        end
    end


    % Judge whether problem_options.minb is a integer greater or equal to 0.
    if isfield(problem_options, ProblemOptionKey.MINB.value)
        if ~isintegerscalar(problem_options.(ProblemOptionKey.MINB.value)) || problem_options.(ProblemOptionKey.MINB.value) < 0
            error("MATLAB:checkValidityProblemOptions:minbNotValid", "The option `minb` should be a integer greater or equal to 0.");
        end
    end
    % Judge whether problem_options.maxb is a integer greater or equal to 0, or equal to Inf.
    if isfield(problem_options, ProblemOptionKey.MAXB.value)
        if (~isintegerscalar(problem_options.(ProblemOptionKey.MAXB.value)) || problem_options.(ProblemOptionKey.MAXB.value) < 0) && problem_options.(ProblemOptionKey.MAXB.value) ~= Inf
            error("MATLAB:checkValidityProblemOptions:maxbNotValid", "The option `maxb` should be a integer greater or equal to 0, or equal to Inf.");
        end
    end
    % Judge whether problem_options.minb is smaller or equal to problem_options.maxb.
    if isfield(problem_options, ProblemOptionKey.MINB.value) && isfield(problem_options, ProblemOptionKey.MAXB.value)
        if problem_options.(ProblemOptionKey.MINB.value) > problem_options.(ProblemOptionKey.MAXB.value)
            error("MATLAB:checkValidityProblemOptions:maxbSmallerThanminb", "The option `minb` should be smaller or equal to `maxb`.");
        end
    end


    % Judge whether problem_options.minlcon is a integer greater or equal to 0.
    if isfield(problem_options, ProblemOptionKey.MINLCON.value)
        if ~isintegerscalar(problem_options.(ProblemOptionKey.MINLCON.value)) || problem_options.(ProblemOptionKey.MINLCON.value) < 0
            error("MATLAB:checkValidityProblemOptions:minlconNotValid", "The option `minlcon` should be a integer greater or equal to 0.");
        end
    end
    % Judge whether problem_options.maxlcon is a integer greater or equal to 0, or equal to Inf.
    if isfield(problem_options, ProblemOptionKey.MAXLCON.value)
        if (~isintegerscalar(problem_options.(ProblemOptionKey.MAXLCON.value)) || problem_options.(ProblemOptionKey.MAXLCON.value) < 0) && problem_options.(ProblemOptionKey.MAXLCON.value) ~= Inf
            error("MATLAB:checkValidityProblemOptions:maxlconNotValid", "The option `maxlcon` should be a integer greater or equal to 0, or equal to Inf.");
        end
    end
    % Judge whether problem_options.minlcon is smaller or equal to problem_options.maxlcon.
    if isfield(problem_options, ProblemOptionKey.MINLCON.value) && isfield(problem_options, ProblemOptionKey.MAXLCON.value)
        if problem_options.(ProblemOptionKey.MINLCON.value) > problem_options.(ProblemOptionKey.MAXLCON.value)
            error("MATLAB:checkValidityProblemOptions:maxlconSmallerThanminlcon", "The option `minlcon` should be smaller or equal to `maxlcon`.");
        end
    end


    % Judge whether problem_options.minnlcon is a integer greater or equal to 0.
    if isfield(problem_options, ProblemOptionKey.MINNLCON.value)
        if ~isintegerscalar(problem_options.(ProblemOptionKey.MINNLCON.value)) || problem_options.(ProblemOptionKey.MINNLCON.value) < 0
            error("MATLAB:checkValidityProblemOptions:minnlconNotValid", "The option `minnlcon` should be a integer greater or equal to 0.");
        end
    end
    % Judge whether problem_options.maxnlcon is a integer greater or equal to 0, or equal to Inf.
    if isfield(problem_options, ProblemOptionKey.MAXNLCON.value)
        if (~isintegerscalar(problem_options.(ProblemOptionKey.MAXNLCON.value)) || problem_options.(ProblemOptionKey.MAXNLCON.value) < 0) && problem_options.(ProblemOptionKey.MAXNLCON.value) ~= Inf
            error("MATLAB:checkValidityProblemOptions:maxnlconNotValid", "The option `maxnlcon` should be a integer greater or equal to 0, or equal to Inf.");
        end
    end
    % Judge whether problem_options.minnlcon is smaller or equal to problem_options.maxnlcon.
    if isfield(problem_options, ProblemOptionKey.MINNLCON.value) && isfield(problem_options, ProblemOptionKey.MAXNLCON.value)
        if problem_options.(ProblemOptionKey.MINNLCON.value) > problem_options.(ProblemOptionKey.MAXNLCON.value)
            error("MATLAB:checkValidityProblemOptions:maxnlconSmallerThanminnlcon", "The option `minnlcon` should be smaller or equal to `maxnlcon`.");
        end
    end


    % Judge whether problem_options.mincon is a integer greater or equal to 0.
    if isfield(problem_options, ProblemOptionKey.MINCON.value)
        if ~isintegerscalar(problem_options.(ProblemOptionKey.MINCON.value)) || problem_options.(ProblemOptionKey.MINCON.value) < 0
            error("MATLAB:checkValidityProblemOptions:minconNotValid", "The option `mincon` should be a integer greater or equal to 0.");
        end
    end
    % Judge whether problem_options.maxcon is a integer greater or equal to 0, or equal to Inf.
    if isfield(problem_options, ProblemOptionKey.MAXCON.value)
        if (~isintegerscalar(problem_options.(ProblemOptionKey.MAXCON.value)) || problem_options.(ProblemOptionKey.MAXCON.value) < 0) && problem_options.(ProblemOptionKey.MAXCON.value) ~= Inf
            error("MATLAB:checkValidityProblemOptions:maxconNotValid", "The option `maxcon` should be a integer greater or equal to 0, or equal to Inf.");
        end
    end
    % Judge whether problem_options.mincon is smaller or equal to problem_options.maxcon.
    if isfield(problem_options, ProblemOptionKey.MINCON.value) && isfield(problem_options, ProblemOptionKey.MAXCON.value)
        if problem_options.(ProblemOptionKey.MINCON.value) > problem_options.(ProblemOptionKey.MAXCON.value)
            error("MATLAB:checkValidityProblemOptions:maxconSmallerThanmincon", "The option `mincon` should be smaller or equal to `maxcon`.");
        end
    end


    % Judge whether problem_options.excludelist is a cell array of strings or chars.
    if isfield(problem_options, ProblemOptionKey.EXCLUDELIST.value)
        if ischarstr(problem_options.(ProblemOptionKey.EXCLUDELIST.value))
            problem_options.(ProblemOptionKey.EXCLUDELIST.value) = {problem_options.(ProblemOptionKey.EXCLUDELIST.value)};
        end
        if ~iscell(problem_options.(ProblemOptionKey.EXCLUDELIST.value)) || ~all(cellfun(@ischarstr, problem_options.(ProblemOptionKey.EXCLUDELIST.value)))
            error("MATLAB:checkValidityProblemOptions:excludelistNotCellOfcharstr", "The option `excludelist` should be a cell array of strings or chars.");
        end
        % Convert to cell array of chars.
        problem_options.(ProblemOptionKey.EXCLUDELIST.value) = cellfun(@char, problem_options.(ProblemOptionKey.EXCLUDELIST.value), 'UniformOutput', false);
        % Convert to a row vector if not.
        if size(problem_options.(ProblemOptionKey.EXCLUDELIST.value), 1) > 1
            problem_options.(ProblemOptionKey.EXCLUDELIST.value) = problem_options.(ProblemOptionKey.EXCLUDELIST.value)';
        end
    end


    % Judge whether problem_options.problem_names is a cell array of strings or chars.
    if isfield(problem_options, ProblemOptionKey.PROBLEM_NAMES.value)
        if ischarstr(problem_options.(ProblemOptionKey.PROBLEM_NAMES.value))
            problem_options.(ProblemOptionKey.PROBLEM_NAMES.value) = {problem_options.(ProblemOptionKey.PROBLEM_NAMES.value)};
        end
        if ~iscell(problem_options.(ProblemOptionKey.PROBLEM_NAMES.value)) || ~all(cellfun(@ischarstr, problem_options.(ProblemOptionKey.PROBLEM_NAMES.value)))
            error("MATLAB:checkValidityProblemOptions:problem_namesNotCellOfcharstr", "The option `problem_names` should be a cell array of strings or chars.");
        end
        % Convert to cell array of chars.
        problem_options.(ProblemOptionKey.PROBLEM_NAMES.value) = cellfun(@char, problem_options.(ProblemOptionKey.PROBLEM_NAMES.value), 'UniformOutput', false);
        % Convert to a row vector if not.
        if size(problem_options.(ProblemOptionKey.PROBLEM_NAMES.value), 1) > 1
            problem_options.(ProblemOptionKey.PROBLEM_NAMES.value) = problem_options.(ProblemOptionKey.PROBLEM_NAMES.value)';
        end
    end
end
