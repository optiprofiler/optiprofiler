function problem_options = checkValidityProblemOptions(problem_options)
%CHECKVALIDITYCUTESTOPTIONS Check the validity of the options in cuteset_options

    % Judge whether problem_options.plibs is a char or a string or a cell array of strings or chars.
    % Judge whether all the members of problem_options.plibs belong to the subfolder names under the relative path '../../../problems/' to this file.
    if isfield(problem_options, ProblemOptionKey.PLIBS.value)
        if ~iscell(problem_options.(ProblemOptionKey.PLIBS.value))
            problem_options.(ProblemOptionKey.PLIBS.value) = {problem_options.(ProblemOptionKey.PLIBS.value)};
        end
        mydir = fileparts(mfilename('fullpath'));
        problems_dir = fullfile(mydir, '../../../../problems');
        subfolder_info = dir(problems_dir);
        subfolder_names = {subfolder_info([subfolder_info.isdir] & ~ismember({subfolder_info.name}, {'.', '..'})).name};
        if isempty(problem_options.(ProblemOptionKey.PLIBS.value)) || ~all(cellfun(@ischarstr, problem_options.(ProblemOptionKey.PLIBS.value))) || ~all(ismember(problem_options.(ProblemOptionKey.PLIBS.value), subfolder_names))
            error("MATLAB:checkValidityProblemOptions:plibsNotValid", "The field 'plibs' of options should be a nonempty cell array of strings or chars selected from '%s'.", strjoin(subfolder_names, "', '"));
        end
        % Convert to cell array of chars.
        problem_options.(ProblemOptionKey.PLIBS.value) = cellfun(@char, problem_options.(ProblemOptionKey.PLIBS.value), 'UniformOutput', false);
        % Convert to a row vector if not.
        if size(problem_options.(ProblemOptionKey.PLIBS.value), 1) > 1
            problem_options.(ProblemOptionKey.PLIBS.value) = problem_options.(ProblemOptionKey.PLIBS.value)';
        end
    end
    % If the OS is not Linux, then the field 'plibs' of options cannot contain 'matcutest'.
    if (~isunix || ismac) && isfield(problem_options, ProblemOptionKey.PLIBS.value) && any(ismember(problem_options.(ProblemOptionKey.PLIBS.value), 'matcutest'))
        error("MATLAB:checkValidityProblemOptions:plibsNotLinux", "The field 'plibs' of options cannot contain 'matcutest' on non-Linux OS.");
    end


    % Judge whether problem_options.ptype is among all possible problem types.
    if isfield(problem_options, ProblemOptionKey.PTYPE.value)
        if ~ischarstr(problem_options.(ProblemOptionKey.PTYPE.value))
            error("MATLAB:checkValidityProblemOptions:ptypeNotcharstr", "The field 'ptype' of options should be a char or a string.");
        else
            % Convert to lower case CHAR.
            problem_options.(ProblemOptionKey.PTYPE.value) = lower(char(problem_options.(ProblemOptionKey.PTYPE.value)));
            % Check whether the problem type belongs to four types ('u', 'b', 'l', 'n') and their combinations.
            if ~all(ismember(problem_options.(ProblemOptionKey.PTYPE.value), 'ubln'))
                error("MATLAB:checkValidityProblemOptions:ptypeNotubln", "The field 'ptype' of options should be a string containing only 'u', 'b', 'l', 'n'.");
            end
        end
    end


    % Judge whether problem_options.mindim is a integer greater or equal to 1.
    if isfield(problem_options, ProblemOptionKey.MINDIM.value)
        if ~isintegerscalar(problem_options.(ProblemOptionKey.MINDIM.value)) || problem_options.(ProblemOptionKey.MINDIM.value) < 1
            error("MATLAB:checkValidityProblemOptions:mindimNotValid", "The field of 'mindim' of options should be a integer greater or equal to 1.");
        end
    end
    % Judge whether problem_options.maxdim is a integer greater or equal to 1, or equal to Inf.
    if isfield(problem_options, ProblemOptionKey.MAXDIM.value)
        if (~isintegerscalar(problem_options.(ProblemOptionKey.MAXDIM.value)) || problem_options.(ProblemOptionKey.MAXDIM.value) < 1) && problem_options.(ProblemOptionKey.MAXDIM.value) ~= Inf
            error("MATLAB:checkValidityProblemOptions:maxdimNotValid", "The field 'maxdim' of options should be a integer greater or equal to 1, or equal to Inf.");
        end
    end
    % Judge whether problem_options.mindim is smaller or equal to problem_options.maxdim.
    if isfield(problem_options, ProblemOptionKey.MINDIM.value) && isfield(problem_options, ProblemOptionKey.MAXDIM.value)
        if problem_options.(ProblemOptionKey.MINDIM.value) > problem_options.(ProblemOptionKey.MAXDIM.value)
            error("MATLAB:checkValidityProblemOptions:maxdimSmallerThanmindim", "The field 'mindim' of options should be smaller or equal to problem_options.maxdim.");
        end
    end


    % Judge whether problem_options.minb is a integer greater or equal to 0.
    if isfield(problem_options, ProblemOptionKey.MINB.value)
        if ~isintegerscalar(problem_options.(ProblemOptionKey.MINB.value)) || problem_options.(ProblemOptionKey.MINB.value) < 0
            error("MATLAB:checkValidityProblemOptions:minbNotValid", "The field 'minb' of options should be a integer greater or equal to 0.");
        end
    end
    % Judge whether problem_options.maxb is a integer greater or equal to 0, or equal to Inf.
    if isfield(problem_options, ProblemOptionKey.MAXB.value)
        if (~isintegerscalar(problem_options.(ProblemOptionKey.MAXB.value)) || problem_options.(ProblemOptionKey.MAXB.value) < 0) && problem_options.(ProblemOptionKey.MAXB.value) ~= Inf
            error("MATLAB:checkValidityProblemOptions:maxbNotValid", "The field 'maxb' of options should be a integer greater or equal to 0, or equal to Inf.");
        end
    end
    % Judge whether problem_options.minb is smaller or equal to problem_options.maxb.
    if isfield(problem_options, ProblemOptionKey.MINB.value) && isfield(problem_options, ProblemOptionKey.MAXB.value)
        if problem_options.(ProblemOptionKey.MINB.value) > problem_options.(ProblemOptionKey.MAXB.value)
            error("MATLAB:checkValidityProblemOptions:maxbSmallerThanminb", "The field 'minb' of options should be smaller or equal to problem_options.maxb.");
        end
    end


    % Judge whether problem_options.mincon is a integer greater or equal to 0.
    if isfield(problem_options, ProblemOptionKey.MINCON.value)
        if ~isintegerscalar(problem_options.(ProblemOptionKey.MINCON.value)) || problem_options.(ProblemOptionKey.MINCON.value) < 0
            error("MATLAB:checkValidityProblemOptions:minconNotValid", "The field 'mincon' of options should be a integer greater or equal to 0.");
        end
    end
    % Judge whether problem_options.maxcon is a integer greater or equal to 0, or equal to Inf.
    if isfield(problem_options, ProblemOptionKey.MAXCON.value)
        if (~isintegerscalar(problem_options.(ProblemOptionKey.MAXCON.value)) || problem_options.(ProblemOptionKey.MAXCON.value) < 0) && problem_options.(ProblemOptionKey.MAXCON.value) ~= Inf
            error("MATLAB:checkValidityProblemOptions:maxconNotValid", "The field 'maxcon' of options should be a integer greater or equal to 0, or equal to Inf.");
        end
    end
    % Judge whether problem_options.mincon is smaller or equal to problem_options.maxcon.
    if isfield(problem_options, ProblemOptionKey.MINCON.value) && isfield(problem_options, ProblemOptionKey.MAXCON.value)
        if problem_options.(ProblemOptionKey.MINCON.value) > problem_options.(ProblemOptionKey.MAXCON.value)
            error("MATLAB:checkValidityProblemOptions:maxconSmallerThanmincon", "The field 'mincon' of options should be smaller or equal to problem_options.maxcon.");
        end
    end


    % Judge whether problem_options.excludelist is a cell array of strings or chars.
    if isfield(problem_options, ProblemOptionKey.EXCLUDELIST.value)
        if ischarstr(problem_options.(ProblemOptionKey.EXCLUDELIST.value))
            problem_options.(ProblemOptionKey.EXCLUDELIST.value) = {problem_options.(ProblemOptionKey.EXCLUDELIST.value)};
        end
        if ~iscell(problem_options.(ProblemOptionKey.EXCLUDELIST.value)) || ~all(cellfun(@ischarstr, problem_options.(ProblemOptionKey.EXCLUDELIST.value)))
            error("MATLAB:checkValidityProblemOptions:excludelistNotCellOfcharstr", "The field 'excludelist' of options should be a cell array of strings or chars.");
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
            error("MATLAB:checkValidityProblemOptions:problem_namesNotCellOfcharstr", "The field 'problem_names' of options should be a cell array of strings or chars.");
        end
        % Convert to cell array of chars.
        problem_options.(ProblemOptionKey.PROBLEM_NAMES.value) = cellfun(@char, problem_options.(ProblemOptionKey.PROBLEM_NAMES.value), 'UniformOutput', false);
        % Convert to a row vector if not.
        if size(problem_options.(ProblemOptionKey.PROBLEM_NAMES.value), 1) > 1
            problem_options.(ProblemOptionKey.PROBLEM_NAMES.value) = problem_options.(ProblemOptionKey.PROBLEM_NAMES.value)';
        end
    end
end