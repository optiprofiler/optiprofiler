function cutest_options = checkValidityCutestOptions(cutest_options)
%CHECKVALIDITYCUTESTOPTIONS Check the validity of the options in cuteset_options

    % Judge whether cutest_options.ptype is among all possible problem types.
    if isfield(cutest_options, CutestOptionKey.PTYPE.value)
        if ~ischarstr(cutest_options.(CutestOptionKey.PTYPE.value))
            error("MATLAB:checkValidityCutestOptions:ptypeNotcharstr", "The field 'ptype' of options should be a char or a string.");
        else
            % Convert to lower case CHAR.
            cutest_options.(CutestOptionKey.PTYPE.value) = lower(char(cutest_options.(CutestOptionKey.PTYPE.value)));
            % Check whether the problem type belongs to four types ('u', 'b', 'l', 'n') and their combinations.
            if ~all(ismember(cutest_options.(CutestOptionKey.PTYPE.value), 'ubln'))
                error("MATLAB:checkValidityCutestOptions:ptypeNotubln", "The field 'ptype' of options should be a string containing only 'u', 'b', 'l', 'n'.");
            end
        end
    end
    % Judge whether cutest_options.mindim is a integer greater or equal to 1.
    if isfield(cutest_options, CutestOptionKey.MINDIM.value)
        if ~isintegerscalar(cutest_options.(CutestOptionKey.MINDIM.value)) || cutest_options.(CutestOptionKey.MINDIM.value) < 1
            error("MATLAB:checkValidityCutestOptions:mindimNotValid", "The field of 'mindim' of options should be a integer greater or equal to 1.");
        end
    end
    % Judge whether cutest_options.maxdim is a integer greater or equal to 1, or equal to Inf.
    if isfield(cutest_options, CutestOptionKey.MAXDIM.value)
        if ~isintegerscalar(cutest_options.(CutestOptionKey.MAXDIM.value)) || (cutest_options.(CutestOptionKey.MAXDIM.value) < 1 && cutest_options.(CutestOptionKey.MAXDIM.value) ~= Inf)
            error("MATLAB:checkValidityCutestOptions:maxdimNotValid", "The field 'maxdim' of options should be a integer greater or equal to 1, or equal to Inf.");
        end
    end
    % Judge whether cutest_options.mindim is smaller or equal to cutest_options.maxdim.
    if isfield(cutest_options, CutestOptionKey.MINDIM.value) && isfield(cutest_options, CutestOptionKey.MAXDIM.value)
        if cutest_options.(CutestOptionKey.MINDIM.value) > cutest_options.(CutestOptionKey.MAXDIM.value)
            error("MATLAB:checkValidityCutestOptions:maxdimSmallerThanmindim", "The field 'mindim' of options should be smaller or equal to cutest_options.maxdim.");
        end
    end
    % Judge whether cutest_options.mincon is a integer greater or equal to 0.
    if isfield(cutest_options, CutestOptionKey.MINCON.value)
        if ~isintegerscalar(cutest_options.(CutestOptionKey.MINCON.value)) || cutest_options.(CutestOptionKey.MINCON.value) < 0
            error("MATLAB:checkValidityCutestOptions:minconNotValid", "The field 'mincon' of options should be a integer greater or equal to 0.");
        end
    end
    % Judge whether cutest_options.maxcon is a integer greater or equal to 0, or equal to Inf.
    if isfield(cutest_options, CutestOptionKey.MAXCON.value)
        if ~isintegerscalar(cutest_options.(CutestOptionKey.MAXCON.value)) || (cutest_options.(CutestOptionKey.MAXCON.value) < 0 && cutest_options.(CutestOptionKey.MAXCON.value) ~= Inf)
            error("MATLAB:checkValidityCutestOptions:maxconNotValid", "The field 'maxcon' of options should be a integer greater or equal to 0, or equal to Inf.");
        end
    end
    % Judge whether cutest_options.mincon is smaller or equal to cutest_options.maxcon.
    if isfield(cutest_options, CutestOptionKey.MINCON.value) && isfield(cutest_options, CutestOptionKey.MAXCON.value)
        if cutest_options.(CutestOptionKey.MINCON.value) > cutest_options.(CutestOptionKey.MAXCON.value)
            error("MATLAB:checkValidityCutestOptions:maxconSmallerThanmincon", "The field 'mincon' of options should be smaller or equal to cutest_options.maxcon.");
        end
    end
    % Judge whether cutest_options.excludelist is a cell array of strings or chars.
    if isfield(cutest_options, CutestOptionKey.EXCLUDELIST.value)
        if ischarstr(cutest_options.(CutestOptionKey.EXCLUDELIST.value))
            cutest_options.(CutestOptionKey.EXCLUDELIST.value) = {cutest_options.(CutestOptionKey.EXCLUDELIST.value)};
        end
        if ~iscell(cutest_options.(CutestOptionKey.EXCLUDELIST.value)) || ~all(cellfun(@ischarstr, cutest_options.(CutestOptionKey.EXCLUDELIST.value)))
            error("MATLAB:checkValidityCutestOptions:excludelistNotCellOfcharstr", "The field 'excludelist' of options should be a cell array of strings or chars.");
        end
        % Convert to cell array of chars.
        cutest_options.(CutestOptionKey.EXCLUDELIST.value) = cellfun(@char, cutest_options.(CutestOptionKey.EXCLUDELIST.value), 'UniformOutput', false);
        % Convert to a row vector if not.
        if size(cutest_options.(CutestOptionKey.EXCLUDELIST.value), 1) > 1
            cutest_options.(CutestOptionKey.EXCLUDELIST.value) = cutest_options.(CutestOptionKey.EXCLUDELIST.value)';
        end
    end

end