function problem_options = getDefaultProblemOptions(problem_options)

    % Set default values for the problem_options.plibs if it is not set.
    if ~isfield(problem_options, ProblemOptionKey.PLIBS.value)
        problem_options.(ProblemOptionKey.PLIBS.value) = {'s2mpj'};
    end

    % Set default values for the problem_options.ptype if it is not set.
    if ~isfield(problem_options, ProblemOptionKey.PTYPE.value)
        problem_options.(ProblemOptionKey.PTYPE.value) = 'u';
    end

    % Set default values for the problem_options.mindim if it is not set.
    if ~isfield(problem_options, ProblemOptionKey.MINDIM.value)
        problem_options.(ProblemOptionKey.MINDIM.value) = 1;
    end

    % Set default values for the problem_options.maxdim if it is not set.
    if ~isfield(problem_options, ProblemOptionKey.MAXDIM.value)
        problem_options.(ProblemOptionKey.MAXDIM.value) = problem_options.(ProblemOptionKey.MINDIM.value) + 1;
    end

    % Set default values for the problem_options.minb if it is not set.
    if ~isfield(problem_options, ProblemOptionKey.MINB.value)
        problem_options.(ProblemOptionKey.MINB.value) = 0;
    end

    % Set default values for the problem_options.maxb if it is not set.
    if ~isfield(problem_options, ProblemOptionKey.MAXB.value)
        problem_options.(ProblemOptionKey.MAXB.value) = problem_options.(ProblemOptionKey.MINB.value) + 10;
    end

    % Set default values for the problem_options.minlcon if it is not set.
    if ~isfield(problem_options, ProblemOptionKey.MINLCON.value)
        problem_options.(ProblemOptionKey.MINLCON.value) = 0;
    end

    % Set default values for the problem_options.maxlcon if it is not set.
    if ~isfield(problem_options, ProblemOptionKey.MAXLCON.value)
        problem_options.(ProblemOptionKey.MAXLCON.value) = problem_options.(ProblemOptionKey.MINLCON.value) + 10;
    end

    % Set default values for the problem_options.minnlcon if it is not set.
    if ~isfield(problem_options, ProblemOptionKey.MINNLCON.value)
        problem_options.(ProblemOptionKey.MINNLCON.value) = 0;
    end

    % Set default values for the problem_options.maxnlcon if it is not set.
    if ~isfield(problem_options, ProblemOptionKey.MAXNLCON.value)
        problem_options.(ProblemOptionKey.MAXNLCON.value) = problem_options.(ProblemOptionKey.MINNLCON.value) + 10;
    end

    % Set default values for the problem_options.maxcon if it is not set.
    if ~isfield(problem_options, ProblemOptionKey.MINCON.value)
        problem_options.(ProblemOptionKey.MINCON.value) = min(problem_options.(ProblemOptionKey.MINLCON.value), problem_options.(ProblemOptionKey.MINNLCON.value));
    end

    % Set default values for the problem_options.maxcon if it is not set.
    if ~isfield(problem_options, ProblemOptionKey.MAXCON.value)
        problem_options.(ProblemOptionKey.MAXCON.value) = max(problem_options.(ProblemOptionKey.MAXLCON.value), problem_options.(ProblemOptionKey.MAXNLCON.value));
    end

    % Set default values for the problem_options.problem_names if it is not set.
    if ~isfield(problem_options, ProblemOptionKey.PROBLEM_NAMES.value)
        problem_options.(ProblemOptionKey.PROBLEM_NAMES.value) = {};
    end

    % Set default values for the problem_options.exclude_list if it is not set.
    if ~isfield(problem_options, ProblemOptionKey.EXCLUDELIST.value)
        problem_options.(ProblemOptionKey.EXCLUDELIST.value) = {};
    end
end