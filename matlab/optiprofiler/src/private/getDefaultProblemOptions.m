function problem_options = getDefaultProblemOptions(problem_options)

    % Set default values for the problem_options.plibs if it is not set.
    if ~isfield(problem_options, ProblemOptionKey.PLIBS.value)
        problem_options.(ProblemOptionKey.PLIBS.value) = {'s2mpj'};
    end

    % If problem_options contains any field from {'ptype', 'mindim', 'maxdim', 'minb', 'maxb', 'minlcon', 'maxlcon', 'minnlcon', 'maxnlcon', 'mincon', 'maxcon', 'excludelist'}, then we will set the default values for all of them if they are not set and return.
    if isfield(problem_options, ProblemOptionKey.PTYPE.value) || isfield(problem_options, ProblemOptionKey.MINDIM.value) || isfield(problem_options, ProblemOptionKey.MAXDIM.value) || isfield(problem_options, ProblemOptionKey.MINB.value) || isfield(problem_options, ProblemOptionKey.MAXB.value) || isfield(problem_options, ProblemOptionKey.MINCON.value) || isfield(problem_options, ProblemOptionKey.MAXCON.value)
        if ~isfield(problem_options, ProblemOptionKey.PTYPE.value)
            problem_options.(ProblemOptionKey.PTYPE.value) = 'u';
        end
        if ~isfield(problem_options, ProblemOptionKey.MINDIM.value)
            problem_options.(ProblemOptionKey.MINDIM.value) = 1;
        end
        if ~isfield(problem_options, ProblemOptionKey.MAXDIM.value)
            problem_options.(ProblemOptionKey.MAXDIM.value) = problem_options.(ProblemOptionKey.MINDIM.value) + 1;
        end
        if ~isfield(problem_options, ProblemOptionKey.MINB.value)
            problem_options.(ProblemOptionKey.MINB.value) = 0;
        end
        if ~isfield(problem_options, ProblemOptionKey.MAXB.value)
            problem_options.(ProblemOptionKey.MAXB.value) = problem_options.(ProblemOptionKey.MINB.value) + 10;
        end
        if ~isfield(problem_options, ProblemOptionKey.MINLCON.value)
            problem_options.(ProblemOptionKey.MINLCON.value) = 0;
        end
        if ~isfield(problem_options, ProblemOptionKey.MAXLCON.value)
            problem_options.(ProblemOptionKey.MAXLCON.value) = problem_options.(ProblemOptionKey.MINLCON.value) + 10;
        end
        if ~isfield(problem_options, ProblemOptionKey.MINNLCON.value)
            problem_options.(ProblemOptionKey.MINNLCON.value) = 0;
        end
        if ~isfield(problem_options, ProblemOptionKey.MAXNLCON.value)
            problem_options.(ProblemOptionKey.MAXNLCON.value) = problem_options.(ProblemOptionKey.MINNLCON.value) + 10;
        end
        if ~isfield(problem_options, ProblemOptionKey.MINCON.value)
            problem_options.(ProblemOptionKey.MINCON.value) = min(problem_options.(ProblemOptionKey.MINLCON.value), problem_options.(ProblemOptionKey.MINNLCON.value));
        end
        if ~isfield(problem_options, ProblemOptionKey.MAXCON.value)
            problem_options.(ProblemOptionKey.MAXCON.value) = max(problem_options.(ProblemOptionKey.MAXLCON.value), problem_options.(ProblemOptionKey.MAXNLCON.value));
        end
        if ~isfield(problem_options, ProblemOptionKey.EXCLUDELIST.value)
            problem_options.(ProblemOptionKey.EXCLUDELIST.value) = {};
        end
        return;
    end
    
    % If problem_options now only contains the field 'plibs', then we will set the default values for {'ptype', 'mindim', 'maxdim', 'minb', 'maxb', 'mincon', 'maxcon', 'excludelist'} if they are not set.
    if numel(fieldnames(problem_options)) == 1
        problem_options.(ProblemOptionKey.PTYPE.value) = 'u';
        problem_options.(ProblemOptionKey.MINDIM.value) = 1;
        problem_options.(ProblemOptionKey.MAXDIM.value) = 2;
        problem_options.(ProblemOptionKey.MINB.value) = 0;
        problem_options.(ProblemOptionKey.MAXB.value) = 0;
        problem_options.(ProblemOptionKey.MINCON.value) = 0;
        problem_options.(ProblemOptionKey.MAXCON.value) = 0;
        problem_options.(ProblemOptionKey.EXCLUDELIST.value) = {};
    end 
end