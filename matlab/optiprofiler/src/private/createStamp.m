function stamp = createStamp(solver_names, problem_options, feature_stamp, time_stamp)
%CREATESTAMP creates a stamp for the current experiment.
% The stamp is a string that contains the solver names, problem settings, feature stamp, and time stamp.

    % Get the solver stamp.
    for i = 1:length(solver_names)
        solver_names{i} = regexprep(solver_names{i}, '[^a-zA-Z0-9_]', '_');
    end
    solver_stamp = strjoin(solver_names, '_');

    % Get the problem stamp.
    % If options.ptype is 'u', then it is in the form 'ptype_mindim_maxdim'.
    % If options.ptype contains others, then it is in the form 'ptype_mindim_maxdim_mincon_maxcon'.
    ptype = problem_options.(ProblemOptionKey.PTYPE.value);
    mindim = problem_options.(ProblemOptionKey.MINDIM.value);
    maxdim = problem_options.(ProblemOptionKey.MAXDIM.value);
    mincon = problem_options.(ProblemOptionKey.MINCON.value);
    maxcon = problem_options.(ProblemOptionKey.MAXCON.value);
    if strcmp(ptype, 'u')
        problem_stamp = sprintf('%s_%d_%d', ptype, mindim, maxdim);
    else
        problem_stamp = sprintf('%s_%d_%d_%d_%d', ptype, mindim, maxdim, mincon, maxcon);
    end

    % Create the final stamp.
    stamp = sprintf('%s_%s_%s_%s', solver_stamp, problem_stamp, feature_stamp, time_stamp);
end