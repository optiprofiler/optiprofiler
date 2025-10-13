function stamp = createStamp(solver_names, problem_options, feature_stamp, time_stamp, path_out)
%CREATESTAMP creates a stamp for the current experiment.
% The stamp is a string that contains the solver names, problem settings, feature stamp, and time stamp.

    % Set the max length of the stamp.
    max_length = 100;
    if ispc
        max_dir_length = 250; % Windows max path length (260) minus some buffer
    else
        max_dir_length = 4000; % Unix max path length (4096) minus some buffer
    end
    % We want to avoid the path 'path_out/stamp/summary_stamp.pdf' exceeding the max path length.
    max_length = min(max_length, floor((max_dir_length - length(path_out) - length('summary_') - length('.pdf') - 1) / 2));

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
    stamp = time_stamp;
    if length(feature_stamp) + length(stamp) + 1 <= max_length
        stamp = sprintf('%s_%s', feature_stamp, stamp);
    end
    if length(problem_stamp) + length(stamp) + 1 <= max_length
        stamp = sprintf('%s_%s', problem_stamp, stamp);
    end
    if length(solver_stamp) + length(stamp) + 1 <= max_length
        stamp = sprintf('%s_%s', solver_stamp, stamp);
    end
end