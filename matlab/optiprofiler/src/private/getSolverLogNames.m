function solver_log_names = getSolverLogNames(solver_names, problem_name_width)
%GETSOLVERLOGNAMES returns short solver names for aligned terminal logs.

    log_width = 104;
    solver_name_width = max(cellfun(@length, solver_names));
    if standardConstrainedResultLogLength(problem_name_width, solver_name_width) <= log_width
        solver_log_names = solver_names;
        return;
    end

    solver_log_names = cell(size(solver_names));
    for i_solver = 1:numel(solver_names)
        solver_log_names{i_solver} = sprintf('solver%d', i_solver);
    end
end

function line_length = standardConstrainedResultLogLength(problem_name_width, solver_name_width)
    prefix = sprintf('[%-7s] ', 'INFO');
    problem_name = repmat('P', 1, max(1, problem_name_width));
    solver_name = repmat('S', 1, max(1, solver_name_width));
    format_info = sprintf('Output result for %%-%ds with %%-%ds (run %%2d/%%2d): f = %%10.4e, maxcv = %%10.4e.', problem_name_width, solver_name_width);
    line = [prefix, sprintf(format_info, problem_name, solver_name, 1, 1, -9.5, 1.1102e-16)];
    line_length = numel(line);
end
