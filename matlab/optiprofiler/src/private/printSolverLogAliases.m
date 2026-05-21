function printSolverLogAliases(solver_names, solver_log_names)
%PRINTSOLVERLOGALIASES prints solver aliases used in terminal logs.

    if isequal(solver_names, solver_log_names)
        return;
    end

    printOptiProfilerMessage('INFO', 'Solver aliases used in the log:');
    max_alias_length = max(cellfun(@length, solver_log_names));
    format_alias = sprintf('%%-%ds = %%s', max_alias_length);
    for i_solver = 1:numel(solver_names)
        printOptiProfilerMessage('INFO', sprintf(format_alias, solver_log_names{i_solver}, solver_names{i_solver}));
    end
end
