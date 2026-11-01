function bundled = isBundledProblemLibrary(name)
%ISBUNDLEDPROBLEMLIBRARY whether the library NAME is part of OptiProfiler.
%   A library is bundled when its role, in matlab_problem_libraries.lock or
%   for the bundled example, is 'bundled-default' or 'bundled-example'.

    [lock, repository_root] = readMatlabProblemLibraryLock();
    entries = lock.libraries;
    names = cell(1, numel(entries));
    roles = cell(1, numel(entries));
    for i_entry = 1:numel(entries)
        if iscell(entries)
            entry = entries{i_entry};
        else
            entry = entries(i_entry);
        end
        names{i_entry} = char(entry.name);
        roles{i_entry} = char(entry.role);
    end
    example = bundledExampleProblemLibrary(repository_root);
    names{end + 1} = example.name;
    roles{end + 1} = example.role;
    bundled = any(strcmp(names, name) & ismember(roles, {'bundled-default', 'bundled-example'}));
end
