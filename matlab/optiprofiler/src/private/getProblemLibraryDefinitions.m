function definitions = getProblemLibraryDefinitions()
%GETPROBLEMLIBRARYDEFINITIONS returns locked, bundled, and local entries.

    [lock, repository_root] = readMatlabProblemLibraryLock();
    registered = loadRegisteredProblemLibraries(getProblemLibraryRegistryFile());
    definitions = emptyDefinitions();

    locked = lock.libraries;
    for i_library = 1:numel(locked)
        entry = getEntry(locked, i_library);
        definition = struct( ...
            'name', char(entry.name), ...
            'api_version', double(entry.api_version), ...
            'role', char(entry.role), ...
            'root', '', ...
            'select_function', char(entry.select_function), ...
            'load_function', char(entry.load_function), ...
            'collect_info_function', optionalString(entry.collect_info_function), ...
            'check_available_function', '', ...
            'platforms', {cellstr(entry.platforms)});

        matching = find(strcmp({registered.name}, definition.name), 1);
        if strcmp(definition.role, 'bundled-default')
            % Bundled defaults are always resolved from the locked gitlink.
            % A manually edited registry must not redirect their source root.
            definition.root = fullfile(repository_root, char(entry.gitlink));
        elseif ~isempty(matching)
            registered_entry = registered(matching);
            if strcmp(registered_entry.select_function, definition.select_function) && ...
                    strcmp(registered_entry.load_function, definition.load_function)
                definition.root = registered_entry.root;
                definition.check_available_function = registered_entry.check_available_function;
            end
        end
        if isempty(definition.root)
            definitions(end + 1) = definition;
        else
            definitions(end + 1) = normalizeProblemLibraryRegistration(definition);
        end
    end

    definitions(end + 1) = normalizeProblemLibraryRegistration(bundledExampleProblemLibrary(repository_root));

    locked_names = {definitions.name};
    for i_registration = 1:numel(registered)
        if ~ismember(registered(i_registration).name, locked_names)
            definitions(end + 1) = registered(i_registration);
        end
    end
end


function entry = getEntry(entries, index)
    if iscell(entries)
        entry = entries{index};
    else
        entry = entries(index);
    end
end


function value = optionalString(value)
    if isempty(value)
        value = '';
    else
        value = char(value);
    end
end


function definitions = emptyDefinitions()
    definitions = repmat(struct( ...
        'name', '', ...
        'api_version', 1, ...
        'role', '', ...
        'root', '', ...
        'select_function', '', ...
        'load_function', '', ...
        'collect_info_function', '', ...
        'check_available_function', '', ...
        'platforms', {{}}), 0, 1);
end
