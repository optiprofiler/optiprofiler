function registration = unregisterProblemLibrary(name)
%UNREGISTERPROBLEMLIBRARY removes one persistent MATLAB library registration.
%
%   REGISTRATION = UNREGISTERPROBLEMLIBRARY(NAME) removes the external or
%   custom registration named NAME and removes its root from the current and
%   persisted MATLAB path. It does not delete source files, runtimes, caches,
%   compiled artifacts, or other user data.
%
%   Bundled S2MPJ and the bundled custom example are part of OptiProfiler and
%   cannot be unregistered.

    name = normalizeName(name);
    if ismember(name, {'s2mpj', 'custom'})
        error("MATLAB:unregisterProblemLibrary:bundledProvider", ...
            "The bundled problem library '%s' cannot be unregistered.", name);
    end

    registry_file = getProblemLibraryRegistryFile();
    registrations = readProblemLibraryRegistry(registry_file);
    names = storedRegistrationNames(registrations, registry_file);
    matching = find(strcmp(names, name));
    if isempty(matching)
        error("MATLAB:unregisterProblemLibrary:notRegistered", ...
            "Problem library '%s' has no persistent registration.", name);
    elseif numel(matching) > 1
        error("MATLAB:problemLibraryRegistry:invalidFile", ...
            "The problem-library registry '%s' contains duplicate entries for '%s'.", ...
            registry_file, name);
    end

    registration = registrations(matching);
    if ~isfield(registration, 'root') || ...
            ~ischarstr(registration.root) || ...
            ~isscalar(string(registration.root))
        error("MATLAB:problemLibraryRegistry:invalidFile", ...
            "The problem-library registry '%s' contains an invalid root for '%s'.", ...
            registry_file, name);
    end
    root = strtrim(char(registration.root));
    if isempty(root)
        error("MATLAB:problemLibraryRegistry:invalidFile", ...
            "The problem-library registry '%s' contains an empty root for '%s'.", ...
            registry_file, name);
    end

    registrations(matching) = [];
    writeProblemLibraryRegistry(registry_file, registrations);
    removeProblemLibraryPath(root);
end


function name = normalizeName(name)
    if ~ischarstr(name) || ~isscalar(string(name))
        error("MATLAB:unregisterProblemLibrary:nameNotValid", ...
            "The problem-library name must be a char or string scalar.");
    end
    name = lower(strtrim(char(name)));
    if isempty(regexp(name, '^[a-z][a-z0-9_]*$', 'once'))
        error("MATLAB:unregisterProblemLibrary:nameNotValid", ...
            "Invalid problem-library name '%s'.", name);
    end
end


function names = storedRegistrationNames(registrations, registry_file)
    if isempty(registrations)
        names = {};
        return
    end
    if ~isfield(registrations, 'name')
        error("MATLAB:problemLibraryRegistry:invalidFile", ...
            "The problem-library registry '%s' contains an entry without a name.", ...
            registry_file);
    end
    names = cell(1, numel(registrations));
    for i_registration = 1:numel(registrations)
        stored_name = registrations(i_registration).name;
        if ~ischarstr(stored_name) || ~isscalar(string(stored_name))
            error("MATLAB:problemLibraryRegistry:invalidFile", ...
                "The problem-library registry '%s' contains an invalid name.", ...
                registry_file);
        end
        names{i_registration} = lower(strtrim(char(stored_name)));
    end
end
