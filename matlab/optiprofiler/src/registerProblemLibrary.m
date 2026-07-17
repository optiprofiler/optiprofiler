function library = registerProblemLibrary(registration)
%REGISTERPROBLEMLIBRARY registers a MATLAB problem-library implementation.
%
%   LIBRARY = REGISTERPROBLEMLIBRARY(REGISTRATION) validates and persists a
%   scalar structure with the fields NAME, ROOT, SELECT_FUNCTION, and
%   LOAD_FUNCTION. Optional fields are API_VERSION, ROLE,
%   COLLECT_INFO_FUNCTION, CHECK_AVAILABLE_FUNCTION, and PLATFORMS.
%
%   The registry stores locations and canonical function names. OptiProfiler
%   therefore does not require an optional library to live inside its own
%   source tree. The environment variable
%   OPTIPROFILER_MATLAB_PROBLEM_LIBRARY_REGISTRY can override the registry
%   file location, which is useful for isolated tests and deployments.

    registration = normalizeProblemLibraryRegistration(registration);
    registration = applyLockedProblemLibraryContract(registration);
    if strcmp(registration.name, 'custom') || ...
            ismember(registration.role, {'bundled-default', 'bundled-example'})
        error("MATLAB:problemLibraryRegistry:bundledRegistration", ...
            "The bundled problem library '%s' cannot be registered externally.", ...
            registration.name);
    end

    registry_file = getProblemLibraryRegistryFile();
    registrations = readProblemLibraryRegistry(registry_file);
    names = storedRegistrationNames(registrations, registry_file);
    matching = find(strcmp(names, registration.name));
    if numel(matching) > 1
        error("MATLAB:problemLibraryRegistry:invalidFile", ...
            "The problem-library registry '%s' contains duplicate entries for '%s'.", ...
            registry_file, registration.name);
    elseif numel(matching) == 1
        stored = registrations(matching);
        try
            stored = normalizeProblemLibraryRegistration(stored);
            stored = applyLockedProblemLibraryContract(stored);
        catch ME
            if strcmp(ME.identifier, 'MATLAB:problemLibraryRegistry:rootNotFound')
                error("MATLAB:registerProblemLibrary:alreadyRegistered", ...
                    "Problem library '%s' is registered at a missing root. " + ...
                    "Call unregisterProblemLibrary('%s') before registering it again.", ...
                    registration.name, registration.name);
            end
            rethrow(ME)
        end
        if isequaln(stored, registration)
            library = materializeProblemLibrary(registration);
            return
        end
        error("MATLAB:registerProblemLibrary:alreadyRegistered", ...
            "Problem library '%s' is already registered with different metadata. " + ...
            "Call unregisterProblemLibrary('%s') before registering a replacement.", ...
            registration.name, registration.name);
    end

    was_on_path = isPathEntry(registration.root);
    library = materializeProblemLibrary(registration);
    registrations(end + 1) = registration;
    try
        writeProblemLibraryRegistry(registry_file, registrations);
    catch ME
        if ~was_on_path && isPathEntry(registration.root)
            rmpath(registration.root);
        end
        rethrow(ME)
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
        name = registrations(i_registration).name;
        if ~ischarstr(name) || ~isscalar(string(name))
            error("MATLAB:problemLibraryRegistry:invalidFile", ...
                "The problem-library registry '%s' contains an invalid name.", ...
                registry_file);
        end
        names{i_registration} = lower(strtrim(char(name)));
    end
end


function tf = isPathEntry(root)
    tf = any(strcmp(strsplit(path, pathsep), root));
end
