function registrations = loadRegisteredProblemLibraries(registry_file)
%LOADREGISTEREDPROBLEMLIBRARIES reads locally registered MATLAB libraries.

    registrations = emptyProblemLibraryRegistrations();
    if exist(registry_file, 'file') ~= 2
        return
    end
    registry = load(registry_file, 'registrations');
    if ~isfield(registry, 'registrations') || ~isstruct(registry.registrations)
        error("MATLAB:problemLibraryRegistry:invalidFile", ...
            "The problem-library registry '%s' is invalid.", registry_file);
    end
    stored = registry.registrations;
    for i_registration = 1:numel(stored)
        try
            registrations(end + 1) = normalizeProblemLibraryRegistration(stored(i_registration));
        catch ME
            % A removed optional checkout must not make bundled S2MPJ
            % unavailable. Treat stale roots as uninstalled; preserve errors
            % for malformed registry metadata.
            if strcmp(ME.identifier, 'MATLAB:problemLibraryRegistry:rootNotFound')
                continue
            end
            rethrow(ME)
        end
    end
end


function registrations = emptyProblemLibraryRegistrations()
    registrations = repmat(struct( ...
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
