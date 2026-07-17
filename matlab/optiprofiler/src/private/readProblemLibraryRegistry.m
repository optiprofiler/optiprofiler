function registrations = readProblemLibraryRegistry(registry_file)
%READPROBLEMLIBRARYREGISTRY reads raw registrations without resolving roots.

    registrations = emptyRegistrations();
    if exist(registry_file, 'file') ~= 2
        return
    end

    registry = load(registry_file);
    if ~isfield(registry, 'registrations') || ~isstruct(registry.registrations)
        error("MATLAB:problemLibraryRegistry:invalidFile", ...
            "The problem-library registry '%s' is invalid.", registry_file);
    end
    registrations = registry.registrations;
    expected_fields = fieldnames(emptyRegistration());
    actual_fields = fieldnames(registrations);
    if ~isempty(setxor(expected_fields, actual_fields))
        error("MATLAB:problemLibraryRegistry:invalidFile", ...
            "The problem-library registry '%s' has an invalid schema.", ...
            registry_file);
    end
end


function registrations = emptyRegistrations()
    registrations = repmat(emptyRegistration(), 0, 1);
end


function registration = emptyRegistration()
    registration = struct( ...
        'name', '', ...
        'api_version', 1, ...
        'role', '', ...
        'root', '', ...
        'select_function', '', ...
        'load_function', '', ...
        'collect_info_function', '', ...
        'check_available_function', '', ...
        'platforms', {{}});
end
