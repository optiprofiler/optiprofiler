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
    library = materializeProblemLibrary(registration);

    registry_file = getProblemLibraryRegistryFile();
    registrations = loadRegisteredProblemLibraries(registry_file);
    matching = strcmp({registrations.name}, registration.name);
    registrations(matching) = [];
    registrations(end + 1) = registration;

    registry_dir = fileparts(registry_file);
    if exist(registry_dir, 'dir') ~= 7
        [status, message] = mkdir(registry_dir);
        if ~status
            error("MATLAB:registerProblemLibrary:registryDirectory", ...
                "Cannot create the problem-library registry directory '%s': %s", ...
                registry_dir, message);
        end
    end

    temporary_file = [tempname(registry_dir), '.mat'];
    cleanup = onCleanup(@() deleteIfExists(temporary_file));
    save(temporary_file, 'registrations');
    [status, message] = movefile(temporary_file, registry_file, 'f');
    if ~status
        error("MATLAB:registerProblemLibrary:registryWrite", ...
            "Cannot update the problem-library registry '%s': %s", ...
            registry_file, message);
    end
    clear cleanup
end


function deleteIfExists(filename)
    if exist(filename, 'file') == 2
        delete(filename);
    end
end
