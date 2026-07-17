function writeProblemLibraryRegistry(registry_file, registrations)
%WRITEPROBLEMLIBRARYREGISTRY atomically persists raw registrations.

    if ~isstruct(registrations)
        error("MATLAB:problemLibraryRegistry:invalidFile", ...
            "Problem-library registrations must be stored as a structure array.");
    end

    registry_dir = fileparts(registry_file);
    if exist(registry_dir, 'dir') ~= 7
        [status, message] = mkdir(registry_dir);
        if ~status
            error("MATLAB:problemLibraryRegistry:registryDirectory", ...
                "Cannot create the problem-library registry directory '%s': %s", ...
                registry_dir, message);
        end
    end

    temporary_file = [tempname(registry_dir), '.mat'];
    cleanup = onCleanup(@() deleteIfExists(temporary_file));
    try
        save(temporary_file, 'registrations');
        [status, message] = movefile(temporary_file, registry_file, 'f');
    catch ME
        error("MATLAB:problemLibraryRegistry:registryWrite", ...
            "Cannot update the problem-library registry '%s': %s", ...
            registry_file, ME.message);
    end
    if ~status
        error("MATLAB:problemLibraryRegistry:registryWrite", ...
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
