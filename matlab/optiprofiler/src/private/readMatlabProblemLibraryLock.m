function [lock, repository_root] = readMatlabProblemLibraryLock()
%READMATLABPROBLEMLIBRARYLOCK reads the repository MATLAB integration lock.

    private_dir = fileparts(mfilename('fullpath'));
    repository_root = fullfile(private_dir, '..', '..', '..', '..');
    [status, attributes] = fileattrib(repository_root);
    if status
        repository_root = attributes.Name;
    end
    lock_path = fullfile(repository_root, 'matlab_problem_libraries.lock');
    if exist(lock_path, 'file') ~= 2
        error("MATLAB:problemLibraryRegistry:lockNotFound", ...
            "Cannot find the MATLAB problem-library lock at '%s'.", lock_path);
    end
    lock = jsondecode(fileread(lock_path));
    if ~isstruct(lock) || ~isfield(lock, 'problem_library_api_version') || ...
            lock.problem_library_api_version ~= 1 || ~isfield(lock, 'libraries')
        error("MATLAB:problemLibraryRegistry:lockNotValid", ...
            "The MATLAB problem-library lock is invalid.");
    end
end
