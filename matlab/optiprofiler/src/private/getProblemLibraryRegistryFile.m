function registry_file = getProblemLibraryRegistryFile()
%GETPROBLEMLIBRARYREGISTRYFILE returns the local MATLAB registry file.

    registry_file = getenv('OPTIPROFILER_MATLAB_PROBLEM_LIBRARY_REGISTRY');
    if isempty(registry_file)
        registry_file = fullfile(prefdir, 'optiprofiler', 'problem_libraries.mat');
    end
    registry_file = char(registry_file);
    if isempty(fileparts(registry_file))
        registry_file = fullfile(pwd, registry_file);
    end
end
