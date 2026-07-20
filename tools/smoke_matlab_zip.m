function smoke_matlab_zip(archive_root)
%SMOKE_MATLAB_ZIP Test setup and the default benchmark from a clean ZIP.

    arguments
        archive_root (1, :) char
    end

    archive_root = char(java.io.File(archive_root).getCanonicalPath());
    assert(exist(fullfile(archive_root, 'setup.m'), 'file') == 2, ...
        'The extracted MATLAB archive does not contain setup.m.');
    assert(exist(fullfile(archive_root, 'python'), 'dir') ~= 7, ...
        'The MATLAB-only archive unexpectedly contains Python sources.');
    assert(exist(fullfile(archive_root, '.github'), 'dir') ~= 7, ...
        'The MATLAB-only archive unexpectedly contains workflows.');

    original_directory = pwd;
    original_path = path;
    cleanup_directory = onCleanup(@() cd(original_directory));
    cleanup_path = onCleanup(@() path(original_path));
    restoredefaultpath;
    registry_root = tempname;
    library_root = tempname;
    mkdir(registry_root);
    mkdir(library_root);
    cleanup_registry = onCleanup(@() remove_directory(registry_root));
    cleanup_libraries = onCleanup(@() remove_directory(library_root));
    setenv('OPTIPROFILER_MATLAB_PROBLEM_LIBRARY_REGISTRY', ...
        fullfile(registry_root, 'problem_libraries.mat'));

    cd(archive_root);
    setup(struct('install_matcutest', false, 'install_solar', false, ...
        'problem_library_root', library_root));
    assert(startsWith(which('benchmark'), archive_root), ...
        'benchmark was not loaded from the clean archive.');
    assert(startsWith(which('s2mpj_load'), archive_root), ...
        'S2MPJ was not loaded from the clean archive.');

    library = resolveProblemLibrary('s2mpj');
    problem = library.load('BEALE');
    assert(strcmp(problem.name, 'BEALE'));
    assert(isequal(problem.x0, [1; 1]));
    assert(isfinite(problem.fun(problem.x0)));
    testOptiProfiler();

    setup uninstall;
    archive_paths = strsplit(path, pathsep);
    assert(~any(startsWith(archive_paths, archive_root)), ...
        'setup uninstall did not remove the extracted engine path.');
end


function remove_directory(directory)
    if exist(directory, 'dir') == 7
        rmdir(directory, 's');
    end
end
