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
    original_registry = getenv( ...
        'OPTIPROFILER_MATLAB_PROBLEM_LIBRARY_REGISTRY');
    original_pathdef = getenv( ...
        'OPTIPROFILER_MATLAB_PROBLEM_LIBRARY_PATHDEF');
    original_startup = getenv( ...
        'OPTIPROFILER_MATLAB_PROBLEM_LIBRARY_STARTUP');
    cleanup_environment = onCleanup(@() restore_environment( ...
        original_registry, original_pathdef, original_startup));
    setenv('OPTIPROFILER_MATLAB_PROBLEM_LIBRARY_REGISTRY', ...
        fullfile(registry_root, 'problem_libraries.mat'));
    setenv('OPTIPROFILER_MATLAB_PROBLEM_LIBRARY_PATHDEF', ...
        fullfile(registry_root, 'pathdef.m'));
    setenv('OPTIPROFILER_MATLAB_PROBLEM_LIBRARY_STARTUP', ...
        fullfile(registry_root, 'startup.m'));

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
    smoke_eval_report(problem);

    setup uninstall;
    archive_paths = strsplit(path, pathsep);
    assert(~any(startsWith(archive_paths, archive_root)), ...
        'setup uninstall did not remove the extracted engine path.');
end


function smoke_eval_report(problem)
%SMOKE_EVAL_REPORT The opt-in report must work from the extracted ZIP alone.
% The archive ships no tests and no schemas, so this checks the contract
% facts directly: two JSON files, a completed status, a companion receipt
% whose SHA256 matches the file, UTF-8 bytes, and no figure or archive output
% in score_only mode.
    report_root = tempname;
    mkdir(report_root);
    cleanup_report = onCleanup(@() remove_directory(report_root)); %#ok<NASGU>
    options = struct('problem', problem, 'score_only', true, 'silent', true, ...
        'n_runs', 1, 'max_eval_factor', 5, 'solver_names', {{'fminsearch_1', 'fminsearch_2'}}, ...
        'report_path', fullfile(report_root, 'smoke.json'));
    scores = evalc_benchmark(options);
    assert(numel(scores) == 2 && all(isfinite(scores)), 'The report run must return two finite scores.');
    produced = dir(fullfile(report_root, '*'));
    produced = {produced(~[produced.isdir]).name};
    assert(isequal(sort(produced), {'smoke.json', 'smoke.plot_data.json'}), ...
        'score_only must produce exactly the two requested JSON files.');
    report = jsondecode(read_utf8(fullfile(report_root, 'smoke.json')));
    assert(strcmp(report.schema, 'optiprofiler.eval_report/1') && strcmp(report.status, 'completed'));
    assert(strcmp(report.producer.language, 'matlab') && strcmp(report.problems.library, 'user'));
    assert(strcmp(report.problems.name, 'BEALE') && numel(report.problems.runs) == 2);
    companion = fullfile(report_root, report.plot_data.path);
    assert(strcmp(report.plot_data.status, 'completed'), 'The companion receipt must be complete when hashing is available.');
    assert(strcmp(report.plot_data.sha256, sha256_of(companion)), 'The companion SHA256 receipt must match the file.');
    detail = jsondecode(read_utf8(companion));
    assert(strcmp(detail.evaluation_id, report.evaluation_id) && numel(detail.histories) == 2);
end


function scores = evalc_benchmark(options)
    % Keep the smoke log readable: the benchmark banner is not evidence here.
    [~, scores] = evalc('benchmark({@smoke_solver_coarse, @smoke_solver_fine}, options)');
end


function x = smoke_solver_coarse(fun, x0)
    x = fminsearch(fun, x0, optimset('MaxFunEvals', 20, 'Display', 'off'));
end


function x = smoke_solver_fine(fun, x0)
    x = fminsearch(fun, x0, optimset('MaxFunEvals', 60, 'TolX', 1e-8, 'TolFun', 1e-8, 'Display', 'off'));
end


function text = read_utf8(path)
    fid = fopen(path, 'rb');
    assert(fid >= 0, 'Cannot open %s', path);
    guard = onCleanup(@() fclose(fid)); %#ok<NASGU>
    text = native2unicode(reshape(fread(fid, '*uint8'), 1, []), 'UTF-8');
end


function value = sha256_of(path)
    digest = java.security.MessageDigest.getInstance('SHA-256');
    fid = fopen(path, 'rb');
    assert(fid >= 0, 'Cannot open %s', path);
    guard = onCleanup(@() fclose(fid)); %#ok<NASGU>
    while ~feof(fid)
        block = fread(fid, 1048576, '*uint8');
        if ~isempty(block), digest.update(block); end
    end
    value = lower(reshape(dec2hex(typecast(digest.digest(), 'uint8'), 2).', 1, []));
end


function restore_environment(registry, pathdef, startup)
    setenv('OPTIPROFILER_MATLAB_PROBLEM_LIBRARY_REGISTRY', registry);
    setenv('OPTIPROFILER_MATLAB_PROBLEM_LIBRARY_PATHDEF', pathdef);
    setenv('OPTIPROFILER_MATLAB_PROBLEM_LIBRARY_STARTUP', startup);
end


function remove_directory(directory)
    if exist(directory, 'dir') == 7
        rmdir(directory, 's');
    end
end
