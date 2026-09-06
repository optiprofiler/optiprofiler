function run_matlab_runtime_robustness(repository_root, output_root)
%RUN_MATLAB_RUNTIME_ROBUSTNESS Exact-SHA history and bundled S2MPJ/setup gate.
% Reloads synthetic histories and checks setup; no optimizer or pool runs.
% All selected cases are mandatory on this platform: an unexpected assumption
% or incomplete method fails the gate, even when its Failed flag is false.
    old_dir=pwd; old_path=path;
    variables={'OPTIPROFILER_MATLAB_PROBLEM_LIBRARY_REGISTRY', ...
        'OPTIPROFILER_MATLAB_PROBLEM_LIBRARY_PATHDEF', ...
        'OPTIPROFILER_MATLAB_PROBLEM_LIBRARY_STARTUP'};
    old_values=cellfun(@getenv,variables,'UniformOutput',false);
    state=tempname; mkdir(state);
    cleanup=onCleanup(@() restoreState(old_dir,old_path,variables,old_values,state));
    setenv(variables{1},fullfile(state,'registry.mat'));
    setenv(variables{2},fullfile(state,'pathdef.m'));
    setenv(variables{3},fullfile(state,'startup.m'));
    if ~isfolder(output_root), mkdir(output_root); end
    [ok,attributes]=fileattrib(output_root); assert(ok); output_root=attributes.Name;
    cd(repository_root);
    [git_status,git_sha]=system('git rev-parse HEAD');
    [dirty_status,dirty]=system('git status --porcelain --untracked-files=no');
    expected_sha=getenv('GITHUB_SHA');
    if git_status~=0, git_sha=''; end
    metadata=struct('github_sha',expected_sha,'git_sha',strtrim(git_sha), ...
        'git_available',git_status==0,'tracked_state',strtrim(dirty), ...
        'platform',computer('arch'),'matlab_release',version('-release'), ...
        'exact_ci_checkout',~isempty(expected_sha) && git_status==0 && ...
        strcmp(strtrim(git_sha),expected_sha) && dirty_status==0 && isempty(strtrim(dirty)));
    writeJson(fullfile(output_root,'source.json'),metadata);
    if ~isempty(expected_sha)
        assert(git_status==0 && strcmp(strtrim(git_sha),expected_sha), ...
            'The checked-out commit differs from the triggering SHA.');
        assert(dirty_status==0 && isempty(strtrim(dirty)), ...
            'The exact-SHA CI checkout has modified tracked files or submodules.');
    end

    source=fullfile(repository_root,'matlab','optiprofiler','src');
    tests=fullfile(repository_root,'matlab','optiprofiler','tests','unit_tests');
    addpath(source,tests);
    suite=[testsuite(fullfile(tests,'TestRuntimeRobustness.m')), ...
        testsuite(fullfile(tests,'TestHistoryExtremes.m'))];
    setup_suite=testsuite(fullfile(tests,'TestSetupPathOwnership.m'));
    % These setup methods are cross-platform and isolate their own state.
    % Default empty-userpath fallback requires a read-only default pathdef;
    % its environment-dependent test is validated separately, not skipped here.
    methods={'testBorrowedPathsAndRepeatedSetup', ...
        'testFallbackPreservesUnrelatedStartupBytes','testUnknownOwnershipDoesNotGuess', ...
        'testPersistenceTargetChangeFails','testCorruptLedgerFailsBeforeCleanup', ...
        'testRelativeOwnershipLocationIsRejected','testPersistedUserPathSurvivesDifferentSessionPath', ...
        'testFailedStartupAppendCanRetry','testExistingUserStartupEntryIsNotClaimed', ...
        'testSidecarWriteFailureDoesNotMutatePaths','testRemovalWriteFailurePreservesBytesAndCanRetry', ...
        'testPersistedPathIsBorrowedEvenWhenNotCurrentlyLoaded', ...
        'testExplicitTargetCannotDegradeToSessionOnly','testSetupOwnsS2RuntimePaths', ...
        'testSetupBorrowsExistingS2RuntimePaths'};
    if ispc
        methods=[methods,{'testWindowsRelativePathsRejected','testWindowsAbsolutePathForms'}];
    else
        methods=[methods,{'testPosixRejectsWindowsDrivePaths'}];
    end
    for k=1:numel(methods)
        selected=setup_suite(strcmp({setup_suite.Name},['TestSetupPathOwnership/',methods{k}]));
        assert(numel(selected)==1,'A mandatory setup method is missing: %s',methods{k});
        suite=[suite,selected]; %#ok<AGROW>
    end
    writeJson(fullfile(output_root,'selected-tests.json'),{suite.Name});
    runner=matlab.unittest.TestRunner.withTextOutput( ...
        'OutputDetail',matlab.unittest.Verbosity.Detailed);
    runner.addPlugin(matlab.unittest.plugins.XMLPlugin.producingJUnitFormat( ...
        fullfile(output_root,'junit.xml')));
    % benchmark changes diary targets; capture the runner output independently.
    log=evalc('results=runner.run(suite);');
    fprintf('%s',log);
    writeText(fullfile(output_root,'matlab.log'),log);
    disp(table(results));
    save(fullfile(output_root,'results.mat'),'results');
    records=arrayfun(@(r) struct('name',r.Name,'passed',r.Passed, ...
        'failed',r.Failed,'incomplete',r.Incomplete,'duration',r.Duration),results);
    writeJson(fullfile(output_root,'results.json'),records);
    [post_status,post_dirty]=system('git status --porcelain --untracked-files=no');
    writeJson(fullfile(output_root,'post-source.json'), ...
        struct('git_available',post_status==0,'tracked_state',strtrim(post_dirty)));
    if ~isempty(expected_sha)
        assert(post_status==0 && isempty(strtrim(post_dirty)), ...
            'Tests modified tracked source files or submodules in the CI checkout.');
    end
    assert(numel(results)==numel(suite) && all([results.Passed]) && ...
        ~any([results.Failed]) && ~any([results.Incomplete]), ...
        'Every selected runtime/history/setup case must complete and pass.');
end

function writeJson(filename,value)
    writeText(filename,jsonencode(value));
end

function writeText(filename,value)
    fid=fopen(filename,'w'); assert(fid>=0,'Cannot write evidence: %s',filename);
    cleanup=onCleanup(@() fclose(fid));
    fprintf(fid,'%s\n',value);
end

function restoreState(old_dir,old_path,variables,values,state)
    cd(old_dir); path(old_path);
    for k=1:numel(variables), setenv(variables{k},values{k}); end
    if isfolder(state), rmdir(state,'s'); end
end
