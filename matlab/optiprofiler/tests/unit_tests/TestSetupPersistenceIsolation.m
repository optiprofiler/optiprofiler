classdef TestSetupPersistenceIsolation < matlab.unittest.TestCase
    methods (Test)

        function testSetupAndUninstallUsePathdefOverride(testCase)
            state = isolateSetupState(testCase);
            pathdef_file = fullfile(state.root, 'pathdef.m');
            startup_file = fullfile(state.root, 'startup.m');
            setenv('OPTIPROFILER_MATLAB_PROBLEM_LIBRARY_PATHDEF', pathdef_file);
            setenv('OPTIPROFILER_MATLAB_PROBLEM_LIBRARY_STARTUP', startup_file);

            runSetup(state.repository_root, state.library_root);

            testCase.verifyEqual(exist(pathdef_file, 'file'), 2);
            testCase.verifyEqual(exist(startup_file, 'file'), 0);
            verifyDefaultStartupUnchanged(testCase, state);

            runUninstall(state.repository_root);
            path(state.original_path);

            testCase.verifyFalse(contains( ...
                fileread(pathdef_file), state.repository_root));
            verifyDefaultStartupUnchanged(testCase, state);
        end

        function testSetupAndUninstallUseStartupOverride(testCase)
            state = isolateSetupState(testCase);
            pathdef_file = fullfile(state.root, 'missing', 'pathdef.m');
            startup_file = fullfile(state.root, 'startup.m');
            setenv('OPTIPROFILER_MATLAB_PROBLEM_LIBRARY_PATHDEF', pathdef_file);
            setenv('OPTIPROFILER_MATLAB_PROBLEM_LIBRARY_STARTUP', startup_file);

            runSetup(state.repository_root, state.library_root);

            testCase.verifyEqual(exist(startup_file, 'file'), 2);
            testCase.verifySubstring( ...
                fileread(startup_file), state.repository_root);
            verifyDefaultStartupUnchanged(testCase, state);

            runUninstall(state.repository_root);
            path(state.original_path);

            testCase.verifyFalse(contains( ...
                fileread(startup_file), state.repository_root));
            verifyDefaultStartupUnchanged(testCase, state);
        end

    end
end


function state = isolateSetupState(testCase)
    test_dir = fileparts(mfilename('fullpath'));
    [status, attributes] = fileattrib( ...
        fullfile(test_dir, '..', '..', '..', '..'));
    assert(status, 'Cannot resolve the repository root.');
    state.repository_root = attributes.Name;
    state.root = tempname;
    mkdir(state.root);
    state.library_root = fullfile(state.root, 'libraries');
    mkdir(state.library_root);
    state.original_directory = pwd;
    state.original_path = path;
    state.original_registry = getenv( ...
        'OPTIPROFILER_MATLAB_PROBLEM_LIBRARY_REGISTRY');
    state.original_pathdef = getenv( ...
        'OPTIPROFILER_MATLAB_PROBLEM_LIBRARY_PATHDEF');
    state.original_startup = getenv( ...
        'OPTIPROFILER_MATLAB_PROBLEM_LIBRARY_STARTUP');
    state.default_startup = getDefaultStartupFile();
    [state.default_startup_exists, state.default_startup_text] = ...
        readOptionalFile(state.default_startup);

    setenv('OPTIPROFILER_MATLAB_PROBLEM_LIBRARY_REGISTRY', ...
        fullfile(state.root, 'problem_libraries.mat'));
    testCase.addTeardown(@() restoreSetupState(state));
end


function runSetup(repository_root, library_root)
    original_directory = pwd;
    cleanup = onCleanup(@() cd(original_directory));
    cd(repository_root);
    options = struct( ...
        'install_matcutest', false, ...
        'install_solar', false, ...
        'problem_library_root', library_root);
    setup(options);
    clear cleanup
end


function runUninstall(repository_root)
    original_directory = pwd;
    cleanup = onCleanup(@() cd(original_directory));
    cd(repository_root);
    setup('uninstall');
    clear cleanup
end


function verifyDefaultStartupUnchanged(testCase, state)
    [exists_now, text_now] = readOptionalFile(state.default_startup);
    testCase.verifyEqual(exists_now, state.default_startup_exists);
    testCase.verifyEqual(text_now, state.default_startup_text);
end


function startup_file = getDefaultStartupFile()
    matlab_userpath = userpath;
    if isempty(matlab_userpath)
        startup_file = '';
    else
        startup_file = fullfile(matlab_userpath, 'startup.m');
    end
end


function [file_exists, text] = readOptionalFile(filename)
    file_exists = ~isempty(filename) && exist(filename, 'file') == 2;
    if file_exists
        text = fileread(filename);
    else
        text = '';
    end
end


function restoreSetupState(state)
    cd(state.original_directory);
    path(state.original_path);
    setenv('OPTIPROFILER_MATLAB_PROBLEM_LIBRARY_REGISTRY', ...
        state.original_registry);
    setenv('OPTIPROFILER_MATLAB_PROBLEM_LIBRARY_PATHDEF', ...
        state.original_pathdef);
    setenv('OPTIPROFILER_MATLAB_PROBLEM_LIBRARY_STARTUP', ...
        state.original_startup);
    if exist(state.root, 'dir') == 7
        rmdir(state.root, 's');
    end
end
