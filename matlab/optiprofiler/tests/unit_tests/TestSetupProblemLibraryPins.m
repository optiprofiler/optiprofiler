classdef TestSetupProblemLibraryPins < matlab.unittest.TestCase
% setup brings a problem library to its locked commit or says that it cannot.
% Everything happens in a sandbox: a copy of the package with its own lock, local
% fixture repositories that stand in for GitHub, a private home and an isolated
% registry, pathdef and startup file. There is no network access, and the real
% MATLAB path and startup files are not touched.
    properties (TestParameter)
        % A problem_library_root that is a symbolic link, or that lies below one.
        linked_root = {'itself', 'below'};
        % The only function of a MatCUTEst of the user on the MATLAB path.
        partial_runtime = {'matcutest', 'macup', 'secup'};
    end

    methods (TestClassSetup)
        function requireGit(testCase)
            testCase.assumeTrue(isunix);
            [status, ~] = system('git --version');
            testCase.assumeEqual(status, 0);
        end
    end

    methods (Test)
        function branchOnlyPinIsFetchedAndCheckedOut(testCase)
            state = sandbox(testCase);
            base = make_remote(state, 'solar_matlab', solar_files());
            clone(state, 'solar_matlab', 'solar');
            pin = branch_commit(state, 'solar_matlab');
            testCase.assertNotEqual(pin, base);
            lock_commits(state, struct('solar', pin));
            run_setup(state, struct('install_solar', true));
            testCase.verifyEqual(head(state, 'solar'), pin);
            testCase.verifyTrue(is_registered(state, 'solar'));
            testCase.verifyTrue(on_path(fullfile(state.libraries, 'solar')));
        end

        function pinHeldLocallyNeedsNoRemote(testCase)
            state = sandbox(testCase);
            make_remote(state, 'solar_matlab', solar_files());
            pin = branch_commit(state, 'solar_matlab');
            clone(state, 'solar_matlab', 'solar');
            rmdir(fullfile(state.remotes, 'optiprofiler', 'solar_matlab.git'), 's');
            lock_commits(state, struct('solar', pin));
            output = run_setup(state, struct('install_solar', true));
            testCase.verifySubstring(output, 'checked out at pinned commit');
            testCase.verifyEqual(head(state, 'solar'), pin);
        end

        function repeatedSetupStaysAtThePin(testCase)
            state = sandbox(testCase);
            pin = make_remote(state, 'solar_matlab', solar_files());
            clone(state, 'solar_matlab', 'solar');
            lock_commits(state, struct('solar', pin));
            run_setup(state, struct('install_solar', true));
            output = run_setup(state, struct('install_solar', true));
            testCase.verifySubstring(output, 'already at pinned commit');
            testCase.verifyEqual(head(state, 'solar'), pin);
        end

        function unreachablePinOfTheBundledLibraryIsAnErrorAfterTheRest(testCase)
            state = sandbox(testCase);
            make_remote(state, 's2mpj_matlab', {'README.md', 'fixture'; ...
                fullfile('src', 'matlab_problems', 'keep.txt'), 'fixture'});
            s2 = fullfile(state.package, 'matlab', 'optiprofiler', 'problem_libs', 's2mpj');
            rmdir(s2, 's');
            git(state.root, sprintf('clone -q "%s" "%s"', remote_path(state, 's2mpj_matlab'), s2));
            before = rev_parse(s2);
            pin = make_remote(state, 'solar_matlab', solar_files());
            lock_commits(state, struct('s2mpj', unknown_commit(), 'solar', pin));
            [output, failure] = run_failing_setup(state, struct('install_solar', true));
            testCase.verifyEqual(failure.identifier, 'setup:PinnedCommitUnavailable');
            testCase.verifySubstring(output, ['WARNING: ', failure.identifier, ': ', failure.message]);
            testCase.verifyEqual(rev_parse(s2), before);
            verify_core_configured(testCase, state);
            testCase.verifyTrue(is_registered(state, 'solar'));
        end

        function otherFailureOfTheBundledLibraryIsNotCalledAPin(testCase)
            state = sandbox(testCase);
            make_remote(state, 's2mpj_matlab', {'README.md', 'fixture'; ...
                fullfile('src', 'matlab_problems', 'keep.txt'), 'fixture'});
            s2 = fullfile(state.package, 'matlab', 'optiprofiler', 'problem_libs', 's2mpj');
            rmdir(s2, 's');
            git(state.root, sprintf('clone -q "%s" "%s"', remote_path(state, 's2mpj_matlab'), s2));
            % A function beside setup.m takes precedence over the built-in one.
            testCase.applyFixture(matlab.unittest.fixtures.SuppressedWarningsFixture('MATLAB:dispatcher:nameConflict'));
            write_text(fullfile(state.package, 'system.m'), sprintf([ ...
                'function varargout = system(command)\n' ...
                'if contains(command, ''s2mpj'') && contains(command, ''rev-parse'')\n' ...
                '    error(''Fixture:Other'', ''a failure of another kind'');\n' ...
                'end\n' ...
                '[varargout{1:max(nargout, 1)}] = builtin(''system'', command);\n' ...
                'end\n']));
            [output, failure] = run_failing_setup(state, struct());
            testCase.verifyEqual(failure.identifier, 'Fixture:Other');
            testCase.verifySubstring(output, 'WARNING: Fixture:Other: a failure of another kind');
            testCase.verifyFalse(contains(output, 'pinned commit'));
            verify_core_configured(testCase, state);
        end

        function unreachablePinOfAnOptionalLibraryIsReportedAndNotRegistered(testCase)
            state = sandbox(testCase);
            base = make_remote(state, 'solar_matlab', solar_files());
            clone(state, 'solar_matlab', 'solar');
            lock_commits(state, struct('solar', unknown_commit()));
            output = run_setup(state, struct('install_solar', true));
            verify_core_configured(testCase, state);
            verify_reported(testCase, output, 'SOLAR MATLAB adapter', 'Cannot check out the pinned commit');
            testCase.verifyEqual(head(state, 'solar'), base);
            testCase.verifyFalse(is_registered(state, 'solar'));
            testCase.verifyFalse(on_path(fullfile(state.libraries, 'solar')));
        end

        function failedMatcutestCloneLeavesTheRestConfigured(testCase)
            testCase.assumeFalse(ismac);
            state = sandbox(testCase);
            pin = make_remote(state, 'solar_matlab', solar_files());
            lock_commits(state, struct('solar', pin));
            output = run_setup(state, struct('install_matcutest', true, 'install_solar', true));
            verify_core_configured(testCase, state);
            verify_reported(testCase, output, 'MatCUTEst', 'Failed to clone');
            testCase.verifyFalse(is_registered(state, 'matcutest'));
            testCase.verifyFalse(on_path(fullfile(state.libraries, 'matcutest')));
            testCase.verifyTrue(is_registered(state, 'solar'));
            testCase.verifyTrue(on_path(fullfile(state.libraries, 'solar')));
        end

        function adapterWithoutManagedEntryPointIsReported(testCase)
            testCase.assumeFalse(ismac);
            state = sandbox(testCase);
            files = matcutest_files('');
            pin = make_remote(state, 'matcutest', files(1:2, :));
            lock_commits(state, struct('matcutest', pin));
            output = run_setup(state, struct('install_matcutest', true));
            verify_core_configured(testCase, state);
            verify_reported(testCase, output, 'MatCUTEst', 'no managed runtime entry point');
            testCase.verifyFalse(is_registered(state, 'matcutest'));
            testCase.verifyFalse(on_path(fullfile(state.libraries, 'matcutest')));
            % The clone stays for the next attempt; it is simply not used.
            testCase.verifyEqual(head(state, 'matcutest'), pin);
        end

        function helpOfASimilarNameIsNotARuntime(testCase)
            testCase.assumeFalse(ismac);
            state = sandbox(testCase);
            % What help('matcutest') returns where the Bioinformatics Toolbox
            % is installed. A function beside setup.m takes precedence over help.
            write_text(fullfile(state.package, 'help.m'), sprintf(['function text = help(~)\n' ...
                'text = ''--- matcutest not found. Showing help for mattest instead. ---'';\nend\n']));
            pin = make_remote(state, 'matcutest', matcutest_files(''));
            lock_commits(state, struct('matcutest', pin));
            output = run_setup(state, struct('install_matcutest', true));
            testCase.verifyFalse(contains(output, 'user-managed'));
            testCase.verifySubstring(output, 'MatCUTEst managed runtime prepared');
            testCase.verifyTrue(on_path(fullfile(state.libraries, 'runtime', 'matcutest', 'mtools', 'src')));
        end

        function runtimeOnThePathIsUsedAndNotClaimed(testCase)
            testCase.assumeFalse(ismac);
            state = sandbox(testCase);
            tools = user_runtime(state, {'macup', 'secup', 'matcutest'});
            pin = make_remote(state, 'matcutest', matcutest_files(''));
            lock_commits(state, struct('matcutest', pin));
            output = run_setup(state, struct('install_matcutest', true));
            testCase.verifySubstring(output, 'existing user-managed MatCUTEst runtime');
            testCase.verifyTrue(is_registered(state, 'matcutest'));
            testCase.verifyTrue(on_path(tools));
            data = load([state.registry, '.setup-paths.mat']);
            testCase.verifyFalse(any(contains([data.records.owned_paths], 'user_runtime')));
        end

        function runtimeWithoutItsEntryPointsIsNotUsed(testCase, partial_runtime)
            testCase.assumeFalse(ismac);
            state = sandbox(testCase);
            user_runtime(state, {partial_runtime});
            pin = make_remote(state, 'matcutest', matcutest_files(''));
            lock_commits(state, struct('matcutest', pin));
            output = run_setup(state, struct('install_matcutest', true));
            testCase.verifyFalse(contains(output, 'user-managed'));
            testCase.verifySubstring(output, 'MatCUTEst managed runtime prepared');
        end

        function failedManagedRuntimeIsReported(testCase)
            testCase.assumeFalse(ismac);
            state = sandbox(testCase);
            pin = make_remote(state, 'matcutest', matcutest_files('error(''Fixture:Runtime'', ''no runtime here'');'));
            lock_commits(state, struct('matcutest', pin));
            output = run_setup(state, struct('install_matcutest', true));
            verify_core_configured(testCase, state);
            verify_reported(testCase, output, 'MatCUTEst', 'no runtime here');
            testCase.verifyFalse(is_registered(state, 'matcutest'));
        end

        function workingMatcutestAdapterIsRegisteredWithItsRuntime(testCase)
            testCase.assumeFalse(ismac);
            state = sandbox(testCase);
            pin = make_remote(state, 'matcutest', matcutest_files(''));
            lock_commits(state, struct('matcutest', pin));
            output = run_setup(state, struct('install_matcutest', true));
            verify_core_configured(testCase, state);
            testCase.verifyFalse(contains(output, 'is not set up'));
            testCase.verifyEqual(head(state, 'matcutest'), pin);
            testCase.verifyTrue(is_registered(state, 'matcutest'));
            testCase.verifyTrue(on_path(fullfile(state.libraries, 'matcutest')));
            testCase.verifyTrue(on_path(fullfile(state.libraries, 'runtime', 'matcutest', 'mtools', 'src')));
        end

        function failedSolarCloneIsReportedWithManualInstructions(testCase)
            state = sandbox(testCase);
            output = run_setup(state, struct('install_solar', true));
            verify_core_configured(testCase, state);
            verify_reported(testCase, output, 'SOLAR MATLAB adapter', 'Failed to clone');
            testCase.verifySubstring(output, 'manually clone or download');
            testCase.verifyFalse(is_registered(state, 'solar'));
            testCase.verifyFalse(on_path(fullfile(state.libraries, 'solar')));
        end

        function unusableLibraryRootIsReported(testCase)
            state = sandbox(testCase);
            occupied = fullfile(state.root, 'a_file');
            write_text(occupied, 'not a directory');
            output = run_setup(state, struct('install_solar', true, ...
                'problem_library_root', fullfile(occupied, 'libraries')));
            verify_core_configured(testCase, state);
            verify_reported(testCase, output, 'SOLAR MATLAB adapter', 'optional problem-library root');
            testCase.verifyFalse(is_registered(state, 'solar'));
        end

        function failedRegistrationLeavesThePathsAndTheOtherRegistrationAlone(testCase)
            state = sandbox(testCase);
            pin = make_remote(state, 'solar_matlab', solar_files());
            clone(state, 'solar_matlab', 'solar');
            lock_commits(state, struct('solar', pin));
            run_setup(state, struct('install_solar', true));
            % A registration of the same name at another root is not replaced.
            moved = fullfile(state.root, 'other_libraries');
            mkdir(moved);
            git(state.root, sprintf('clone -q "%s" "%s"', remote_path(state, 'solar_matlab'), fullfile(moved, 'solar')));
            output = run_setup(state, struct('install_solar', true, 'problem_library_root', moved));
            verify_core_configured(testCase, state);
            verify_reported(testCase, output, 'SOLAR MATLAB adapter', 'already registered');
            testCase.verifyFalse(on_path(fullfile(moved, 'solar')));
            testCase.verifyFalse(contains(fileread(fullfile(state.root, 'pathdef.m')), fullfile(moved, 'solar')));
            testCase.verifyTrue(is_registered(state, 'solar'));
            data = load(state.registry);
            testCase.verifyEqual(data.registrations.root, fullfile(state.libraries, 'solar'));
        end

        function libraryThatLeavesItsPinBehindIsWithdrawn(testCase)
            state = sandbox(testCase);
            pin = make_remote(state, 'solar_matlab', solar_files());
            lock_commits(state, struct('solar', pin));
            run_setup(state, struct('install_solar', true));
            testCase.assertTrue(is_registered(state, 'solar'));
            lock_commits(state, struct('solar', unknown_commit()));
            output = run_setup(state, struct('install_solar', true));
            verify_reported(testCase, output, 'SOLAR MATLAB adapter', 'Cannot check out the pinned commit');
            verify_core_configured(testCase, state);
            testCase.verifyFalse(is_registered(state, 'solar'));
            testCase.verifyFalse(on_path(fullfile(state.libraries, 'solar')));
            testCase.verifyFalse(contains(fileread(fullfile(state.root, 'pathdef.m')), fullfile(state.libraries, 'solar')));
            % The checkout stays for a later setup.
            testCase.verifyEqual(head(state, 'solar'), pin);
        end

        function relativeLibraryRootIsResolvedSoThatItsPathsCanBeWithdrawn(testCase)
            state = sandbox(testCase);
            pin = make_remote(state, 'solar_matlab', solar_files());
            lock_commits(state, struct('solar', pin));
            relative = struct('install_solar', true, 'problem_library_root', fullfile('..', 'libraries'));
            run_setup(state, relative);
            solar = fullfile(state.libraries, 'solar');
            testCase.assertTrue(is_registered(state, 'solar'));
            testCase.verifyTrue(on_path(solar));
            lock_commits(state, struct('solar', unknown_commit()));
            run_setup(state, relative);
            testCase.verifyFalse(is_registered(state, 'solar'));
            testCase.verifyFalse(on_path(solar));
            testCase.verifyFalse(contains(fileread(fullfile(state.root, 'pathdef.m')), solar));
        end

        function linkedRootIsRecordedUnderItsResolvedName(testCase, linked_root)
            state = throughALinkedRoot(testCase, linked_root, sandbox(testCase));
            run_setup(state, state.options);
            solar = fullfile(state.libraries, 'solar');
            testCase.verifyTrue(on_path(solar));
            data = load([state.registry, '.setup-paths.mat']);
            testCase.verifyTrue(any(strcmp({data.records.owner}, solar)));
            testCase.verifyFalse(any(strcmp({data.records.owner}, fullfile(state.link, 'solar'))));
            lock_commits(state, struct('solar', unknown_commit()));
            run_setup(state, state.options);
            testCase.verifyFalse(is_registered(state, 'solar'));
            verify_not_on_path(testCase, state);
        end

        function recordThroughALinkedRootIsWithdrawn(testCase, linked_root)
            state = recordedThroughALinkedRoot(testCase, linked_root);
            lock_commits(state, struct('solar', unknown_commit()));
            output = run_setup(state, state.options);
            verify_reported(testCase, output, 'SOLAR MATLAB adapter', 'Cannot check out the pinned commit');
            testCase.verifyFalse(is_registered(state, 'solar'));
            verify_not_on_path(testCase, state);
            data = load([state.registry, '.setup-paths.mat']);
            record = data.records(strcmp({data.records.owner}, fullfile(state.link, 'solar')));
            testCase.verifyEmpty(record.owned_paths);
        end

        function unregisteredRecordThroughALinkedRootIsWithdrawn(testCase, linked_root)
            state = recordedThroughALinkedRoot(testCase, linked_root);
            delete(state.registry);
            lock_commits(state, struct('solar', unknown_commit()));
            output = run_setup(state, state.options);
            verify_reported(testCase, output, 'SOLAR MATLAB adapter', 'Cannot check out the pinned commit');
            verify_not_on_path(testCase, state);
        end

        function startupLineThroughALinkedRootIsWithdrawn(testCase, linked_root)
            state = sandbox(testCase);
            setenv('OPTIPROFILER_MATLAB_PROBLEM_LIBRARY_PATHDEF', fullfile(state.root, 'no_directory', 'pathdef.m'));
            state = recordedThroughALinkedRoot(testCase, linked_root, state);
            startup = fullfile(state.root, 'startup.m');
            testCase.assertSubstring(fileread(startup), fullfile(state.link, 'solar'));
            lock_commits(state, struct('solar', unknown_commit()));
            run_setup(state, state.options);
            testCase.verifyFalse(contains(fileread(startup), fullfile(state.link, 'solar')));
            testCase.verifyFalse(on_path(fullfile(state.libraries, 'solar')));
        end

        function libraryKeptThroughALinkedRootCanBeUnregistered(testCase, linked_root)
            state = recordedThroughALinkedRoot(testCase, linked_root);
            output = run_setup(state, state.options);
            testCase.verifyFalse(contains(output, 'is not set up'));
            testCase.assertTrue(is_registered(state, 'solar'));
            unregisterProblemLibrary('solar');
            testCase.verifyFalse(is_registered(state, 'solar'));
            verify_not_on_path(testCase, state);
        end

        function setupRecoversOnceTheLibraryIsReachable(testCase)
            state = sandbox(testCase);
            run_setup(state, struct('install_solar', true));
            testCase.verifyFalse(is_registered(state, 'solar'));
            pin = make_remote(state, 'solar_matlab', solar_files());
            lock_commits(state, struct('solar', pin));
            output = run_setup(state, struct('install_solar', true));
            testCase.verifyFalse(contains(output, 'is not set up'));
            testCase.verifyEqual(head(state, 'solar'), pin);
            testCase.verifyTrue(is_registered(state, 'solar'));
        end

        function localChangesAreLeftAlone(testCase)
            state = sandbox(testCase);
            base = make_remote(state, 'solar_matlab', solar_files());
            clone(state, 'solar_matlab', 'solar');
            pin = branch_commit(state, 'solar_matlab');
            edited = fullfile(state.libraries, 'solar', 'solar_load.m');
            write_text(edited, sprintf('function p = solar_load(~)\n%% my local edit\np = [];\nend\n'));
            lock_commits(state, struct('solar', pin));
            output = run_setup(state, struct('install_solar', true));
            testCase.verifySubstring(output, 'has local changes');
            testCase.verifyEqual(head(state, 'solar'), base);
            testCase.verifySubstring(fileread(edited), 'my local edit');
            testCase.verifyTrue(is_registered(state, 'solar'));
        end

        function directoryThatIsNotARepositoryIsUsedAsIs(testCase)
            state = sandbox(testCase);
            root = fullfile(state.libraries, 'solar');
            mkdir(root);
            files = solar_files();
            for k = 1:size(files, 1)
                write_text(fullfile(root, files{k, 1}), files{k, 2});
            end
            lock_commits(state, struct('solar', unknown_commit()));
            output = run_setup(state, struct('install_solar', true));
            testCase.verifySubstring(output, 'not a git repository');
            testCase.verifyTrue(is_registered(state, 'solar'));
        end
    end
end


function state = throughALinkedRoot(testCase, form, state)
% The library at its pin, with options that reach its root through a symbolic link.
    pin = make_remote(state, 'solar_matlab', solar_files());
    lock_commits(state, struct('solar', pin));
    if strcmp(form, 'itself')
        target = state.libraries;
        link = fullfile(state.root, 'linked_libraries');
        state.link = link;
    else
        target = fullfile(state.root, 'real');
        state.libraries = fullfile(target, 'libraries');
        mkdir(state.libraries);
        link = fullfile(state.root, 'link');
        state.link = fullfile(link, 'libraries');
    end
    [status, output] = system(sprintf('ln -s "%s" "%s"', target, link));
    testCase.assertEqual(status, 0, output);
    state.options = struct('install_solar', true, 'problem_library_root', state.link);
end

function state = recordedThroughALinkedRoot(testCase, form, state)
% The state that a setup leaves that does not resolve problem_library_root: the
% MATLAB path and the registry hold the resolved name, and the ownership record
% holds the name as it was given.
    if nargin < 3
        state = sandbox(testCase);
    end
    state = throughALinkedRoot(testCase, form, state);
    setup_file = fullfile(state.package, 'setup.m');
    current = fileread(setup_file);
    unresolved = regexprep(current, 'root = resolve\w*\(root\);', '');
    testCase.assertNotEqual(unresolved, current);
    write_text(setup_file, unresolved);
    run_setup(state, state.options);
    write_text(setup_file, current);
    data = load([state.registry, '.setup-paths.mat']);
    testCase.assertTrue(any(strcmp({data.records.owner}, fullfile(state.link, 'solar'))));
    testCase.assertTrue(on_path(fullfile(state.libraries, 'solar')));
end

function verify_not_on_path(testCase, state)
    testCase.verifyFalse(on_path(fullfile(state.libraries, 'solar')));
    testCase.verifyFalse(contains(fileread(fullfile(state.root, 'pathdef.m')), 'solar'));
end

function state = sandbox(testCase)
    state.root = tempname;
    mkdir(state.root);
    [~, attributes] = fileattrib(state.root);
    state.root = attributes.Name;
    state.home = fullfile(state.root, 'home');
    state.remotes = fullfile(state.root, 'remotes');
    state.libraries = fullfile(state.root, 'libraries');
    state.package = fullfile(state.root, 'package');
    state.registry = fullfile(state.root, 'registry.mat');
    mkdir(state.home);
    mkdir(fullfile(state.remotes, 'optiprofiler'));
    mkdir(state.libraries);

    old_path = path;
    old_dir = pwd;
    variables = {'HOME', 'XDG_CONFIG_HOME', 'GIT_CONFIG_NOSYSTEM', ...
        'OPTIPROFILER_MATLAB_PROBLEM_LIBRARY_REGISTRY', ...
        'OPTIPROFILER_MATLAB_PROBLEM_LIBRARY_PATHDEF', ...
        'OPTIPROFILER_MATLAB_PROBLEM_LIBRARY_STARTUP'};
    values = cellfun(@getenv, variables, 'UniformOutput', false);
    testCase.addTeardown(@() restore(state.root, old_path, old_dir, variables, values));
    setenv('HOME', state.home);
    setenv('XDG_CONFIG_HOME', fullfile(state.home, '.config'));
    setenv('GIT_CONFIG_NOSYSTEM', '1');
    setenv('OPTIPROFILER_MATLAB_PROBLEM_LIBRARY_REGISTRY', state.registry);
    setenv('OPTIPROFILER_MATLAB_PROBLEM_LIBRARY_PATHDEF', fullfile(state.root, 'pathdef.m'));
    setenv('OPTIPROFILER_MATLAB_PROBLEM_LIBRARY_STARTUP', fullfile(state.root, 'startup.m'));
    write_text(fullfile(state.home, '.gitconfig'), sprintf([ ...
        '[user]\n\tname = OptiProfiler Test\n\temail = test@example.org\n' ...
        '[url "file://%s/"]\n\tinsteadOf = https://github.com/\n' ...
        '[protocol "file"]\n\tallow = always\n' ...
        '[uploadpack]\n\tallowAnySHA1InWant = true\n' ...
        '[advice]\n\tdetachedHead = false\n'], state.remotes));

    test_dir = fileparts(mfilename('fullpath'));
    [~, attributes] = fileattrib(fullfile(test_dir, '..', '..', '..', '..'));
    repository = attributes.Name;
    src = fullfile(state.package, 'matlab', 'optiprofiler', 'src');
    mkdir(fileparts(src));
    copyfile(fullfile(repository, 'matlab', 'optiprofiler', 'src'), src);
    copyfile(fullfile(repository, 'setup.m'), state.package);
    copyfile(fullfile(repository, 'matlab_problem_libraries.lock'), state.package);
    s2 = fullfile(state.package, 'matlab', 'optiprofiler', 'problem_libs', 's2mpj');
    mkdir(fullfile(s2, 'src', 'matlab_problems'));
    write_text(fullfile(s2, 'README.md'), 'fixture');
end

function restore(root, old_path, old_dir, variables, values)
    path(old_path);
    cd(old_dir);
    clear setup
    for k = 1:numel(variables)
        setenv(variables{k}, values{k});
    end
    if isfolder(root)
        rmdir(root, 's');
    end
end

function output = run_setup(state, options)
    defaults = struct('install_matcutest', false, 'install_solar', false, ...
        'problem_library_root', state.libraries);
    for field = reshape(fieldnames(options), 1, [])
        defaults.(field{1}) = options.(field{1});
    end
    cd(state.package);
    clear setup
    output = evalc('setup(defaults)');
end

function [output, failure] = run_failing_setup(state, options)
% evalc keeps nothing of a call that raises.
    defaults = struct('install_matcutest', false, 'install_solar', false, ...
        'problem_library_root', state.libraries);
    for field = reshape(fieldnames(options), 1, [])
        defaults.(field{1}) = options.(field{1});
    end
    cd(state.package);
    clear setup
    failure = [];
    output = evalc('try, setup(defaults); catch failure, end');
    assert(~isempty(failure), 'setup did not raise.');
end

function files = solar_files()
    files = { ...
        'solar_select.m', sprintf('function names = solar_select(~)\nnames = {};\nend\n'); ...
        'solar_load.m', sprintf('function p = solar_load(~)\np = [];\nend\n'); ...
        'solar_collect_info.m', sprintf('function info = solar_collect_info()\ninfo = {};\nend\n')};
end

function files = matcutest_files(extra)
    body = sprintf(['function r = matcutest_setup(parent)\n%s\n' ...
        'if ~isfolder(parent), mkdir(parent); end\n' ...
        '[~, a] = fileattrib(parent);\n' ...
        'root = fullfile(a.Name, ''matcutest'');\n' ...
        'tools = fullfile(root, ''mtools'', ''src'');\n' ...
        'if ~isfolder(tools), mkdir(tools); end\n' ...
        'r = struct(''schema_version'', 1, ''provider'', ''matcutest'', ' ...
        '''runtime_root'', root, ''runtime_paths'', {{tools}});\nend\n'], extra);
    files = { ...
        'matcutest_select.m', sprintf('function names = matcutest_select(~)\nnames = {};\nend\n'); ...
        'matcutest_load.m', sprintf('function p = matcutest_load(~)\np = [];\nend\n'); ...
        'matcutest_setup.m', body};
end

function tools = user_runtime(state, names)
% A MatCUTEst of the user on the MATLAB path, with the functions NAMES.
    tools = fullfile(state.root, 'user_runtime', 'mtools', 'src');
    mkdir(tools);
    for name = names
        write_text(fullfile(tools, [name{1}, '.m']), sprintf('function %s\nend\n', name{1}));
    end
    addpath(tools);
end

function verify_core_configured(testCase, state)
    src = fullfile(state.package, 'matlab', 'optiprofiler', 'src');
    s2 = fullfile(state.package, 'matlab', 'optiprofiler', 'problem_libs', 's2mpj');
    pathdef = fileread(fullfile(state.root, 'pathdef.m'));
    for directory = {src, s2}
        testCase.verifyTrue(on_path(directory{1}), directory{1});
        testCase.verifySubstring(pathdef, directory{1});
    end
end

function verify_reported(testCase, output, library, reason)
    % MATLAB wraps the lines of a warning.
    output = regexprep(output, '\s+', ' ');
    testCase.verifySubstring(output, sprintf('The optional problem library %s is not set up', library));
    testCase.verifySubstring(output, reason);
    testCase.verifySubstring(output, 'The package is ready to use');
end

function path_name = remote_path(state, name)
    path_name = fullfile(state.remotes, 'optiprofiler', [name, '.git']);
end

function sha = make_remote(state, name, files)
    remote = remote_path(state, name);
    mkdir(remote);
    git(remote, 'init -q');
    git(remote, 'symbolic-ref HEAD refs/heads/main');
    for k = 1:size(files, 1)
        target = fullfile(remote, files{k, 1});
        if ~isfolder(fileparts(target)), mkdir(fileparts(target)); end
        write_text(target, files{k, 2});
    end
    git(remote, 'add -A');
    git(remote, 'commit -q -m base');
    sha = rev_parse(remote);
end

function sha = branch_commit(state, name)
% A commit that only a side branch reaches, as a pin on an unmerged branch is.
    remote = remote_path(state, name);
    git(remote, 'checkout -q -b pinned');
    write_text(fullfile(remote, 'pinned.txt'), 'pinned');
    git(remote, 'add -A');
    git(remote, 'commit -q -m pinned');
    sha = rev_parse(remote);
    git(remote, 'checkout -q main');
end

function clone(state, name, directory)
    target = fullfile(state.libraries, directory);
    git(state.root, sprintf('clone -q "https://github.com/optiprofiler/%s.git" "%s"', name, target));
end

function lock_commits(state, commits)
    lock_file = fullfile(state.package, 'matlab_problem_libraries.lock');
    text = fileread(lock_file);
    decoded = jsondecode(text);
    for name = reshape(fieldnames(commits), 1, [])
        for k = 1:numel(decoded.libraries)
            entry = decoded.libraries;
            if iscell(entry), entry = entry{k}; else, entry = entry(k); end
            if strcmp(entry.name, name{1})
                text = strrep(text, entry.commit, commits.(name{1}));
            end
        end
    end
    write_text(lock_file, text);
end

function sha = unknown_commit()
    sha = '0123456789abcdef0123456789abcdef01234567';
end

function sha = head(state, directory)
    sha = rev_parse(fullfile(state.libraries, directory));
end

function sha = rev_parse(directory)
    sha = strtrim(git(directory, 'rev-parse HEAD'));
end

function output = git(directory, arguments)
    [status, output] = system(sprintf('git -C "%s" %s', directory, arguments));
    assert(status == 0, 'git %s failed in %s: %s', arguments, directory, output);
end

function tf = is_registered(state, name)
    tf = false;
    if isfile(state.registry)
        data = load(state.registry);
        tf = any(strcmp({data.registrations.name}, name));
    end
end

function tf = on_path(directory)
    tf = any(strcmp(strsplit(path, pathsep), directory));
end

function write_text(filename, text)
    fid = fopen(filename, 'w');
    assert(fid >= 0, 'Cannot write %s.', filename);
    cleanup = onCleanup(@() fclose(fid));
    fwrite(fid, text, 'char');
end
