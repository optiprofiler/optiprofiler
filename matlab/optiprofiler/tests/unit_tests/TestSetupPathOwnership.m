classdef TestSetupPathOwnership < matlab.unittest.TestCase
    methods (Test)
        function testBorrowedPathsAndRepeatedSetup(testCase)
            state = isolate(testCase);
            borrowed = fullfile(state.root, 'borrowed');
            owned = fullfile(state.root, 'owned');
            mkdir(borrowed); mkdir(owned);
            addpath(borrowed);
            setupPathOwnership('add', state.root, fullfile(state.root, 'provider'), {borrowed, owned});
            setupPathOwnership('add', state.root, fullfile(state.root, 'provider'), {borrowed, owned});
            data = load([state.registry, '.setup-paths.mat']);
            testCase.verifyEqual(data.records.owned_paths, {owned});
            testCase.verifyTrue(setupPathOwnership('remove-owner', fullfile(state.root, 'provider')));
            testCase.verifyTrue(on_path(borrowed));
            testCase.verifyFalse(on_path(owned));
            testCase.verifyTrue(setupPathOwnership('remove-context', state.root));
            testCase.verifyTrue(on_path(borrowed));
            testCase.verifyTrue(isfolder(owned));
        end

        function testFallbackPreservesUnrelatedStartupBytes(testCase)
            state = isolate(testCase);
            setenv('OPTIPROFILER_MATLAB_PROBLEM_LIBRARY_PATHDEF', ...
                fullfile(state.root, 'missing', 'pathdef.m'));
            before = unicode2native(sprintf('%% user mentions %s and %% optiprofiler\r\n%% keep trailing bytes', state.root), 'UTF-8');
            write_bytes(state.startup, before);
            owned = fullfile(state.root, 'owned');
            mkdir(owned);
            setupPathOwnership('add', state.root, fullfile(state.root, 'provider'), {owned});
            installed = read_bytes(state.startup);
            setupPathOwnership('add', state.root, fullfile(state.root, 'provider'), {owned});
            testCase.verifyEqual(read_bytes(state.startup), installed);
            suffix = uint8(sprintf('\r\n%% unrelated later user text'));
            write_bytes(state.startup, [installed, suffix]);
            setupPathOwnership('remove-context', state.root);
            testCase.verifyEqual(read_bytes(state.startup), [before, suffix]);
        end

        function testUnknownOwnershipDoesNotGuess(testCase)
            state = isolate(testCase);
            addpath(state.root);
            before = uint8('% optiprofiler user comment');
            write_bytes(state.startup, before);
            testCase.verifyFalse(setupPathOwnership('remove-context', state.root));
            testCase.verifyTrue(on_path(state.root));
            testCase.verifyEqual(read_bytes(state.startup), before);
        end

        function testPersistenceTargetChangeFails(testCase)
            state = isolate(testCase);
            owned = fullfile(state.root, 'owned');
            mkdir(owned);
            setupPathOwnership('add', state.root, fullfile(state.root, 'provider'), {owned});
            setenv('OPTIPROFILER_MATLAB_PROBLEM_LIBRARY_STARTUP', ...
                fullfile(state.root, 'other-startup.m'));
            testCase.verifyError(@() setupPathOwnership('add', ...
                state.root, fullfile(state.root, 'provider'), {owned}), 'OptiProfiler:PathOwnership');
        end

        function testCorruptLedgerFailsBeforeCleanup(testCase)
            state = isolate(testCase);
            addpath(state.root);
            unrelated = true; %#ok<NASGU>
            save([state.registry, '.setup-paths.mat'], 'unrelated');
            testCase.verifyError(@() setupPathOwnership('remove-context', ...
                state.root), 'OptiProfiler:PathOwnership');
            testCase.verifyTrue(on_path(state.root));
        end

        function testRelativeOwnershipLocationIsRejected(testCase)
            state = isolate(testCase);
            setenv('OPTIPROFILER_MATLAB_PROBLEM_LIBRARY_REGISTRY', 'relative.mat');
            testCase.verifyError(@() setupPathOwnership('add', state.root, ...
                fullfile(state.root, 'provider'), {state.root}), 'OptiProfiler:PathOwnership');
        end

        function testWindowsRelativePathsRejected(testCase)
            % This is an actual Windows test, not a simulated ispc branch.
            testCase.assumeTrue(ispc);
            state = isolate(testCase);
            [~, unique] = fileparts(state.root);
            invalid = {['\', unique, '\registry.mat'], ['/', unique, '/registry.mat'], ...
                ['C:', unique, '\registry.mat'], [unique, '\registry.mat']};
            for k = 1:numel(invalid)
                setenv('OPTIPROFILER_MATLAB_PROBLEM_LIBRARY_REGISTRY', invalid{k});
                % Read-only removal detects bad registry locations without
                % creating a file at a root-relative destination on failure.
                testCase.verifyError(@() setupPathOwnership('remove-context', ...
                    state.root), 'OptiProfiler:PathOwnership');
            end
            setenv('OPTIPROFILER_MATLAB_PROBLEM_LIBRARY_REGISTRY', state.registry);
            for k = 1:numel(invalid)
                testCase.verifyError(@() setupPathOwnership('add', ...
                    invalid{k}, state.root, {state.root}), 'OptiProfiler:PathOwnership');
            end
            testCase.verifyFalse(isfile([state.registry, '.setup-paths.mat']));
        end

        function testWindowsAbsolutePathForms(testCase)
            testCase.assumeTrue(ispc);
            state = isolate(testCase);
            forms = {'Z:\project', 'Z:/project', '\\example.invalid\share\project'};
            owned = fullfile(state.root, 'owned'); mkdir(owned);
            for k = 1:numel(forms)
                % Context/owner are labels only: all I/O stays in state.root.
                % In particular, the test never accesses the dummy UNC share.
                testCase.verifyTrue(all(setupPathOwnership('add', ...
                    forms{k}, forms{k}, {owned})));
                testCase.verifyTrue(setupPathOwnership('remove-context', forms{k}));
                testCase.verifyFalse(on_path(owned));
            end
        end

        function testPosixRejectsWindowsDrivePaths(testCase)
            testCase.assumeTrue(isunix);
            state = isolate(testCase);
            for value = {'C:\registry.mat', 'C:/registry.mat', 'C:registry.mat'}
                setenv('OPTIPROFILER_MATLAB_PROBLEM_LIBRARY_REGISTRY', value{1});
                testCase.verifyError(@() setupPathOwnership('remove-context', ...
                    state.root), 'OptiProfiler:PathOwnership');
            end
        end

        function testPersistedUserPathSurvivesDifferentSessionPath(testCase)
            state = isolate(testCase);
            user_dir = fullfile(state.root, 'user');
            owned = fullfile(state.root, 'owned');
            mkdir(user_dir); mkdir(owned);
            addpath(user_dir);
            setupPathOwnership('add', state.root, fullfile(state.root, 'provider'), {owned});
            rmpath(user_dir);
            setupPathOwnership('remove-context', state.root);
            text = fileread(fullfile(state.root, 'pathdef.m'));
            testCase.verifySubstring(text, user_dir);
            testCase.verifyFalse(contains(text, owned));
            testCase.verifyFalse(on_path(user_dir));
        end

        function testFailedStartupAppendCanRetry(testCase)
            state = isolate(testCase);
            setenv('OPTIPROFILER_MATLAB_PROBLEM_LIBRARY_PATHDEF', fullfile(state.root, 'missing', 'pathdef.m'));
            startup_dir = fullfile(state.root, 'not-yet-created');
            startup = fullfile(startup_dir, 'startup.m');
            setenv('OPTIPROFILER_MATLAB_PROBLEM_LIBRARY_STARTUP', startup);
            owned = fullfile(state.root, 'owned');
            mkdir(owned);
            callback = @() setupPathOwnership('add', state.root, fullfile(state.root, 'provider'), {owned});
            testCase.verifyError(callback, 'OptiProfiler:PathPersistence');
            mkdir(startup_dir);
            callback();
            testCase.assertTrue(isfile(startup), 'Retry must perform the previously failed write.');
            testCase.verifySubstring(fileread(startup), owned);
            setupPathOwnership('remove-context', state.root);
            testCase.verifyEmpty(read_bytes(startup));
        end

        function testExistingUserStartupEntryIsNotClaimed(testCase)
            state = isolate(testCase);
            setenv('OPTIPROFILER_MATLAB_PROBLEM_LIBRARY_PATHDEF', fullfile(state.root, 'missing', 'pathdef.m'));
            owned = fullfile(state.root, 'owned');
            mkdir(owned);
            original = uint8(sprintf('\naddpath(''%s''); %% optiprofiler setup-owned', owned));
            write_bytes(state.startup, original);
            setupPathOwnership('add', state.root, fullfile(state.root, 'provider'), {owned});
            setupPathOwnership('remove-context', state.root);
            testCase.verifyEqual(read_bytes(state.startup), original);
        end

        function testSidecarWriteFailureDoesNotMutatePaths(testCase)
            state = isolate(testCase);
            mkdir([state.registry, '.setup-paths.mat']);
            before = path;
            testCase.verifyError(@() setupPathOwnership('add', state.root, ...
                fullfile(state.root, 'provider'), {state.root}), 'OptiProfiler:PathPersistence');
            testCase.verifyEqual(path, before);
        end

        function testRemovalWriteFailurePreservesBytesAndCanRetry(testCase)
            state = isolate(testCase);
            setenv('OPTIPROFILER_MATLAB_PROBLEM_LIBRARY_PATHDEF', fullfile(state.root, 'missing', 'pathdef.m'));
            owned = fullfile(state.root, 'owned');
            mkdir(owned);
            setupPathOwnership('add', state.root, fullfile(state.root, 'provider'), {owned});
            original = read_bytes(state.startup);
            assert(fileattrib(state.startup, '-w'));
            testCase.verifyError(@() setupPathOwnership('remove-context', state.root), 'OptiProfiler:PathPersistence');
            testCase.verifyEqual(read_bytes(state.startup), original);
            records = load([state.registry, '.setup-paths.mat']);
            testCase.verifyEqual(records.records.owned_paths, {owned});
            assert(fileattrib(state.startup, '+w'));
            setupPathOwnership('remove-context', state.root);
            testCase.verifyEmpty(read_bytes(state.startup));
        end

        function testPersistedPathIsBorrowedEvenWhenNotCurrentlyLoaded(testCase)
            state = isolate(testCase);
            user_dir = fullfile(state.root, 'user');
            mkdir(user_dir);
            addpath(user_dir);
            assert(savepath(fullfile(state.root, 'pathdef.m')) == 0);
            rmpath(user_dir);
            setupPathOwnership('add', state.root, fullfile(state.root, 'provider'), {user_dir});
            records = load([state.registry, '.setup-paths.mat']);
            testCase.verifyEmpty(records.records.owned_paths);
            setupPathOwnership('remove-context', state.root);
            testCase.verifyTrue(on_path(user_dir));
            testCase.verifySubstring(fileread(fullfile(state.root, 'pathdef.m')), user_dir);
        end

        function testInvalidRuntimeReceiptCannotAddOutsidePath(testCase)
            testCase.assumeTrue(isunix && ~ismac);
            state = isolate(testCase);
            adapter = fullfile(state.root, 'adapter');
            runtime_parent = fullfile(state.root, 'runtime');
            runtime_root = fullfile(runtime_parent, 'matcutest');
            mkdir(adapter); mkdir(fullfile(runtime_root, 'mtools', 'src'));
            source = sprintf(['function r = matcutest_setup(~)\n' ...
                'r = struct(''schema_version'',1,''provider'',''matcutest'',''runtime_root'',''%s'',''runtime_paths'',{{''%s''}});\nend\n'], ...
                runtime_root, state.root);
            write_bytes(fullfile(adapter, 'matcutest_setup.m'), uint8(source));
            before = path;
            testCase.verifyError(@() prepareMatcutestRuntime(adapter, runtime_parent), ...
                'OptiProfiler:InvalidMatcutestReceipt');
            testCase.verifyEqual(path, before);
        end

        function testStartupPosixModesPreserved(testCase)
            testCase.assumeTrue(isunix);
            for mode = {'600', '640', '644'}
                state = isolate(testCase);
                setenv('OPTIPROFILER_MATLAB_PROBLEM_LIBRARY_PATHDEF', fullfile(state.root, 'missing', 'pathdef.m'));
                original = uint8(sprintf('%% private startup\r\n'));
                write_bytes(state.startup, original);
                set_mode(state.startup, mode{1});
                owned = fullfile(state.root, 'owned'); mkdir(owned);
                for repeat = 1:2
                    setupPathOwnership('add', state.root, state.root, {owned});
                    testCase.verifyEqual(get_mode(state.startup), mode{1});
                end
                setupPathOwnership('remove-context', state.root);
                testCase.verifyEqual(get_mode(state.startup), mode{1});
                testCase.verifyEqual(read_bytes(state.startup), original);
            end
        end

        function testPathdefPosixModesPreserved(testCase)
            testCase.assumeTrue(isunix);
            for mode = {'600', '640', '644'}
                state = isolate(testCase);
                filename = fullfile(state.root, 'pathdef.m');
                assert(savepath(filename) == 0);
                set_mode(filename, mode{1});
                owned = fullfile(state.root, 'owned'); mkdir(owned);
                for repeat = 1:2
                    setupPathOwnership('add', state.root, state.root, {owned});
                    testCase.verifyEqual(get_mode(filename), mode{1});
                end
                setupPathOwnership('remove-context', state.root);
                testCase.verifyEqual(get_mode(filename), mode{1});
                testCase.verifyFalse(contains(fileread(filename), owned));
            end
        end

        function testNewPersistenceFilesArePrivate(testCase)
            testCase.assumeTrue(isunix);
            for fallback = [false, true]
                state = isolate(testCase);
                filename = fullfile(state.root, 'pathdef.m');
                if fallback
                    setenv('OPTIPROFILER_MATLAB_PROBLEM_LIBRARY_PATHDEF', fullfile(state.root, 'missing', 'pathdef.m'));
                    filename = state.startup;
                end
                owned = fullfile(state.root, 'owned'); mkdir(owned);
                setupPathOwnership('add', state.root, state.root, {owned});
                testCase.verifyEqual(get_mode(filename), '600');
                setupPathOwnership('remove-context', state.root);
                testCase.verifyEqual(get_mode(filename), '600');
            end
        end

        function testPrivatePathdefStagingFailureAndRetry(testCase)
            testCase.assumeTrue(isunix);
            state = isolate(testCase);
            % Persist the fixture directory before creating the shadow: setup
            % deliberately merges persisted entries ahead of session entries.
            stub = fullfile(state.root, 'stub'); mkdir(stub); addpath(stub, '-begin');
            filename = fullfile(state.root, 'pathdef.m');
            assert(savepath(filename) == 0);
            set_mode(filename, '640');
            original = read_bytes(filename);
            owned = fullfile(state.root, 'owned'); mkdir(owned);
            source = sprintf(['function result = savepath(filename)\n' ...
                'assertPrivate(fileparts(filename));\n' ...
                'fid=fopen(filename,''w''); c=onCleanup(@() fclose(fid));\n' ...
                'fprintf(fid,''%%s'',''partial private path data'');\n' ...
                'assertPrivate(fileparts(filename));\n' ...
                'error(''OptiProfilerTest:PartialPathdef'',''Injected partial write.'');\n' ...
                'result=1;\nend\n%s'], private_assertion());
            write_bytes(fullfile(stub, 'savepath.m'), uint8(source));
            addpath(stub, '-begin');
            callback = @() setupPathOwnership('add', state.root, state.root, {owned});
            testCase.verifyError(callback, 'OptiProfilerTest:PartialPathdef');
            testCase.verifyEqual(read_bytes(filename), original);
            testCase.verifyEqual(get_mode(filename), '640');
            testCase.verifyEmpty(dir(fullfile(state.root, '.optiprofiler-setup-*')));
            delete(fullfile(stub, 'savepath.m')); rmpath(stub); clear savepath
            callback();
            testCase.verifyEqual(get_mode(filename), '640');
            setupPathOwnership('remove-context', state.root);
            testCase.verifyEqual(get_mode(filename), '640');
        end

        function testPrivateStartupStagingFailureAndRetry(testCase)
            testCase.assumeTrue(isunix);
            state = isolate(testCase);
            setenv('OPTIPROFILER_MATLAB_PROBLEM_LIBRARY_PATHDEF', fullfile(state.root, 'missing', 'pathdef.m'));
            original = uint8('% private startup');
            write_bytes(state.startup, original); set_mode(state.startup, '600');
            owned = fullfile(state.root, 'owned'); mkdir(owned);
            stub = fullfile(state.root, 'stub'); mkdir(stub);
            source = sprintf(['function result = fwrite(fid,data,varargin)\n' ...
                'filename=fopen(fid);\n' ...
                'if contains(filename,''.optiprofiler-setup-'')\n' ...
                'assertPrivate(fileparts(filename));\n' ...
                'builtin(''fwrite'',fid,data(1),varargin{:});\n' ...
                'assertPrivate(fileparts(filename));\n' ...
                'error(''OptiProfilerTest:PartialStartup'',''Injected partial write.'');\nend\n' ...
                'result=builtin(''fwrite'',fid,data,varargin{:});\nend\n%s'], private_assertion());
            write_bytes(fullfile(stub, 'fwrite.m'), uint8(source));
            addpath(stub, '-begin');
            callback = @() setupPathOwnership('add', state.root, state.root, {owned});
            testCase.verifyError(callback, 'OptiProfilerTest:PartialStartup');
            testCase.verifyEqual(read_bytes(state.startup), original);
            testCase.verifyEqual(get_mode(state.startup), '600');
            testCase.verifyEmpty(dir(fullfile(state.root, '.optiprofiler-setup-*')));
            data = load([state.registry, '.setup-paths.mat']);
            testCase.verifyEmpty(data.records.startup_additions);
            rmpath(stub); clear fwrite
            callback();
            testCase.verifyEqual(get_mode(state.startup), '600');
            setupPathOwnership('remove-context', state.root);
            testCase.verifyEqual(get_mode(state.startup), '600');
            testCase.verifyEqual(read_bytes(state.startup), original);
        end

        function testReadOnlyPathdefRemovalPreservesModeAndCanRetry(testCase)
            testCase.assumeTrue(isunix);
            state = isolate(testCase);
            filename = fullfile(state.root, 'pathdef.m');
            owned = fullfile(state.root, 'owned'); mkdir(owned);
            setupPathOwnership('add', state.root, state.root, {owned});
            original = read_bytes(filename); set_mode(filename, '400');
            testCase.verifyError(@() setupPathOwnership('remove-context', state.root), 'OptiProfiler:PathPersistence');
            testCase.verifyEqual(read_bytes(filename), original);
            testCase.verifyEqual(get_mode(filename), '400');
            data = load([state.registry, '.setup-paths.mat']);
            testCase.verifyEqual(data.records.owned_paths, {owned});
            set_mode(filename, '600');
            setupPathOwnership('remove-context', state.root);
            testCase.verifyEqual(get_mode(filename), '600');
        end

        function testPublicationFailureRetainsPrivateCompleteCandidate(testCase)
            testCase.assumeTrue(isunix);
            state = isolate(testCase);
            stub = fullfile(state.root, 'stub'); mkdir(stub); addpath(stub, '-begin');
            filename = fullfile(state.root, 'pathdef.m');
            assert(savepath(filename) == 0); set_mode(filename, '640');
            original = read_bytes(filename);
            backup = fullfile(state.root, 'original.m');
            owned = fullfile(state.root, 'owned'); mkdir(owned);
            source = sprintf(['function result = savepath(filename)\n' ...
                'assertPrivate(fileparts(filename));\n' ...
                'fid=fopen(filename,''w''); fprintf(fid,''%%s'',''completed private candidate''); fclose(fid);\n' ...
                'movefile(''%s'',''%s''); mkdir(''%s''); result=0;\nend\n%s'], ...
                strrep(filename, '''', ''''''), strrep(backup, '''', ''''''), ...
                strrep(filename, '''', ''''''), private_assertion());
            write_bytes(fullfile(stub, 'savepath.m'), uint8(source));
            addpath(stub, '-begin');
            callback = @() setupPathOwnership('add', state.root, state.root, {owned});
            testCase.verifyError(callback, 'OptiProfiler:PathPersistence');
            entries = dir(fullfile(state.root, '.optiprofiler-setup-*'));
            testCase.assertNumElements(entries, 1);
            staging = fullfile(state.root, entries(1).name);
            candidate = fullfile(staging, 'candidate.m');
            testCase.verifyEqual(get_mode(staging), '700');
            testCase.verifyEqual(get_mode(candidate), '640');
            testCase.verifyEqual(read_bytes(candidate), uint8('completed private candidate'));
            testCase.verifyEqual(read_bytes(backup), original);
            testCase.verifyEqual(get_mode(backup), '640');
            delete(fullfile(stub, 'savepath.m')); rmpath(stub); clear savepath
            rmdir(filename); movefile(backup, filename);
            callback();
            testCase.verifyEqual(get_mode(filename), '640');
            setupPathOwnership('remove-context', state.root);
        end

        function testSidecarPosixModesPreserved(testCase)
            testCase.assumeTrue(isunix);
            for mode = {'new', '600', '640', '644'}
                state = isolate(testCase);
                owned = fullfile(state.root, 'owned'); mkdir(owned);
                filename = [state.registry, '.setup-paths.mat'];
                callback = @() setupPathOwnership('add', state.root, state.root, {owned});
                expected = mode{1};
                if strcmp(expected, 'new')
                    expected = '600';
                else
                    callback(); set_mode(filename, expected);
                end
                callback();
                testCase.verifyEqual(get_mode(filename), expected);
                callback();
                testCase.verifyEqual(get_mode(filename), expected);
                setupPathOwnership('remove-context', state.root);
                testCase.verifyEqual(get_mode(filename), expected);
            end
        end

        function testSymlinkTargetsAreNotReplaced(testCase)
            testCase.assumeTrue(isunix);
            for kind = {'startup', 'pathdef', 'sidecar'}
                for dangling = [false, true]
                    state = isolate(testCase);
                    owned = fullfile(state.root, 'owned'); mkdir(owned);
                    callback = @() setupPathOwnership('add', state.root, state.root, {owned});
                    real_dir = fullfile(state.root, 'real'); mkdir(real_dir);
                    if strcmp(kind{1}, 'startup')
                        filename = state.startup;
                        setenv('OPTIPROFILER_MATLAB_PROBLEM_LIBRARY_PATHDEF', fullfile(state.root, 'missing', 'pathdef.m'));
                    elseif strcmp(kind{1}, 'pathdef')
                        filename = fullfile(state.root, 'pathdef.m');
                    else
                        filename = [state.registry, '.setup-paths.mat'];
                    end
                    [~, name, ext] = fileparts(filename);
                    target = fullfile(real_dir, [name, ext]);
                    if ~dangling
                        if strcmp(kind{1}, 'sidecar')
                            callback(); movefile(filename, target);
                        elseif strcmp(kind{1}, 'pathdef')
                            assert(savepath(target) == 0);
                        else
                            write_bytes(target, uint8('% user dotfile'));
                        end
                        original = read_bytes(target);
                    end
                    assert(system(['/bin/ln -s ', quotePosix(target), ' ', quotePosix(filename)]) == 0);
                    testCase.verifyError(callback, 'OptiProfiler:SymlinkPersistence');
                    testCase.verifyEqual(system(['/bin/test -L ', quotePosix(filename)]), 0);
                    if dangling
                        testCase.verifyFalse(isfile(target));
                    else
                        testCase.verifyEqual(read_bytes(target), original);
                    end
                end
            end
        end

        function testSymlinkParentDirectoryIsAllowed(testCase)
            testCase.assumeTrue(isunix);
            state = isolate(testCase);
            real_dir = fullfile(state.root, 'real'); mkdir(real_dir);
            link_dir = fullfile(state.root, 'parent');
            assert(system(['/bin/ln -s ', quotePosix(real_dir), ' ', quotePosix(link_dir)]) == 0);
            filename = fullfile(link_dir, 'pathdef.m');
            setenv('OPTIPROFILER_MATLAB_PROBLEM_LIBRARY_PATHDEF', filename);
            owned = fullfile(state.root, 'owned'); mkdir(owned);
            setupPathOwnership('add', state.root, state.root, {owned});
            testCase.verifyEqual(get_mode(filename), '600');
            setupPathOwnership('remove-context', state.root);
            testCase.verifyEqual(system(['/bin/test -L ', quotePosix(link_dir)]), 0);
            testCase.verifyEqual(get_mode(filename), '600');
        end
    end
end

function state = isolate(testCase)
    state.root = tempname;
    mkdir(state.root);
    state.registry = fullfile(state.root, 'registry.mat');
    state.startup = fullfile(state.root, 'startup.m');
    old_path = path;
    old_dir = pwd;
    variables = {'OPTIPROFILER_MATLAB_PROBLEM_LIBRARY_REGISTRY', ...
        'OPTIPROFILER_MATLAB_PROBLEM_LIBRARY_PATHDEF', ...
        'OPTIPROFILER_MATLAB_PROBLEM_LIBRARY_STARTUP'};
    values = cellfun(@getenv, variables, 'UniformOutput', false);
    testCase.addTeardown(@() restore(state.root, old_path, old_dir, variables, values));
    setenv(variables{1}, state.registry);
    setenv(variables{2}, fullfile(state.root, 'pathdef.m'));
    setenv(variables{3}, state.startup);
    cd(fullfile(fileparts(mfilename('fullpath')), '..', '..', 'src', 'private'));
end

function restore(root, old_path, old_dir, variables, values)
    path(old_path); cd(old_dir);
    for k = 1:numel(variables)
        setenv(variables{k}, values{k});
    end
    if isfolder(root)
        rmdir(root, 's');
    end
end

function result = on_path(directory)
    result = any(strcmp(strsplit(path, pathsep), directory));
end

function write_bytes(filename, bytes)
    fid = fopen(filename, 'wb');
    cleanup = onCleanup(@() fclose(fid)); %#ok<NASGU>
    assert(fid >= 0);
    assert(fwrite(fid, bytes, 'uint8') == numel(bytes));
end

function bytes = read_bytes(filename)
    fid = fopen(filename, 'rb');
    cleanup = onCleanup(@() fclose(fid)); %#ok<NASGU>
    assert(fid >= 0);
    bytes = fread(fid, Inf, '*uint8').';
end

function set_mode(filename, mode)
    assert(system(['/bin/chmod ', mode, ' ', quotePosix(filename)]) == 0);
end

function mode = get_mode(filename)
    if ismac
        command = '/usr/bin/stat -f %Lp ';
    else
        command = '/usr/bin/stat -c %a ';
    end
    [status, mode] = system([command, quotePosix(filename)]);
    assert(status == 0);
    mode = strtrim(mode);
end

function source = private_assertion()
    source = sprintf(['function assertPrivate(directory)\n' ...
        '[ok,a]=fileattrib(directory);\n' ...
        'assert(ok && a.UserRead && a.UserWrite && a.UserExecute && ' ...
        '~a.GroupRead && ~a.GroupWrite && ~a.GroupExecute && ' ...
        '~a.OtherRead && ~a.OtherWrite && ~a.OtherExecute,' ...
        '''OptiProfilerTest:UnsafeStaging'',''Staging must be private before any write.'');\nend\n']);
end
