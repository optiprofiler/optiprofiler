classdef TestProblemLibraryRegistry < matlab.unittest.TestCase
    methods (Test)

        function testBundledS2mpj(testCase)
            cleanup = isolateRegistry(testCase); %#ok<NASGU>
            library = resolveProblemLibrary('s2mpj');

            testCase.verifyEqual(library.name, 's2mpj');
            testCase.verifyEqual(library.role, 'bundled-default');
            testCase.verifyEqual(func2str(library.select), 's2mpj_select');
            testCase.verifyEqual(func2str(library.load), 's2mpj_load');
        end

        function testExplicitRegistrationPersists(testCase)
            cleanup = isolateRegistry(testCase); %#ok<NASGU>
            provider_root = createProvider(testCase, 'sample', true);
            registration = providerRegistration('sample', provider_root);

            registered = registerProblemLibrary(registration);
            resolved = resolveProblemLibrary('sample');

            testCase.verifyEqual(registered.root, resolved.root);
            testCase.verifyEqual(resolved.select(struct()), {'DEMO'});
            testCase.verifyEqual(resolved.load('DEMO'), 'DEMO');
            testCase.verifyEqual(func2str(resolved.select), 'sample_select');
        end

        function testUnregisterRemovesRegistrationAndPathButPreservesSource(testCase)
            cleanup = isolateRegistry(testCase); %#ok<NASGU>
            provider_root = createProvider(testCase, 'sample', true);
            registerProblemLibrary(providerRegistration('sample', provider_root));
            startup_file = getenv('OPTIPROFILER_MATLAB_PROBLEM_LIBRARY_STARTUP');
            writeText(startup_file, sprintf([ ...
                'addpath(''%s'');\t%% optiprofiler\n', ...
                'disp(''keep this line'');\n'], provider_root));

            removed = unregisterProblemLibrary('sample');

            testCase.verifyEqual(removed.root, provider_root);
            testCase.verifyError(@() resolveProblemLibrary('sample'), ...
                "MATLAB:resolveProblemLibrary:notRegistered");
            testCase.verifyFalse(any(strcmp(strsplit(path, pathsep), provider_root)));
            testCase.verifyEqual(exist(provider_root, 'dir'), 7);
            pathdef_file = getenv('OPTIPROFILER_MATLAB_PROBLEM_LIBRARY_PATHDEF');
            testCase.verifyEqual(exist(pathdef_file, 'file'), 2);
            testCase.verifyFalse(contains(fileread(pathdef_file), provider_root));
            startup_text = fileread(startup_file);
            testCase.verifyFalse(contains(startup_text, provider_root));
            testCase.verifySubstring(startup_text, 'keep this line');
        end

        function testUnregisterStaleRootPreservesOtherStaleEntries(testCase)
            cleanup = isolateRegistry(testCase); %#ok<NASGU>
            first_root = createProvider(testCase, 'first_stale', true);
            second_root = createProvider(testCase, 'second_stale', true);
            registerProblemLibrary(providerRegistration('first_stale', first_root));
            registerProblemLibrary(providerRegistration('second_stale', second_root));
            cleanupProvider(first_root);
            cleanupProvider(second_root);

            unregisterProblemLibrary('first_stale');

            registry = load(getenv('OPTIPROFILER_MATLAB_PROBLEM_LIBRARY_REGISTRY'), ...
                'registrations');
            testCase.verifyEqual({registry.registrations.name}, {'second_stale'});
            testCase.verifyEqual(registry.registrations.root, second_root);
        end

        function testBundledProvidersCannotBeUnregistered(testCase)
            cleanup = isolateRegistry(testCase); %#ok<NASGU>
            testCase.verifyError(@() unregisterProblemLibrary('s2mpj'), ...
                "MATLAB:unregisterProblemLibrary:bundledProvider");
            testCase.verifyError(@() unregisterProblemLibrary('custom'), ...
                "MATLAB:unregisterProblemLibrary:bundledProvider");
        end

        function testLockedExternalProviderCanBeUnregistered(testCase)
            cleanup = isolateRegistry(testCase); %#ok<NASGU>
            provider_root = createProvider(testCase, 'solar', true);
            registerProblemLibrary(providerRegistration('solar', provider_root));

            unregisterProblemLibrary('solar');

            testCase.verifyError(@() resolveProblemLibrary('solar'), ...
                "MATLAB:resolveProblemLibrary:notInstalled");
        end

        function testRegistrationIsIdempotentButRejectsReplacement(testCase)
            cleanup = isolateRegistry(testCase); %#ok<NASGU>
            first_root = createProvider(testCase, 'sample', true);
            registration = providerRegistration('sample', first_root);
            first = registerProblemLibrary(registration);
            repeated = registerProblemLibrary(registration);
            testCase.verifyEqual(first.root, repeated.root);

            second_root = createProvider(testCase, 'sample', true);
            replacement = providerRegistration('sample', second_root);
            testCase.verifyError(@() registerProblemLibrary(replacement), ...
                "MATLAB:registerProblemLibrary:alreadyRegistered");
            testCase.verifyEqual(fileparts(which('sample_load')), first_root);
            testCase.verifyFalse(any(strcmp(strsplit(path, pathsep), second_root)));
        end

        function testBundledRoleCannotBeClaimedByExternalRegistration(testCase)
            cleanup = isolateRegistry(testCase); %#ok<NASGU>
            provider_root = createProvider(testCase, 'pretend', true);
            registration = providerRegistration('pretend', provider_root);
            registration.role = 'bundled-example';
            testCase.verifyError(@() registerProblemLibrary(registration), ...
                "MATLAB:problemLibraryRegistry:bundledRegistration");
        end

        function testUnregisterErrorsAreStable(testCase)
            cleanup = isolateRegistry(testCase); %#ok<NASGU>
            testCase.verifyError(@() unregisterProblemLibrary('bad-name'), ...
                "MATLAB:unregisterProblemLibrary:nameNotValid");
            testCase.verifyError(@() unregisterProblemLibrary('missing'), ...
                "MATLAB:unregisterProblemLibrary:notRegistered");

            registry_file = getenv('OPTIPROFILER_MATLAB_PROBLEM_LIBRARY_REGISTRY');
            invalid = 1; %#ok<NASGU>
            save(registry_file, 'invalid');
            testCase.verifyError(@() unregisterProblemLibrary('missing'), ...
                "MATLAB:problemLibraryRegistry:invalidFile");
        end

        function testMissingRequiredFunction(testCase)
            cleanup = isolateRegistry(testCase); %#ok<NASGU>
            provider_root = createProvider(testCase, 'missing', false);
            registration = providerRegistration('missing', provider_root);

            testCase.verifyError(@() registerProblemLibrary(registration), ...
                "MATLAB:problemLibraryRegistry:functionNotFound");
        end

        function testUnsupportedPlatform(testCase)
            cleanup = isolateRegistry(testCase); %#ok<NASGU>
            provider_root = createProvider(testCase, 'platform', true);
            registration = providerRegistration('platform', provider_root);
            if ispc
                registration.platforms = {'linux'};
            else
                registration.platforms = {'windows'};
            end
            registerProblemLibrary(registration);

            testCase.verifyError(@() resolveProblemLibrary('platform'), ...
                "MATLAB:resolveProblemLibrary:unsupportedPlatform");
        end

        function testStaleOptionalRootDoesNotBreakBundledDefault(testCase)
            cleanup = isolateRegistry(testCase); %#ok<NASGU>
            provider_root = createProvider(testCase, 'stale', true);
            registerProblemLibrary(providerRegistration('stale', provider_root));
            cleanupProvider(provider_root);

            library = resolveProblemLibrary('s2mpj');
            testCase.verifyEqual(library.name, 's2mpj');
        end

        function testLockedProviderRejectsDifferentFunctions(testCase)
            cleanup = isolateRegistry(testCase); %#ok<NASGU>
            provider_root = createProvider(testCase, 'sample', true);
            registration = providerRegistration('sample', provider_root);
            registration.name = 'solar';

            testCase.verifyError(@() registerProblemLibrary(registration), ...
                "MATLAB:problemLibraryRegistry:lockedContractMismatch");
        end

        function testUnregisteredFolderIsNotDiscovered(testCase)
            cleanup = isolateRegistry(testCase); %#ok<NASGU>
            provider_root = createProvider(testCase, 'unregistered', true);
            addpath(provider_root);
            testCase.addTeardown(@() removePathIfPresent(provider_root));

            testCase.verifyError(@() resolveProblemLibrary('unregistered'), ...
                "MATLAB:resolveProblemLibrary:notRegistered");
        end

        function testResolvedProviderRunsOnParallelWorker(testCase)
            testCase.assumeTrue(license('test', 'Distrib_Computing_Toolbox'));
            testCase.assumeEqual(exist('parpool', 'file'), 2);

            cleanup = isolateRegistry(testCase); %#ok<NASGU>
            provider_root = createProvider(testCase, 'parallel', true);
            registerProblemLibrary(providerRegistration('parallel', provider_root));
            library = resolveProblemLibrary('parallel');

            pool = gcp('nocreate');
            if isempty(pool)
                pool = parpool('local', 2);
                testCase.addTeardown(@() deletePoolIfPresent(pool));
            end
            future = parfeval(pool, library.load, 1, 'DEMO');
            testCase.verifyEqual(fetchOutputs(future), 'DEMO');
        end

    end
end


function cleanup = isolateRegistry(testCase)
    original_registry = getenv('OPTIPROFILER_MATLAB_PROBLEM_LIBRARY_REGISTRY');
    original_pathdef = getenv('OPTIPROFILER_MATLAB_PROBLEM_LIBRARY_PATHDEF');
    original_startup = getenv('OPTIPROFILER_MATLAB_PROBLEM_LIBRARY_STARTUP');
    registry_root = tempname;
    mkdir(registry_root);
    registry_file = fullfile(registry_root, 'problem_libraries.mat');
    pathdef_file = fullfile(registry_root, 'pathdef.m');
    startup_file = fullfile(registry_root, 'startup.m');
    setenv('OPTIPROFILER_MATLAB_PROBLEM_LIBRARY_REGISTRY', registry_file);
    setenv('OPTIPROFILER_MATLAB_PROBLEM_LIBRARY_PATHDEF', pathdef_file);
    setenv('OPTIPROFILER_MATLAB_PROBLEM_LIBRARY_STARTUP', startup_file);
    testCase.addTeardown(@() setenv( ...
        'OPTIPROFILER_MATLAB_PROBLEM_LIBRARY_REGISTRY', original_registry));
    testCase.addTeardown(@() setenv( ...
        'OPTIPROFILER_MATLAB_PROBLEM_LIBRARY_PATHDEF', original_pathdef));
    testCase.addTeardown(@() setenv( ...
        'OPTIPROFILER_MATLAB_PROBLEM_LIBRARY_STARTUP', original_startup));
    testCase.addTeardown(@() removeDirectoryIfPresent(registry_root));
    cleanup = onCleanup(@() rehash);
end


function root = createProvider(testCase, name, include_load)
    root = tempname;
    mkdir(root);
    select_file = fullfile(root, [name, '_select.m']);
    writeText(select_file, sprintf([ ...
        'function names = %s_select(~)\n', ...
        '    names = {''DEMO''};\n', ...
        'end\n'], name));
    if include_load
        load_file = fullfile(root, [name, '_load.m']);
        writeText(load_file, sprintf([ ...
            'function problem = %s_load(name)\n', ...
            '    problem = name;\n', ...
            'end\n'], name));
    end
    collect_info_file = fullfile(root, [name, '_collect_info.m']);
    writeText(collect_info_file, sprintf([ ...
        'function info = %s_collect_info()\n', ...
        '    info = struct();\n', ...
        'end\n'], name));
    testCase.addTeardown(@() cleanupProvider(root));
end


function registration = providerRegistration(name, root)
    registration = struct( ...
        'name', name, ...
        'root', root, ...
        'select_function', [name, '_select'], ...
        'load_function', [name, '_load']);
end


function writeText(filename, contents)
    file_id = fopen(filename, 'w');
    if file_id == -1
        error('Cannot create test provider file %s.', filename);
    end
    cleanup = onCleanup(@() fclose(file_id));
    fprintf(file_id, '%s', contents);
    clear cleanup
end


function removePathIfPresent(root)
    if contains([path, pathsep], [root, pathsep])
        rmpath(root);
    end
end


function removeDirectoryIfPresent(root)
    if exist(root, 'dir') == 7
        rmdir(root, 's');
    end
end


function cleanupProvider(root)
    removePathIfPresent(root);
    removeDirectoryIfPresent(root);
end


function deletePoolIfPresent(pool)
    if ~isempty(pool) && isvalid(pool)
        delete(pool);
    end
end
