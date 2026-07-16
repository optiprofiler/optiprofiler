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
    registry_root = tempname;
    mkdir(registry_root);
    registry_file = fullfile(registry_root, 'problem_libraries.mat');
    setenv('OPTIPROFILER_MATLAB_PROBLEM_LIBRARY_REGISTRY', registry_file);
    testCase.addTeardown(@() setenv( ...
        'OPTIPROFILER_MATLAB_PROBLEM_LIBRARY_REGISTRY', original_registry));
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
