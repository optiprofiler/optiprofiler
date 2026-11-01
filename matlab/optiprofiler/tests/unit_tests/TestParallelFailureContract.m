classdef TestParallelFailureContract < matlab.unittest.TestCase
    methods (Test)
        function poolSetupPrintsNothingWhenSilent(testCase)
            % The toolbox functions are shadowed, so that the pool opens or
            % fails on demand, whatever the licence and the pools of the session.
            root = tempname; mkdir(root);
            old_path = path;
            old_directory = pwd;
            testCase.addTeardown(@() rmdir(root, 's'));
            testCase.addTeardown(@() path(old_path));
            testCase.addTeardown(@() cd(old_directory));
            testCase.addTeardown(@() rmappdata(0, 'OptiProfilerPoolFails'));
            writeFile(fullfile(root, 'gcp.m'), sprintf('function pool=gcp(varargin)\npool=[];\nend\n'));
            writeFile(fullfile(root, 'parcluster.m'), sprintf('function c=parcluster(varargin)\nc=struct(''NumWorkers'',4);\nend\n'));
            writeFile(fullfile(root, 'parpool.m'), sprintf(['function parpool(varargin)\n', ...
                'if getappdata(0,''OptiProfilerPoolFails''), error(''Probe:NoPool'',''no pool today''); end\nend\n']));
            addpath(root);
            cd(fullfile(fileparts(mfilename('fullpath')), '..', '..', 'src', 'private'));

            setappdata(0, 'OptiProfilerPoolFails', false);
            text = evalc('success = setupPool(2, true);');
            testCase.verifyTrue(success);
            testCase.verifyEmpty(text);
            setappdata(0, 'OptiProfilerPoolFails', true);
            text = evalc('success = setupPool(2, true);');
            testCase.verifyFalse(success);
            testCase.verifyEmpty(text);
            text = evalc('success = setupPool(2, false);');
            testCase.verifyFalse(success);
            testCase.verifySubstring(text, 'Failed to open a parallel pool with 2 workers.');

            % Without the toolbox.
            rmpath(root);
            toolbox = fileparts(which('parpool'));
            if ~isempty(toolbox), rmpath(toolbox); end
            text = evalc('success = setupPool(2, true);');
            testCase.verifyFalse(success);
            testCase.verifyEmpty(text);
        end
        function ordinaryFailureAndTaskCancellationKeepBorrowedPool(testCase)
            testCase.assumeTrue(logical(license('test','Distrib_Computing_Toolbox')));
            root = tempname; mkdir(root);
            testCase.addTeardown(@() rmdir(root, 's'));
            keys = {'OPTIPROFILER_MATLAB_PROBLEM_LIBRARY_REGISTRY', ...
                'OPTIPROFILER_MATLAB_PROBLEM_LIBRARY_PATHDEF', ...
                'OPTIPROFILER_MATLAB_PROBLEM_LIBRARY_STARTUP'};
            files = {'registry.mat','pathdef.m','startup.m'};
            for k = 1:numel(keys)
                old = getenv(keys{k});
                testCase.addTeardown(@() setenv(keys{k},old));
                setenv(keys{k},fullfile(root,files{k}));
            end
            writeFile(fullfile(root,'followmat_select.m'), ...
                sprintf('function n=followmat_select(~)\nn={''bad'',''good''};\nend\n'));
            writeFile(fullfile(root,'followmat_load.m'), sprintf([ ...
                'function p=followmat_load(n)\n' ...
                'if strcmp(n,''bad''), error(''FollowMat:LoadFailure'',''ordinary failure''); end\n' ...
                'p=Problem(struct(''fun'',@(x)sum(x.^2),''x0'',[1;2],''name'',n));\nend\n']));
            registerProblemLibrary(struct('name','followmat','root',root, ...
                'select_function','followmat_select','load_function','followmat_load'));
            testCase.addTeardown(@() removeProviderPath(root));
            pool = gcp('nocreate');
            if isempty(pool)
                pool = parpool('local',2);
                testCase.addTeardown(@() delete(pool));
            end
            testCase.assumeGreaterThanOrEqual(pool.NumWorkers,2);
            options = struct('plibs',{{'followmat'}},'n_jobs',2,'score_only',true, ...
                'silent',true,'run_plain',false,'n_runs',1,'max_eval_factor',1,'max_tol_order',1);
            scores = benchmark({@initialSolver,@initialSolver},options);
            testCase.verifyTrue(all(isfinite(scores)));
            testCase.verifyEqual(gcp('nocreate'),pool);
            future = parfeval(pool,@pause,0,60);
            cancel(future);
            wait(future,'finished',10);
            testCase.verifyEqual(future.State,'finished');
            testCase.verifyEqual(gcp('nocreate'),pool);
            again = benchmark({@initialSolver,@initialSolver},options);
            testCase.verifyEqual(again,scores);
        end
    end
end
function x=initialSolver(fun,x0)
    fun(x0); x=x0;
end
function writeFile(name,text)
    fid=fopen(name,'w'); c=onCleanup(@()fclose(fid)); %#ok<NASGU>
    fprintf(fid,'%s',text);
end
function removeProviderPath(root)
    if contains([path,pathsep],[root,pathsep]), rmpath(root); end
end
