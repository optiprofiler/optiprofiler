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
        function silentRunPrintsNothingOfAFailingSolver(testCase)
            root = tempname; mkdir(root);
            testCase.addTeardown(@() rmdir(root, 's'));
            levels = {struct(), struct('solver_verbose',0), struct('solver_verbose',1), struct('solver_verbose',2)};
            printed = {'', '', silentLevelWarning(), silentLevelWarning()};
            for k = 1:numel(levels)
                name = sprintf('level%d',k);
                [text,report] = runFailingSolver(root,name,failingRunOptions(1,true,levels{k}),'');
                testCase.verifyEqual(printedText(text),printed{k});
                verifyFailureRecorded(testCase,report,root,name);
            end

            % The record holds the effective level, and a replay keeps silent.
            file = dir(fullfile(root,'level1','**','options_refined.mat'));
            path_record = fullfile(file(1).folder,file(1).name);
            saved = load(path_record);
            testCase.verifyEqual(saved.options_refined.solver_verbose,0);
            [text,report] = runFailingSolver(root,'replay',loadBenchmarkOptions(path_record),'');
            testCase.verifyEmpty(text);
            verifyFailureRecorded(testCase,report,root,'replay');
        end
        function silentWorkersPrintNothingOfAFailingSolver(testCase)
            testCase.assumeTrue(logical(license('test','Distrib_Computing_Toolbox')));
            root = tempname; mkdir(root);
            testCase.addTeardown(@() rmdir(root, 's'));
            pool = gcp('nocreate');
            if isempty(pool)
                pool = parpool('local',2);
                testCase.addTeardown(@() delete(pool));
            end
            testCase.assumeGreaterThanOrEqual(pool.NumWorkers,2);
            where = fullfile(root,'solver_process.txt');

            % What a worker prints reaches the text that evalc returns.
            text = runFailingSolver(root,'shown',failingRunOptions(2,false,struct()),where);
            testCase.verifySubstring(text,'An error occurred while solving');

            levels = {struct(), struct('solver_verbose',2)};
            printed = {'', silentLevelWarning()};
            for k = 1:numel(levels)
                name = sprintf('workers%d',k);
                delete(where);
                [text,report] = runFailingSolver(root,name,failingRunOptions(2,true,levels{k}),where);
                testCase.verifyEqual(printedText(text),printed{k});
                testCase.verifyEqual(fileread(where),'1');
                verifyFailureRecorded(testCase,report,root,name);
            end
        end
    end
end
function x=initialSolver(fun,x0)
    fun(x0); x=x0;
end
function x=failingSolver(fun,x0,where) %#ok<STOUT>
    fun(x0);
    if ~isempty(where), writeFile(where,sprintf('%d',~isempty(getCurrentTask()))); end
    fprintf('SOLVER OUTPUT\n');
    error('Probe:SolverFailure','deliberate solver failure');
end
function text=silentLevelWarning()
    text='The option `solver_verbose` will be set to 0 because silent is true.';
end
function text=printedText(text)
    % A batch session marks a warning up as [\bWarning: text]\b.
    text=regexprep(strtrim(erase(text,char(8))),'^\[?Warning: (.*?)\]?$','$1');
end
function options=failingRunOptions(n_jobs,silent,level)
    options=struct('plibs',{{'custom'}},'problem_names',{{'custom1'}},'n_jobs',n_jobs,'silent',silent, ...
        'n_runs',1,'max_eval_factor',1,'max_tol_order',1,'run_plain',false,'draw_hist_plots','none', ...
        'solver_names',{{'failing','stay'}});
    for field=fieldnames(level)'
        options.(field{1})=level.(field{1});
    end
end
function [text,report]=runFailingSolver(root,name,options,where)
    options.savepath=root;
    options.benchmark_id=name;
    options.report_path=fullfile(root,[name,'.json']);
    solvers={@(fun,x0)failingSolver(fun,x0,where),@initialSolver}; %#ok<NASGU>
    state=warning('off','backtrace');
    restore=onCleanup(@()warning(state)); %#ok<NASGU>
    text=evalc('benchmark(solvers,options);');
    report=jsondecode(fileread(options.report_path));
end
function verifyFailureRecorded(testCase,report,root,name)
    testCase.verifyEqual(report.configuration.effective.profile_options.solver_verbose,0);
    problems=report.problems;
    if ~iscell(problems), problems=num2cell(problems); end
    abnormal=[];
    for i=1:numel(problems)
        runs=problems{i}.runs;
        if ~iscell(runs), runs=num2cell(runs); end
        for j=1:numel(runs)
            if runs{j}.abnormal_termination, abnormal(end+1)=runs{j}.solver_index; end %#ok<AGROW>
        end
    end
    testCase.verifyEqual(abnormal,1);
    file=dir(fullfile(root,name,'**','report.txt'));
    section=extractAfter(fileread(fullfile(file(1).folder,file(1).name)),'## Solver runs that terminated abnormally');
    section=extractBefore([section,newline,'## '],[newline,'## ']);
    testCase.verifySubstring(section,'custom1');
end
function writeFile(name,text)
    fid=fopen(name,'w'); c=onCleanup(@()fclose(fid)); %#ok<NASGU>
    fprintf(fid,'%s',text);
end
function removeProviderPath(root)
    if contains([path,pathsep],[root,pathsep]), rmpath(root); end
end
