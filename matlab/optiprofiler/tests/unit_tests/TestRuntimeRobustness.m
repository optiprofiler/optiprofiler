classdef TestRuntimeRobustness < matlab.unittest.TestCase
%TESTRUNTIMEROBUSTNESS Runtime failures must not masquerade as valid results.
    methods (Test)
        function selectorFailureIsNotAnEmptySelection(testCase)
            old = pwd;
            cleanup = onCleanup(@() cd(old));
            cd(fullfile(fileparts(mfilename('fullpath')), '..', '..', 'src', 'private'));
            fixture = fullfile(fileparts(mfilename('fullpath')), '..', 'fixtures', 'selector_failure');
            addpath(fixture);
            path_cleanup = onCleanup(@() rmpath(fixture));
            library = 'broken_selector';
            for names = {{}, {'DO_NOT_RUN'}}
                options = struct('problem_names', {names{1}}, 'excludelist', {{}});
                testCase.verifyError(@() solveAllProblems({}, library, [], options, ...
                    struct('silent', true), false, ''), 'OptiProfiler:ProblemSelectionFailed');
            end
        end

        function loadFailureRestoresWarningState(testCase)
            old = pwd;
            old_warnings = warning;
            old_path = path;
            work = tempname;
            mkdir(work);
            cleanup = onCleanup(@() restoreState(old, old_warnings, old_path, work));
            folder = fullfile(work, 'bench', 'experiment', 'test_log');
            mkdir(folder);
            fid = fopen(fullfile(folder, 'time_stamp_20200101_000000.txt'), 'w'); fclose(fid);
            fid = fopen(fullfile(folder, 'data_for_loading.mat'), 'w');
            fprintf(fid, 'Intentionally not a MAT file.'); fclose(fid);
            addpath(fullfile(fileparts(mfilename('fullpath')), '..', '..', 'src'));
            cd(work);
            warning('on', 'all'); warning('off', 'test:PreserveDisabled');
            before = warning;
            failed = false;
            try
                benchmark(struct('load', '20200101_000000', 'benchmark_id', 'bench', ...
                    'score_only', true, 'silent', true));
            catch
                failed = true;
            end
            testCase.verifyTrue(failed);
            testCase.verifyEqual(warning, before);
        end

        function nestedCallerCopiesOnlyImmediateSource(testCase)
            old=pwd; old_path=path; old_warnings=warning;
            work=tempname; mkdir(work);
            cleanup=onCleanup(@() restoreState(old,old_warnings,old_path,work));
            test_file=[mfilename('fullpath'),'.m'];
            addpath(fullfile(fileparts(test_file),'..','..','src'));
            TestLoadReplot.makeExperiment(work,'bench',2,1); cd(work);
            options=TestLoadReplot.loadOptions('savepath',work,'draw_hist_plots','none', ...
                'n_jobs',1,'summarize_performance_profiles',false, ...
                'summarize_data_profiles',false,'summarize_log_ratio_profiles',false);
            % Nested helpers reproduce the real caller stack, rather than
            % testing copyfile on a preselected scalar path.
            scores=nestedLoadCaller(options);
            copies=dir(fullfile(work,'bench','*','test_log','TestRuntimeRobustness.m'));
            testCase.assertNumElements(copies,1);
            testCase.verifyEqual(fileread(fullfile(copies(1).folder,copies(1).name)),fileread(test_file));
            testCase.verifyTrue(contains(fileread(fullfile(copies(1).folder,'README.txt')),copies(1).name));
            testCase.verifyNumElements(scores,2);
        end

        function scoreOnlyDoesNotAllocateFigures(testCase)
            old = pwd;
            cleanup = onCleanup(@() cd(old));
            cd(fullfile(fileparts(mfilename('fullpath')), '..', '..', 'src', 'private'));
            options = struct('score_only', true, 'semilogx', true);
            curve = struct('perf', {cell(2, 2)}, 'data', {cell(2, 2)}, ...
                'log_ratio', {cell(1, 2)});
            before = findall(groot, 'Type', 'figure');
            [a, b, c, curves] = drawProfiles(reshape([1, 2, 3, 4], 2, 2), ...
                [2; 3], {'one', 'two'}, '0.1', {}, false, false, false, false, options, curve);
            testCase.verifyEmpty(a); testCase.verifyEmpty(b); testCase.verifyEmpty(c);
            testCase.verifyEqual(findall(groot, 'Type', 'figure'), before);
            testCase.verifyNotEmpty(curves.perf{1, 1});
        end

        function exclusiveDirectoryNeverReusesExistingOutput(testCase)
            old=pwd; work=tempname; mkdir(work);
            cleanup=onCleanup(@() restoreDirectory(old,work));
            cd(fullfile(fileparts(mfilename('fullpath')),'..','..','src','private'));
            folder=fullfile(work,'space '' quote ; $ literal');
            testCase.verifyTrue(exclusiveDirectory(folder));
            testCase.verifyFalse(exclusiveDirectory(folder));
            sentinel=fullfile(folder,'original'); writeText(sentinel,'unchanged');
            testCase.verifyFalse(exclusiveDirectory(folder));
            testCase.verifyEqual(fileread(sentinel),'unchanged');
        end

        function atomicReplacePreservesOldFileOnFailure(testCase)
            old=pwd; work=tempname; mkdir(work);
            cleanup=onCleanup(@() restoreDirectory(old,work));
            cd(fullfile(fileparts(mfilename('fullpath')),'..','..','src','private'));
            source=fullfile(work,'new '' source'); target=fullfile(work,'old target');
            writeText(source,'new content'); writeText(target,'old content');
            testCase.verifyError(@() atomicReplaceFile(fullfile(work,'missing'),target),'OptiProfiler:AtomicReplace');
            testCase.verifyEqual(fileread(target),'old content');
            atomicReplaceFile(source,target);
            testCase.verifyEqual(fileread(target),'new content');
            testCase.verifyFalse(isfile(source));
            testCase.verifyEmpty(dir(fullfile(work,'.optiprofiler-replace.*')));
        end

        function renderingFailurePreservesCurvesAndNumericalErrors(testCase)
            old=pwd; cleanup=onCleanup(@() cd(old));
            cd(fullfile(fileparts(mfilename('fullpath')),'..','..','src','private'));
            options=getDefaultProfileOptions({@(f,x) x,@(f,x) x},Feature('plain'), ...
                struct('score_only',false,'silent',true,'n_jobs',1));
            curve=struct('perf',{cell(2,2)},'data',{cell(2,2)},'log_ratio',{cell(1,2)});
            before=findall(groot,'Type','figure');
            % Invalid summary axes fault only the rendering stage. This also
            % covers native allocation failure on no-graphics installations.
            [a,b,c,result,ok]=drawProfiles([1 2;3 4],[2;3],{'one','two'},'0.1', ...
                {struct('not_an_axes',true)},true,true,false,false,options,curve);
            testCase.verifyFalse(ok);
            testCase.verifyEmpty(a); testCase.verifyEmpty(b); testCase.verifyEmpty(c);
            testCase.verifyNotEmpty(result.perf{1,1});
            testCase.verifyEqual(findall(groot,'Type','figure'),before);
            failed=false;
            try
                drawProfiles([1 2;3 4],[],{'one','two'},'0.1',{},false,false,false,false,options,curve);
            catch
                failed=true;
            end
            testCase.verifyTrue(failed,'Numerical curve errors must not become graphics fallbacks.');
        end

        function loadNotesDoNotInferTruthFromCustomStamps(testCase)
            old=pwd; work=tempname; mkdir(work);
            cleanup=onCleanup(@() restoreDirectory(old,work));
            cd(fullfile(fileparts(mfilename('fullpath')),'..','..','src','private'));
            results={struct('plib','synthetic','fun_histories',1,'maxcv_histories',0, ...
                'merit_histories',1,'feature_stamp','arbitrary_user_stamp')};
            report=fullfile(work,'report.txt'); readme=fullfile(work,'README.txt');
            appendRuntimeReportNotes(results,Feature('plain'),true,report,readme);
            text=fileread(report);
            testCase.verifyTrue(contains(text,'does not reevaluate points'));
            testCase.verifyTrue(contains(text,'Older quantized archives'));
            testCase.verifyFalse(contains(text,'ground_truth='));
            testCase.verifyTrue(contains(fileread(readme),'Saved truth channels'));
            % Cross a chunk boundary and include repeated/padded tail entries.
            values=zeros(1,1000004);
            values([1,1000000,1000001,1000002,1000003,1000004])=[1e101,-1e101,NaN,Inf,-Inf,-1e101];
            results{1}.fun_histories=values;
            chunk_report=fullfile(work,'chunk-counts.txt');
            appendRuntimeReportNotes(results,Feature('plain'),true,chunk_report,readme);
            testCase.verifyTrue(contains(fileread(chunk_report), ...
                'fun_histories: clipped finite entries=3; nonfinite placeholders=3'));
            testCase.verifyEqual(results{1}.fun_histories,values);
            for truth=[false,true]
                report=fullfile(work,sprintf('new-%d.txt',truth));
                appendRuntimeReportNotes(results,Feature('quantized',struct('ground_truth',truth)),false,report,readme);
                testCase.verifyTrue(contains(fileread(report),sprintf('ground_truth=%d',truth)));
            end
        end

        function portableSummaryRespectsOptions(testCase)
            old=pwd; work=tempname; mkdir(work);
            cleanup=onCleanup(@() restoreDirectory(old,work));
            cd(fullfile(fileparts(mfilename('fullpath')),'..','..','src','private'));
            record=struct('perf',{repmat({[0 1;0 1]},2,2)}, ...
                'data',{repmat({[0 1;0 1]},2,2)}, ...
                'log_ratio',{{[1 2;-1 -.5],[1 2;-1 -.5]}});
            curves={struct('hist',record,'out',record)};
            options=struct('semilogx',true,'summarize_performance_profiles',false, ...
                'summarize_data_profiles',true,'summarize_log_ratio_profiles',false, ...
                'summarize_output_based_profiles',false);
            exportPortableProfiles(curves,{'one','two'},options,work);
            summary=fileread(fullfile(work,'summary.svg'));
            testCase.verifyTrue(contains(summary,'hist Data profile'));
            testCase.verifyFalse(contains(summary,'Performance profile'));
            testCase.verifyFalse(contains(summary,'Log-ratio profile'));
            testCase.verifyFalse(contains(summary,'out Data profile'));
            testCase.verifyTrue(isfile(fullfile(work,'log_ratio_out_1.svg')));
            index=fileread(fullfile(work,'summary.html'));
            testCase.verifyTrue(contains(index,'href="log_ratio_out_1.svg"'));
            testCase.verifyFalse(contains(index,work));
            folder=fullfile(work,'no-summary'); mkdir(folder);
            options.summarize_data_profiles=false;
            exportPortableProfiles(curves,{'one','two'},options,folder);
            testCase.verifyFalse(isfile(fullfile(folder,'summary.svg')));
            testCase.verifyTrue(isfile(fullfile(folder,'summary.html')));
        end

        function portableSvgEscapesTextAndBreaksNonfiniteLines(testCase)
            old=pwd; work=tempname; mkdir(work);
            cleanup=onCleanup(@() restoreDirectory(old,work));
            cd(fullfile(fileparts(mfilename('fullpath')),'..','..','src','private'));
            panel=struct('title','title <&>', 'xlabel','x', 'ylabel','y', ...
                'curves',{{[1 2 3 4;1 NaN 2 3]}},'labels',{{'solver <&>'}},'note','display note');
            file=fullfile(work,'chart.svg'); writeCurveSvg(file,panel,'heading');
            first=fileread(file); writeCurveSvg(file,panel,'heading');
            testCase.verifyEqual(fileread(file),first);
            testCase.verifyTrue(contains(first,'solver &lt;&amp;&gt;'));
            testCase.verifyFalse(contains(first,'NaN'));
            testCase.verifyEqual(numel(strfind(first,'<polyline')),2);
            panel.curves={[-realmax,realmax;-realmax,realmax]};
            writeCurveSvg(file,panel,'Extreme finite coordinates');
            testCase.verifyFalse(contains(fileread(file),'NaN'));
            testCase.verifyFalse(contains(fileread(file),'Inf'));
            panel.curves={[realmin*eps;realmin*eps]};
            writeCurveSvg(file,panel,'Subnormal constant');
            testCase.verifyFalse(contains(fileread(file),'NaN'));
            testCase.verifyTrue(contains(fileread(file),'<circle'));
            panel.curves={[1 2;2 4]}; panel.include_zero=true;
            panel.labels={repmat('W',1,225)};
            writeCurveSvg(file,panel,'Long legend with a zero reference');
            drawing=fileread(file);
            testCase.verifyTrue(contains(drawing,'stroke-dasharray="4 3"'));
            testCase.verifyFalse(contains(drawing,repmat('W',1,61)));
            testCase.verifyEqual(numel(strfind(drawing,repmat('W',1,60))),3);
            testCase.verifyError(@() writeCurveSvg(work,panel,'Cannot overwrite a directory'),'OptiProfiler:SvgOutput');
        end
    end
end

function scores=nestedLoadCaller(options)
    scores=immediateLoadCaller(options);
end

function scores=immediateLoadCaller(options)
    scores=benchmark(options);
end

function restoreState(old, old_warnings, old_path, work)
    cd(old); warning(old_warnings); path(old_path); rmdir(work, 's');
end

function restoreDirectory(old,work)
    cd(old); rmdir(work,'s');
end

function writeText(file,text)
    fid=fopen(file,'w'); cleanup=onCleanup(@() fclose(fid)); fprintf(fid,'%s',text);
end
