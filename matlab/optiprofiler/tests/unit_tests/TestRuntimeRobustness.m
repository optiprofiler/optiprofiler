classdef TestRuntimeRobustness < matlab.unittest.TestCase
%TESTRUNTIMEROBUSTNESS Runtime failures must not masquerade as valid results.
    methods (Test)
        function selectorFailureIsNotAnEmptySelection(testCase)
            old = pwd;
            cleanup = onCleanup(@() cd(old));
            cd(fullfile(fileparts(mfilename('fullpath')), '..', '..', 'src', 'private'));
            library = struct('name', 'broken_selector', 'select', ...
                @(options) error('test:SelectorFailed', 'Injected selector failure.'));
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

        function summaryProfilesKeepSelectedOrder(testCase)
            old = pwd;
            directory_guard = onCleanup(@() cd(old)); %#ok<NASGU>
            cd(fullfile(fileparts(mfilename('fullpath')), '..', '..', 'src', 'private'));
            options = getDefaultProfileOptions({@(varargin) [], @(varargin) []}, ...
                Feature('plain'), struct('n_jobs', 1, 'silent', true));
            options.ylabel_performance_profile = 'perf %s';
            options.ylabel_data_profile = 'data %s';
            options.ylabel_log_ratio_profile = 'ratio %s';
            labels = {'perf 0.1', 'data 0.1', 'ratio 0.1'};
            for mask = 0:7
                enabled = logical(bitget(mask, 1:3));
                summary = figure('Visible', 'off');
                summary_guard = onCleanup(@() close(summary));
                axs = cell(1, nnz(enabled));
                for row = 1:numel(axs)
                    axs{row} = subplot(numel(axs), 1, row, 'Parent', summary);
                end
                curve = struct('perf', {cell(2, 2)}, 'data', {cell(2, 2)}, ...
                    'log_ratio', {cell(1, 2)});
                [a, b, c, ~, ok] = drawProfiles([1, 2; 3, 4], [2; 3], ...
                    {'one', 'two'}, '0.1', axs, any(enabled), enabled(1), ...
                    enabled(2), enabled(3), options, curve);
                figures_guard = onCleanup(@() close([a, b, c]));
                testCase.verifyTrue(ok, sprintf('Summary selection %d failed.', mask));
                selected_labels = labels(enabled);
                for row = 1:numel(axs)
                    testCase.verifyEqual(axs{row}.YLabel.String, selected_labels{row});
                    testCase.verifyNotEmpty(axs{row}.Children);
                end
                clear figures_guard summary_guard;
            end
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

        function exportFailureThatDestroysFiguresKeepsSvgFallback(testCase)
            % When the graphics process exits during a native export, its
            % figures are gone. Closing them must not abort the SVG fallback.
            old=pwd; old_path=path; old_warnings=warning;
            work=tempname; mkdir(work);
            cleanup=onCleanup(@() restoreState(old,old_warnings,old_path,work));
            addpath(fullfile(fileparts(mfilename('fullpath')),'..','..','src'));
            doubles=fullfile(work,'doubles'); mkdir(doubles);
            writeText(fullfile(doubles,'exportgraphics.m'),sprintf([ ...
                'function exportgraphics(varargin)\n' ...
                '%% Test double: the graphics process exits during an export.\n' ...
                'state = getappdata(groot, ''TestRuntimeRobustnessExport'');\n' ...
                'state.calls = state.calls + 1;\n' ...
                'setappdata(groot, ''TestRuntimeRobustnessExport'', state);\n' ...
                'figures = findall(groot, ''Type'', ''figure'');\n' ...
                'for k = 1:numel(figures)\n' ...
                '    if ~any(figures(k) == state.before), delete(figures(k)); end\n' ...
                'end\n' ...
                'error(''TestRuntimeRobustness:GraphicsExited'', ''Simulated exit of the graphics process.'');\n' ...
                'end\n']));
            before=findall(groot,'Type','figure');
            setappdata(groot,'TestRuntimeRobustnessExport',struct('calls',0,'before',before));
            state_cleanup=onCleanup(@() rmappdata(groot,'TestRuntimeRobustnessExport'));
            addpath(doubles);
            TestLoadReplot.makeExperiment(work,'bench',2,1); cd(work);
            options=TestLoadReplot.loadOptions('savepath',work,'draw_hist_plots','none','max_tol_order',2);
            options.report_path=fullfile(work,'fallback.json');
            output=evalc('scores=benchmark(options);');
            state=getappdata(groot,'TestRuntimeRobustnessExport');
            testCase.assumeGreaterThan(state.calls,0,'Native export is not attempted in this graphics environment.');
            testCase.verifyNumElements(scores,2);
            testCase.verifyTrue(contains(output,'Native profile export failed; using SVG fallback'));
            testCase.verifyNotEmpty(dir(fullfile(work,'bench','**','summary.html')));
            testCase.verifyEqual(findall(groot,'Type','figure'),before);
            % The portable output is a degradation of the requested PDF/FIG
            % output, so the report must not call the rendering complete.
            report=jsondecode(fileread(options.report_path));
            testCase.verifyEqual(report.stages.rendering.status,'partial');
            testCase.verifyEqual(report.status,'partial');
            diagnostics=report.diagnostics;
            if isstruct(diagnostics), diagnostics=num2cell(diagnostics); end
            codes=cellfun(@(d) d.code,diagnostics,'UniformOutput',false);
            testCase.verifyTrue(any(strcmp(codes,'native_profile_export_failed')),strjoin(codes,', '));
        end

        function stalledNativeExportIsBounded(testCase)
            % A native export that waits for an event that never comes ends
            % at the limit: the watchdog deletes the figure and the call fails
            % with OptiProfiler:NativeExportTimeout. No timer is left behind.
            old=pwd; cleanup=onCleanup(@() cd(old));
            cd(fullfile(fileparts(mfilename('fullpath')),'..','..','src','private'));
            fig=figure('Visible','off');
            started=tic;
            testCase.verifyError(@() boundedFigureCall(fig,@() waitfor(fig),1),'OptiProfiler:NativeExportTimeout');
            testCase.verifyLessThan(toc(started),60);
            testCase.verifyFalse(isgraphics(fig));
            testCase.verifyEmpty(exportTimers());
            % A stalled exportgraphics fails once its figure is gone; the
            % caller still sees the timeout, not that secondary error.
            fig=figure('Visible','off');
            testCase.verifyError(@() boundedFigureCall(fig,@() failWhenClosed(fig),1),'OptiProfiler:NativeExportTimeout');
            testCase.verifyFalse(isgraphics(fig));
            testCase.verifyEmpty(exportTimers());
        end

        function nativeExportIsNotAttemptedWithoutWatchdog(testCase)
            % If the watchdog cannot be started (here because of an invalid
            % limit), the export is not run, because it might never return.
            old=pwd; cleanup=onCleanup(@() cd(old));
            cd(fullfile(fileparts(mfilename('fullpath')),'..','..','src','private'));
            fig=figure('Visible','off');
            figure_cleanup=onCleanup(@() delete(fig(isgraphics(fig))));
            testCase.verifyError(@() boundedFigureCall(fig,@() setappdata(fig,'exported',true),-1), ...
                'OptiProfiler:NativeExportUnbounded');
            testCase.verifyFalse(isappdata(fig,'exported'));
            testCase.verifyTrue(isgraphics(fig));
            testCase.verifyEmpty(exportTimers());
        end

        function failingPageReloadKeepsTimeLimit(testCase)
            % The reload attempt and the time limit share one timer, and an
            % error in a timer callback stops that timer. With a test double of
            % webFigureViewRecovery that always fails, the export must still end
            % at the limit.
            old=pwd; work=tempname; mkdir(work);
            cleanup=onCleanup(@() restoreDirectory(old,work));
            copyfile(fullfile(fileparts(mfilename('fullpath')),'..','..','src','private','boundedFigureCall.m'),work);
            % Exercise the recovery callback even on platforms where the real
            % capability gate disables it. This double is isolated to WORK.
            fid=fopen(fullfile(work,'webFigureRecoverySupported.m'),'w');
            fprintf(fid,'%s\n','function supported=webFigureRecoverySupported()', ...
                'supported=true;','end');
            fclose(fid);
            fid=fopen(fullfile(work,'webFigureViewRecovery.m'),'w');
            fprintf(fid,'%s\n','function status=webFigureViewRecovery(~)', ...
                'error(''test:Recovery'',''Recovery failed.'');','end');
            fclose(fid);
            cd(work);
            fig=figure('Visible','off');
            figure_cleanup=onCleanup(@() delete(fig(isgraphics(fig))));
            % Without the limit, this guard ends the wait and the test fails.
            guard=timer('StartDelay',90,'TimerFcn',@(~,~) delete(fig(isgraphics(fig))));
            guard_cleanup=onCleanup(@() removeTimer(guard));
            start(guard);
            started=tic;
            testCase.verifyError(@() boundedFigureCall(fig,@() waitfor(fig),31),'OptiProfiler:NativeExportTimeout');
            testCase.verifyLessThan(toc(started),60);
            testCase.verifyEmpty(exportTimers());
        end

        function nativeExportWatchdogMatchesRecoveryScope(testCase)
            % Outside the validated release, keep the old single-shot path:
            % no early callback and no access to MATLAB's webwindow internals.
            old=pwd; cleanup=onCleanup(@() cd(old));
            cd(fullfile(fileparts(mfilename('fullpath')),'..','..','src','private'));
            fig=figure('Visible','off');
            figure_cleanup=onCleanup(@() delete(fig(isgraphics(fig))));
            expected='singleShot';
            if strcmp(computer('arch'),'glnxa64') && strcmp(version('-release'),'2026a')
                expected='fixedRate';
            end
            boundedFigureCall(fig,@() testCase.verifyEqual(char(string(get(exportTimers(),'ExecutionMode'))),expected),60);
            testCase.verifyEmpty(exportTimers());
        end

        function figureViewRecoveryOnlyConcernsOwnInvisibleFigure(testCase)
            % The page recovery is offered only for invisible figures on the
            % tested Linux release, and never acts on a ready or another figure.
            old=pwd; cleanup=onCleanup(@() cd(old));
            cd(fullfile(fileparts(mfilename('fullpath')),'..','..','src','private'));
            shown=figure('Visible','on'); hidden=figure('Visible','off');
            figure_cleanup=onCleanup(@() delete([shown(isgraphics(shown)),hidden(isgraphics(hidden))]));
            supported=strcmp(computer('arch'),'glnxa64') && strcmp(version('-release'),'2026a');
            testCase.verifyEqual(webFigureRecoverySupported(),supported);
            testCase.verifyEqual(webFigureViewRecovery(shown),'unsupported');
            if ~supported
                testCase.verifyEqual(webFigureViewRecovery(hidden),'unsupported');
                return;
            end
            testCase.verifyEqual(webFigureViewRecovery(groot),'figure gone');
            plot(axes(hidden),1:3);
            waitForFigureView(hidden);
            page=figurePageUrl(hidden);
            testCase.verifyEqual(webFigureViewRecovery(hidden),'view ready');
            testCase.verifyEqual(figurePageUrl(hidden),page);
            delete(hidden);
            testCase.verifyEqual(webFigureViewRecovery(hidden),'figure gone');
        end

        function failedFigurePageRecoversDuringExport(testCase)
            % A figure page whose load failed leaves exportgraphics waiting.
            % The recovery reloads only that figure's page, and the same export
            % then completes; unrelated figures are left as they are. Elsewhere
            % the recovery is inactive, as figureViewRecoveryOnlyConcernsOwnInvisibleFigure
            % checks; this test then has nothing to check and passes, because
            % the runtime CI requires every case to complete.
            if ~(strcmp(computer('arch'),'glnxa64') && strcmp(version('-release'),'2026a'))
                testCase.log(matlab.unittest.Verbosity.Terse, ...
                    'Not applicable: the figure page recovery is only enabled on Linux with R2026a.');
                return;
            end
            old=pwd; work=tempname; mkdir(work);
            cleanup=onCleanup(@() restoreDirectory(old,work));
            cd(fullfile(fileparts(mfilename('fullpath')),'..','..','src','private'));
            other=breakFigurePage(figure('Visible','off'));
            shown=figure('Visible','on'); plot(axes(shown),1:3); waitForFigureView(shown);
            fig=breakFigurePage(figure('Visible','off'));
            figure_cleanup=onCleanup(@() delete([other(isgraphics(other)),shown(isgraphics(shown)),fig(isgraphics(fig))]));
            testCase.assertEqual(char(string(get(fig,'FigureViewReady'))),'off','The test could not break the figure page.');
            other_page=figurePageUrl(other); shown_page=figurePageUrl(shown);
            file=fullfile(work,'recovered.pdf');
            boundedFigureCall(fig,@() exportgraphics(fig,file,'ContentType','vector'),120);
            testCase.verifyTrue(isfile(file));
            testCase.verifyEqual(figurePageUrl(other),other_page);
            testCase.verifyEqual(char(string(get(other,'FigureViewReady'))),'off');
            testCase.verifyEqual(figurePageUrl(shown),shown_page);
            testCase.verifyEmpty(exportTimers());
        end

        function reloadedFigurePageIsReadyOnReturn(testCase)
            % After a reload, the recovery returns only once the reloaded view
            % is ready (or after 20 seconds), so that a waiting export resumes
            % with the layout of a page that loaded normally.
            if ~(strcmp(computer('arch'),'glnxa64') && strcmp(version('-release'),'2026a'))
                testCase.log(matlab.unittest.Verbosity.Terse, ...
                    'Not applicable: the figure page recovery is only enabled on Linux with R2026a.');
                return;
            end
            old=pwd; cleanup=onCleanup(@() cd(old));
            cd(fullfile(fileparts(mfilename('fullpath')),'..','..','src','private'));
            fig=breakFigurePage(figure('Visible','off'));
            figure_cleanup=onCleanup(@() delete(fig(isgraphics(fig))));
            testCase.assertEqual(char(string(get(fig,'FigureViewReady'))),'off','The test could not break the figure page.');
            started=tic;
            testCase.verifyEqual(webFigureViewRecovery(fig),'reloaded');
            testCase.verifyEqual(char(string(get(fig,'FigureViewReady'))),'on');
            testCase.verifyLessThan(toc(started),25);
        end

        function boundedExportKeepsOrdinaryResultsAndErrors(testCase)
            old=pwd; work=tempname; mkdir(work);
            cleanup=onCleanup(@() restoreDirectory(old,work));
            cd(fullfile(fileparts(mfilename('fullpath')),'..','..','src','private'));
            fig=figure('Visible','off');
            figure_cleanup=onCleanup(@() delete(fig(isgraphics(fig))));
            plot(axes(fig),1:3);
            file=fullfile(work,'bounded.pdf');
            exportGraphicsBounded(fig,file,'ContentType','vector');
            testCase.verifyTrue(isfile(file));
            testCase.verifyTrue(isgraphics(fig));
            testCase.verifyError(@() boundedFigureCall(fig,@() error('test:Export','Export failed.'),60),'test:Export');
            testCase.verifyTrue(isgraphics(fig));
            testCase.verifyEmpty(exportTimers());
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

function failWhenClosed(fig)
    % Like a stalled exportgraphics: wait for the figure's view, then fail
    % because the figure is gone.
    waitfor(fig);
    error('TestRuntimeRobustness:FigureClosed','The figure is invalid or has been closed.');
end

function timers=exportTimers()
    timers=timerfindall('Tag','OptiProfilerNativeExportWatchdog');
end

function removeTimer(t)
    stop(t); delete(t);
end

function waitForFigureView(fig)
    waited=tic;
    while toc(waited)<60 && ~strcmp(char(string(get(fig,'FigureViewReady'))),'on')
        pause(0.05);
    end
end

function url=figurePageUrl(fig)
    % Test-only use of MATLAB internals: the address shown by the window of
    % MATLAB's embedded browser that hosts FIG's page.
    url='';
    windows=matlab.internal.cef.webwindowmanager.instance.windowList;
    for k=1:numel(windows)
        if contains(char(windows(k).URL),['channel=/uifigure/' char(string(fig.Uuid))])
            url=char(windows(k).URL);
        end
    end
end

function fig=breakFigurePage(fig)
    % Test-only use of MATLAB internals: before the new figure's page reports
    % ready, navigate that page to an unreachable local address (refused by
    % Chromium), as when one of its own requests fails. MATLAB still associates
    % the window with the figure's page, as after a real failure.
    waited=tic;
    window=[];
    while isempty(window) && toc(waited)<10
        windows=matlab.internal.cef.webwindowmanager.instance.windowList;
        for k=1:numel(windows)
            if contains(char(windows(k).URL),['channel=/uifigure/' char(string(fig.Uuid))])
                window=windows(k);
            end
        end
        pause(0.005);
    end
    script=['(function(){if(location.href.indexOf("cefComponentContainer.html")<0){return "wait";}' ...
        'location.replace("https://127.0.0.1:9/");return "sent";})()'];
    while ~isempty(window) && toc(waited)<10
        try
            if contains(char(string(window.executeJS(script,2))),'sent')
                break;
            end
        catch
        end
        pause(0.01);
    end
    plot(axes(fig),1:3);
    pause(1);
end

function writeText(file,text)
    fid=fopen(file,'w'); cleanup=onCleanup(@() fclose(fid)); fprintf(fid,'%s',text);
end
