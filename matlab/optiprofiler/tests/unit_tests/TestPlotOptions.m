classdef TestPlotOptions < matlab.unittest.TestCase
% User-native properties must reach public saved figures, not affect scores.
    methods (TestMethodSetup)
        function useRepositorySources(testCase)
            root = fileparts(mfilename('fullpath'));
            testCase.applyFixture(matlab.unittest.fixtures.PathFixture(fullfile(root,'..','..','src')));
        end
    end
    methods (Test)
        function nativePropertiesReachPublicFigures(testCase)
            root = tempname; mkdir(root);
            registry = getenv('OPTIPROFILER_MATLAB_PROBLEM_LIBRARY_REGISTRY');
            oldpath = path;
            guard = onCleanup(@() cleanup(root,registry,oldpath)); %#ok<NASGU>
            setenv('OPTIPROFILER_MATLAB_PROBLEM_LIBRARY_REGISTRY',fullfile(root,'registry.mat'));
            library = fullfile(root,'library'); mkdir(library);
            TestReferencePoints.write(fullfile(library,'style_probe_select.m'), sprintf( ...
                'function names=style_probe_select(~)\n names={''P''};\nend\n'));
            TestReferencePoints.write(fullfile(library,'style_probe_load.m'), sprintf( ...
                'function p=style_probe_load(name)\n p=Problem(struct(''fun'',@(x)sum(x.^2),''x0'',2,''name'',name));\nend\n'));
            registerProblemLibrary(struct('name','style_probe','root',library, ...
                'select_function','style_probe_select','load_function','style_probe_load'));
            options = struct('plibs',{{'style_probe'}},'n_jobs',1,'n_runs',1,'max_eval_factor',2, ...
                'max_tol_order',1,'run_plain',false,'draw_hist_plots','sequential','silent',true, ...
                'savepath',root,'benchmark_id','style');
            solvers = {@TestReferencePoints.stay,@TestReferencePoints.zero};
            control = options; control.score_only = true;
            baseline = benchmark(solvers,control);
            options.plot_options = struct('text',struct('FontSize',13),'axes',struct('FontSize',17), ...
                'xlabel',struct('FontSize',21),'ylabel',struct('FontSize',22), ...
                'title',struct('FontSize',23),'legend',struct('FontSize',18), ...
                'line',struct('LineWidth',2));
            styled = benchmark(solvers,options);
            testCase.verifyEqual(styled,baseline);
            files = dir(fullfile(root,'**','perf_hist_1.fig'));
            testCase.assertNumElements(files,1);
            figure = openfig(fullfile(files(1).folder,files(1).name),'invisible');
            figure_guard = onCleanup(@() close(figure)); %#ok<NASGU>
            axes = findall(figure,'Type','axes');
            testCase.verifyEqual(axes.FontSize,17);
            testCase.verifyEqual(axes.XLabel.FontSize,21);
            testCase.verifyEqual(axes.YLabel.FontSize,22);
            legend = findall(figure,'Type','legend');
            testCase.verifyEqual(legend.FontSize,18);
            titles = findall(figure,'Type','text');
            testCase.verifyTrue(any([titles.FontSize] == 23));
            % Inspect every native artifact, including histories and the
            % nested summary: styling must not stop at single profiles.
            files = dir(fullfile(root,'**','*.fig'));
            testCase.verifyGreaterThan(numel(files),6);
            for k = 1:numel(files)
                saved = openfig(fullfile(files(k).folder,files(k).name),'invisible');
                saved_guard = onCleanup(@() close(saved)); %#ok<NASGU>
                all_axes = findall(saved,'Type','axes');
                for j = 1:numel(all_axes)
                    testCase.verifyEqual(all_axes(j).FontSize,17);
                    testCase.verifyEqual(all_axes(j).XLabel.FontSize,21);
                    testCase.verifyEqual(all_axes(j).YLabel.FontSize,22);
                end
                clear saved_guard;
            end
            archive = dir(fullfile(root,'style','**','data_for_loading.mat'));
            testCase.assertNumElements(archive,1);
            archive_path = fullfile(archive.folder,archive.name);
            original_bytes = TestReferencePoints.bytes(archive_path);
            old_directory = pwd;
            directory_guard = onCleanup(@() cd(old_directory)); %#ok<NASGU>
            cd(root);
            loaded = benchmark(struct('load','latest','benchmark_id','style','savepath',fullfile(root,'reload'), ...
                'silent',true,'n_jobs',1,'max_tol_order',1,'run_plain',false, ...
                'plot_options',struct('xlabel',struct('FontSize',25))));
            testCase.verifyEqual(loaded,styled);
            testCase.verifyEqual(TestReferencePoints.bytes(archive_path),original_bytes);
            restyled = dir(fullfile(root,'reload','**','perf_hist_1.fig'));
            testCase.assertNumElements(restyled,1);
            saved = openfig(fullfile(restyled.folder,restyled.name),'invisible');
            saved_guard = onCleanup(@() close(saved)); %#ok<NASGU>
            axes = findall(saved,'Type','axes');
            testCase.verifyEqual(axes.XLabel.FontSize,25);
        end
        function unknownGroupFailsBeforeExecution(testCase)
            problem = Problem(struct('fun',@(x)sum(x.^2),'x0',1));
            options = struct('problem',problem,'n_jobs',1,'score_only',true,'silent',true, ...
                'plot_options',struct('unknown_group',struct()));
            testCase.verifyError(@() benchmark({@TestReferencePoints.stay,@TestReferencePoints.zero},options), ...
                'OptiProfiler:PlotOptions:UnknownGroup');
        end
        function invalidValueFailsBeforeExecution(testCase)
            % Values are checked on unparented native objects before any
            % solver runs; no validation figure is created or left open.
            setappdata(0, 'OptiProfilerPlotOptionsSolverCalls', 0);
            testCase.addTeardown(@() rmappdata(0, 'OptiProfilerPlotOptionsSolverCalls'));
            figures_before = numel(findall(groot, 'Type', 'figure'));
            problem = Problem(struct('fun',@(x)sum(x.^2),'x0',1));
            options = struct('problem',problem,'n_jobs',1,'score_only',true,'silent',true, ...
                'plot_options',struct('xlabel',struct('FontSize','big')));
            counted = @(fun, x0) countedStay(fun, x0);
            testCase.verifyError(@() benchmark({counted, counted}, options), ...
                'OptiProfiler:PlotOptions:InvalidProperty');
            testCase.verifyEqual(getappdata(0, 'OptiProfilerPlotOptionsSolverCalls'), 0);
            testCase.verifyEqual(numel(findall(groot, 'Type', 'figure')), figures_before);
            % Valid values pass the same early check and still run.
            options.plot_options = struct('xlabel',struct('FontSize',18),'legend',struct('Location','southeast'), ...
                'line',struct('LineWidth',2),'axes',struct('XGrid','on'));
            scores = benchmark({counted, counted}, options);
            testCase.verifyEqual(numel(scores), 2);
            testCase.verifyGreaterThan(getappdata(0, 'OptiProfilerPlotOptionsSolverCalls'), 0);
            testCase.verifyEqual(numel(findall(groot, 'Type', 'figure')), figures_before);
        end
        function labelTemplatesSubstituteOnlyTheTolerance(testCase)
            % Labels were formatted by sprintf: \a, \b, \t, \n became control
            % characters and a single % ended the label.
            root = tempname; mkdir(root);
            old_directory = pwd;
            testCase.addTeardown(@() rmdir(root, 's'));
            testCase.addTeardown(@() cd(old_directory));
            cd(root);
            TestLoadReplot.makeExperiment(root, 'bench', 2, 1);
            options = struct('load', '20200101_000000', 'benchmark_id', 'bench', 'max_tol_order', 1, ...
                'silent', true, 'draw_hist_plots', 'none', 'savepath', fullfile(root, 'latex'), ...
                'ylabel_performance_profile', 'Solved $\alpha$ at $%s$, 100\%', ...
                'ylabel_data_profile', 'Data $\\beta$ at $%s$ (100%%)');
            benchmark(options);
            testCase.verifyEqual(savedLabel(fullfile(root, 'latex'), 'perf_hist_1.fig'), 'Solved $\alpha$ at $10^{-1}$, 100\%');
            testCase.verifyEqual(savedLabel(fullfile(root, 'latex'), 'data_hist_1.fig'), 'Data $\beta$ at $10^{-1}$ (100%)');
            testCase.verifyEqual(savedLabel(fullfile(root, 'latex'), 'log-ratio_hist_1.fig'), 'Log-ratio profiles ($\mathrm{tol} = 10^{-1}$)');
            % The same rule without the LaTeX interpreter.
            options.savepath = fullfile(root, 'plain');
            options.ylabel_performance_profile = 'C:\temp\new: 50% solved at %s';
            options.ylabel_data_profile = "tol %s / %s";
            options.plot_options = struct('ylabel', struct('Interpreter', 'none'));
            benchmark(options);
            testCase.verifyEqual(savedLabel(fullfile(root, 'plain'), 'perf_hist_1.fig'), 'C:\temp\new: 50% solved at 10^{-1}');
            testCase.verifyEqual(savedLabel(fullfile(root, 'plain'), 'data_hist_1.fig'), 'tol 10^{-1} / 10^{-1}');
        end

        function unknownPropertyFailsBeforeExecution(testCase)
            problem = Problem(struct('fun',@(x)sum(x.^2),'x0',1));
            options = struct('problem',problem,'n_jobs',1,'score_only',true,'silent',true, ...
                'plot_options',struct('axes',struct('UnknownProperty',1)));
            testCase.verifyError(@() benchmark({@TestReferencePoints.stay,@TestReferencePoints.zero},options), ...
                'OptiProfiler:PlotOptions:UnknownProperty');
        end
    end
end

function cleanup(root,registry,oldpath)
    path(oldpath);
    setenv('OPTIPROFILER_MATLAB_PROBLEM_LIBRARY_REGISTRY',registry);
    if isfolder(root), rmdir(root,'s'); end
end

function label = savedLabel(root, name)
    files = dir(fullfile(root, '**', name));
    assert(isscalar(files), 'Expected one %s under %s.', name, root);
    figure = openfig(fullfile(files(1).folder, files(1).name), 'invisible');
    guard = onCleanup(@() close(figure)); %#ok<NASGU>
    axes = findall(figure, 'Type', 'axes');
    label = char(axes(1).YLabel.String);
end

function x = countedStay(fun, x0)
    setappdata(0, 'OptiProfilerPlotOptionsSolverCalls', getappdata(0, 'OptiProfilerPlotOptionsSolverCalls') + 1);
    fun(x0);
    x = x0;
end
