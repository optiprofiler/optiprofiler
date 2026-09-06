classdef TestHistoryExtremes < matlab.unittest.TestCase
%TESTHISTORYEXTREMES Real-helper display tests; no solver or provider needed.
    properties
        ProcessHistory
        DrawHistory
        PrepareHistory
        ComputeShift
        WorkDir
        SavedPath
    end
    methods (TestMethodSetup)
        function obtainPrivateHelpers(testCase)
            testCase.SavedPath = path;
            source = fullfile(fileparts(mfilename('fullpath')), '..', '..', 'src');
            addpath(source);
            old = pwd;
            restore = onCleanup(@() cd(old));
            cd(fullfile(source, 'private'));
            testCase.ProcessHistory = @processHistYaxes;
            testCase.DrawHistory = @drawHist;
            testCase.PrepareHistory = @prepareHistoryPlotData;
            testCase.ComputeShift = @computeHistoryYShift;
            testCase.WorkDir = tempname;
            mkdir(testCase.WorkDir);
        end
    end
    methods (TestMethodTeardown)
        function restoreState(testCase)
            path(testCase.SavedPath);
            rmdir(testCase.WorkDir, 's');
        end
    end
    methods (Test)
        function maskIsLocalToItsRun(testCase)
            original = reshape([1,1,100,100,2,2,200,200,3,3,Inf,Inf], [2,2,3]);
            expected = original;
            expected(:, 2, 3) = 250;
            actual = testCase.ProcessHistory(original, [1;100]);
            testCase.verifyEqual(actual, expected);
            testCase.verifyTrue(all(isinf(original(:, 2, 3)), 'all'));
        end
        function invalidInitialAndOverflowStayFinite(testCase)
            cases = {[-1e308, Inf, 1e308], [1,Inf,3], [NaN,Inf,-Inf]};
            inits = [0, Inf, NaN];
            for k = 1:numel(cases)
                original = repmat(reshape(cases{k}, 1,1,[]), 2,2,1);
                before = original;
                actual = testCase.ProcessHistory(original, [inits(k);inits(k)]);
                testCase.verifyTrue(all(isfinite(actual), 'all'));
                testCase.verifyTrue(isequaln(original, before));
            end
        end
        function annotationsStateDisplayLimitAndCount(testCase)
            original = reshape([1e300, NaN, 2], 1,1,[]);
            [~, note] = testCase.ProcessHistory(original, 1);
            testCase.verifySubstring(note, 'Display clipped at +/-1e100: 1 entries');
            testCase.verifySubstring(note, 'Nonfinite placeholders: 1 entries');
        end
        function extremeHistoriesExportAndStayInView(testCase)
            cases = {[1,1e290,2], [-1e308,0,1e308], [NaN,Inf,-Inf], ...
                [1e-308,1e-200,1e-150], [1e308,1e308,1e308], ...
                [-1e308,-1e308,-1e308], [0,0,0], [1e-308,1e-308,1e-308]};
            for k = 1:numel(cases)
                for is_cum = [false,true]
                    for errorbar = {'minmax','meanstd'}
                        original = repmat(reshape(cases{k},1,1,[]),2,2,1);
                        testCase.renderAndCheck(original, is_cum, errorbar{1});
                    end
                end
            end
        end
        function singletonDimensions(testCase)
            for shape = {[1,2,3],[2,2,1],[1,2,1]}
                testCase.renderAndCheck(ones(shape{1}), true, 'minmax');
            end
        end
        function meanstdKeepsSampleNormalization(testCase)
            % Unequal runs expose the established N-1 versus Python N policy.
            runs = {[1,4,2], [0,0,0;4,8,2], [1,4,2;3,8,6;5,6,10]};
            means = {[1,4,2], [2,4,1], [3,6,6]};
            deviations = {[0,0,0], sqrt(2)*[2,4,1], [2,2,4]};
            options = struct('errorbar_type','meanstd','hist_aggregation','min');
            for k = 1:numel(runs)
                n_runs = size(runs{k},1);
                y = reshape(runs{k},1,n_runs,3);
                before = y;
                shift = testCase.ComputeShift(y, options);
                expected_shift = 0;
                if k == 2
                    % Mean minus sample std is negative at the second evaluation.
                    expected_shift = (sqrt(2)-1)*4;
                end
                testCase.verifyEqual(shift, expected_shift, 'AbsTol',1e-14);
                for is_cum = [false,true]
                    expected = [means{k}; means{k}-deviations{k}; means{k}+deviations{k}] + shift;
                    if is_cum
                        expected = cummin(expected,2);
                    end
                    [x,m,l,u,count] = testCase.PrepareHistory(y,is_cum,shift,3*ones(1,n_runs),options);
                    testCase.verifyEqual(x{1},1:3);
                    testCase.verifyEqual([m{1};l{1};u{1}],expected,'AbsTol',1e-14);
                    testCase.verifyEqual(count,n_runs);
                end
                testCase.verifyEqual(y,before);
            end
        end
    end
    methods (Access = private)
        function renderAndCheck(testCase, original, is_cum, errorbar)
            before = original;
            fig = figure('Visible','off', 'Position',[50,50,800,500]);
            cleanup = onCleanup(@() close(fig));
            ax = axes(fig);
            opts = struct('errorbar_type',errorbar,'hist_aggregation','min', ...
                'line_colors',[0,0.4,0.7;0.9,0.4,0], 'line_styles',{{'-','--'}}, ...
                'line_widths',[1,1], 'xlabel_data_profile','Evaluations / (n+1)');
            names = {'a','b'};
            n_solvers = size(original,1);
            n_runs = size(original,2);
            n_evals = size(original,3);
            inits = reshape(original(1,:,1), [],1);
            testCase.DrawHistory(original, original, original, inits, inits, inits, ...
                names(1:n_solvers), {ax}, is_cum, 'u', 2, ...
                repmat(n_evals,n_solvers,n_runs), opts, 500);
            file = fullfile(testCase.WorkDir,'history.pdf');
            exportgraphics(fig, file, 'ContentType','vector');
            fid = fopen(file, 'r');
            signature = fread(fid, 4, '*char')';
            fclose(fid);
            testCase.verifyEqual(signature, '%PDF');
            limits = ax.YLim;
            testCase.verifyTrue(all(isfinite(limits)) && limits(1)<limits(2));
            testCase.verifyTrue(all(isfinite(ax.YTick)));
            lines = findall(ax,'Type','line');
            for i = 1:numel(lines)
                displayed = lines(i).YData;
                testCase.verifyTrue(all(isfinite(displayed)));
                testCase.verifyTrue(all(displayed >= limits(1) & displayed <= limits(2)));
            end
            testCase.verifyTrue(isequaln(original, before));
        end
    end
end
