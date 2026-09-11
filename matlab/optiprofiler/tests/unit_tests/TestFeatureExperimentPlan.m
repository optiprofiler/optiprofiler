classdef TestFeatureExperimentPlan < matlab.unittest.TestCase
% Public experiment ownership: a reusable feature never owns repetitions.
    methods (Test)
        function canonicalFeatureUsesTopLevelCount(testCase)
            output = tempname;
            mkdir(output);
            testCase.addTeardown(@() rmdir(output, 's'));
            calls = [0, 0];
            feature = Feature('plain');
            problem = Problem(struct('fun', @(x) sum(x.^2), 'x0', [2; 1], 'name', 'PLAN_PROBE'));
            options = struct('feature', feature, 'n_runs', 3, ...
                'problem', problem, 'solver_isrand', [false, true], ...
                'solver_names', {{'deterministic', 'randomized'}}, ...
                'n_jobs', 1, 'max_eval_factor', 2, 'seed', 17, ...
                'score_only', true, 'draw_hist_plots', 'none', ...
                'silent', true, 'savepath', output);
            scores = benchmark({@(fun, x0) probe(1, fun, x0), ...
                @(fun, x0) probe(2, fun, x0)}, options);
            % These are actual user solver invocations, not planned run counts.
            testCase.verifyEqual(calls, [1, 3]);
            testCase.verifySize(scores, [2, 1]);
            testCase.verifyEmpty(feature.stages);
            testCase.verifyTrue(feature.is_identity);
            testCase.verifyEmpty(dir(fullfile(output, '**', 'data_for_loading.mat')));

            function x = probe(index, fun, x0)
                calls(index) = calls(index) + 1;
                fun(x0);
                x = 0.5 * x0;
                fun(x);
            end
        end
    end
end
