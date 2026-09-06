classdef TestTruthRobustness < matlab.unittest.TestCase
    % Regression seams use actual solver dispatch as well as public oracles.
    methods (Test)
        function testQuantizedInitHistoryAndOutput(testCase)
            src = fileparts(which('Problem'));
            original = pwd;
            cleanup = onCleanup(@() cd(original));
            cd(fullfile(src, 'private'));
            defaults = @getDefaultProfileOptions;
            solveOne = @solveOneProblem;
            cd(original);
            for truth = [false, true]
                for kind = {'u', 'cub', 'ceq', 'mixed', 'linear'}
                    for mesh = {'absolute', 'relative'}
                        p = TestTruthRobustness.makeProblem(kind{1}, .49);
                        feature = Feature('quantized', struct('mesh_size', 1, ...
                            'mesh_type', mesh{1}, 'ground_truth', truth, 'n_runs', 1));
                        q = .49 * ~truth;
                        cv = 0;
                        if ismember(kind{1}, {'cub', 'mixed'}), cv = max(cv, q^2 - .25^2); end
                        if ismember(kind{1}, {'ceq', 'mixed'}), cv = max(cv, q^2); end
                        if strcmp(kind{1}, 'linear'), cv = p.maxcv(p.x0); end
                        solvers = {@TestTruthRobustness.stay};
                        options = defaults(solvers, feature, struct('solver_names', {{'stay'}}, ...
                            'solver_isrand', false, 'project_x0', false, 'max_eval_factor', 1, ...
                            'silent', true, 'seed', 17, 'solver_verbose', 2, ...
                            'score_only', true, 'draw_hist_plots', 'none'));
                        result = solveOne(solvers, p, feature, p.name, length(p.name), options, false, '');
                        for field = {'fun_inits', 'fun_history', 'fun_out'}
                            testCase.verifyEqual(result.(field{1}), q, 'AbsTol', 1e-15);
                        end
                        for field = {'maxcv_inits', 'maxcv_history', 'maxcv_out'}
                            testCase.verifyEqual(result.(field{1}), cv, 'AbsTol', 1e-15);
                        end
                        testCase.verifyEqual(result.n_eval, 1);
                        testCase.verifyFalse(any(result.solver_abnormal_termination, 'all'));
                        testCase.verifyFalse(any(result.solver_output_fallback, 'all'));
                        testCase.verifyEqual(p.x0, .49);
                    end
                end
            end
        end

        function testTruthDoesNotConsumeConstraintBudget(testCase)
            fp = FeaturedProblem(TestTruthRobustness.makeProblem('mixed', .49), ...
                Feature('quantized', struct('mesh_size', 1)), 1, 17);
            fp.fun(fp.x0);
            for i = 1:8
                testCase.verifyEqual(fp.maxcv(1.49), 1);
            end
            % No hidden truth call may consume the solver's first cub/ceq call.
            testCase.verifyEqual(fp.cub(fp.x0), -.25^2);
            testCase.verifyEqual(fp.ceq(fp.x0), 0);
            before = {fp.fun_hist, fp.cub_hist, fp.ceq_hist, fp.maxcv_hist, rng};
            for i = 1:8
                testCase.verifyEqual(fp.maxcv(1.49), 1);
            end
            testCase.verifyEqual({fp.fun_hist, fp.cub_hist, fp.ceq_hist, fp.maxcv_hist, rng}, before);
            testCase.verifyEqual([fp.n_eval_fun, fp.n_eval_cub, fp.n_eval_ceq], [1, 1, 1]);
        end

        function testConstraintCountFailureIsNotUnconstrained(testCase)
            for kind = {'cub', 'ceq'}
                opts = struct('fun', @(x) x(1), 'x0', 1);
                opts.(kind{1}) = @TestTruthRobustness.sometimesFails;
                p = Problem(opts);
                testCase.verifyEqual(p.ptype, 'n');
                p.x0 = -1;
                if strcmp(kind{1}, 'cub')
                    testCase.verifyError(@() p.m_nonlinear_ub, 'OptiProfilerTest:ConstraintUnavailable');
                else
                    testCase.verifyError(@() p.m_nonlinear_eq, 'OptiProfilerTest:ConstraintUnavailable');
                end
                % ptype already has a conservative fallback for unknown counts.
                testCase.verifyEqual(p.ptype, 'n');
                p.x0 = 2;
                testCase.verifyEqual(p.ptype, 'n');
            end
        end

        function testAbsentAndEmptyConstraintsRemainValid(testCase)
            p = Problem(struct('fun', @(x) 0, 'x0', 0));
            testCase.verifyEqual([p.m_nonlinear_ub, p.m_nonlinear_eq], [0, 0]);
            p = Problem(struct('fun', @(x) 0, 'x0', 0, 'cub', @(x) [], 'ceq', @(x) []));
            testCase.verifyEqual(p.ptype, 'u');
        end

        function testAffineCountProbeUsesOriginalCoordinates(testCase)
            p = Problem(struct('fun', @(x) x(1)^2, 'x0', 1, ...
                'cub', @TestTruthRobustness.sometimesFails));
            feature = Feature('custom', struct('mod_affine', @TestTruthRobustness.reverseCoordinates));
            fp = FeaturedProblem(p, feature, 3, 17);
            testCase.verifyEqual(fp.x0, -1);
            testCase.verifyEqual(fp.m_nonlinear_ub, 2);
            testCase.verifyEqual(fp.cub(fp.x0), [1; 1]);
            testCase.verifyEqual(fp.maxcv(fp.x0), 1);
        end

        function testDocumentedConditionFactorMatchesMatrix(testCase)
            for n = [1, 2, 4]
                for rotated = [false, true]
                    p = Problem(struct('fun', @(x) sum(x.^2), 'x0', ones(n, 1)));
                    feature = Feature('linearly_transformed', struct('condition_factor', 2, 'rotated', rotated));
                    [A, ~, invA] = feature.modifier_affine(17, p);
                    expected = 2^sqrt(n);
                    if n == 1, expected = 1; end
                    testCase.verifyEqual(cond(A), expected, 'RelTol', 1e-14);
                    testCase.verifyEqual(A * invA, eye(n), 'AbsTol', 1e-14);
                end
            end
        end

        function testNonquantifiablePreservesNaN(testCase)
            p = Problem(struct('fun', @(x) 0, 'x0', 0, ...
                'cub', @(x) [NaN; -Inf; 0; Inf], 'ceq', @(x) [NaN; -Inf; 0; Inf]));
            feature = Feature('nonquantifiable_constraints');
            testCase.verifyEqual(feature.modifier_cub(p.x0, 17, p, 0), [NaN; 0; 0; 1]);
            testCase.verifyEqual(feature.modifier_ceq(p.x0, 17, p, 0), [NaN; 1; 0; 1]);
            fp = FeaturedProblem(p, feature, 2, 17);
            fp.fun(fp.x0);
            testCase.verifyTrue(isnan(fp.maxcv_init));
            testCase.verifyTrue(all(isnan(fp.maxcv_hist)));
        end
    end

    methods (Static)
        function p = makeProblem(kind, x0)
            opts = struct('fun', @(x) x(1), 'x0', x0, 'name', 'QUANTIZED_TRUTH');
            if ismember(kind, {'cub', 'mixed'}), opts.cub = @(x) x(1)^2 - .25^2; end
            if ismember(kind, {'ceq', 'mixed'}), opts.ceq = @(x) x(1)^2; end
            if strcmp(kind, 'linear')
                opts.xl = .2; opts.xu = .4; opts.aub = 1; opts.bub = .3;
            end
            p = Problem(opts);
        end

        function xout = stay(fun, x0, varargin)
            xout = x0;
            fun(xout);
            if numel(varargin) == 8
                varargin{7}(xout);
                varargin{8}(xout);
            end
            assert(isequal(xout, .49));
        end

        function c = sometimesFails(x)
            if x(1) < 0
                error('OptiProfilerTest:ConstraintUnavailable', 'Constraint evaluation is unavailable.');
            end
            c = [x(1); x(1)^2];
        end

        function [A, b, invA] = reverseCoordinates(~, ~)
            A = -1; b = 0; invA = -1;
        end
    end
end
