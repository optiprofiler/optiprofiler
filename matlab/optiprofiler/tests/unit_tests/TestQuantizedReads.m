classdef TestQuantizedReads < matlab.unittest.TestCase
% Which points a quantized feature reads from the oracle, and what it records.
% Mirrors Python's test_quantized_reads.py.
    properties (Constant)
        Mesh = 0.5
        X0 = 0.3
        SnappedX0 = 0.5
        X = 1.2
        SnappedX = 1
    end
    properties (TestParameter)
        groundTruth = {false, true};
    end

    methods (Test)
        function constructionReadsTheTruthPointOnce(testCase, groundTruth)
            [problem, log] = testCase.oracle();
            featured = FeaturedProblem(problem, testCase.quantized(groundTruth), 10, 1);
            truth = testCase.truthPoint(groundTruth, testCase.X0, testCase.SnappedX0);
            testCase.verifyReads(log, truth, truth, truth);
            testCase.verifyEqual(featured.fun_init, truth^2);
            testCase.verifyEqual(featured.maxcv_init, TestQuantizedReads.violation(truth));
        end

        function objectiveQueryReadsObservationAndTruthOnly(testCase, groundTruth)
            [problem, log] = testCase.oracle();
            featured = FeaturedProblem(problem, testCase.quantized(groundTruth), 10, 1);
            TestQuantizedReads.clear(log);
            testCase.verifyEqual(featured.fun(testCase.X), testCase.SnappedX^2);
            truth = testCase.truthPoint(groundTruth, testCase.X, testCase.SnappedX);
            if groundTruth
                expected = testCase.SnappedX;
            else
                expected = [testCase.SnappedX, testCase.X];
            end
            testCase.verifyReads(log, expected, truth, truth);
            testCase.verifyEqual(featured.fun_hist, truth^2);
            testCase.verifyEqual(featured.maxcv_hist, TestQuantizedReads.violation(truth));
        end

        function constraintQueryReadsObservationAndTruthOnly(testCase, groundTruth)
            for channel = {'cub', 'ceq'}
                [problem, log] = testCase.oracle();
                featured = FeaturedProblem(problem, testCase.quantized(groundTruth), 10, 1);
                TestQuantizedReads.clear(log);
                value = featured.(channel{1})(testCase.X);
                offset = 2;
                if strcmp(channel{1}, 'ceq'), offset = 3; end
                testCase.verifyEqual(value, testCase.SnappedX - offset);
                truth = testCase.truthPoint(groundTruth, testCase.X, testCase.SnappedX);
                if groundTruth
                    expected = testCase.SnappedX;
                else
                    expected = [testCase.SnappedX, testCase.X];
                end
                testCase.verifyEqual(log(channel{1}), expected, channel{1});
                testCase.verifyEmpty(log('fun'));
                testCase.verifyEqual(featured.([channel{1}, '_hist']), truth - offset, channel{1});
            end
        end

        function unrecordedConstraintQueryReadsTheObservationOnly(testCase, groundTruth)
            for channel = {'cub', 'ceq'}
                [problem, log] = testCase.oracle();
                featured = FeaturedProblem(problem, testCase.quantized(groundTruth), 10, 1);
                TestQuantizedReads.clear(log);
                featured.(channel{1})(testCase.X, false);
                testCase.verifyEqual(log(channel{1}), testCase.SnappedX, channel{1});
            end
        end

        function violationReadsEachConstraintAtTheTruthPointOnce(testCase, groundTruth)
            [problem, log] = testCase.oracle();
            featured = FeaturedProblem(problem, testCase.quantized(groundTruth), 10, 1);
            TestQuantizedReads.clear(log);
            truth = testCase.truthPoint(groundTruth, testCase.X, testCase.SnappedX);
            testCase.verifyEqual(featured.maxcv(testCase.X), TestQuantizedReads.violation(truth));
            testCase.verifyReads(log, [], truth, truth);
        end

        function problemWithoutConstraintsReadsNoConstraintOracle(testCase)
            log = TestQuantizedReads.newLog();
            problem = Problem(struct('fun', @(x) TestQuantizedReads.logged(log, 'fun', x(1)^2, x), 'x0', testCase.X0));
            featured = FeaturedProblem(problem, testCase.quantized(true), 10, 1);
            testCase.verifyEqual(log('fun'), testCase.SnappedX0);
            testCase.verifyEmpty(featured.cub(testCase.X));
            testCase.verifyEmpty(featured.ceq(testCase.X));
            testCase.verifyEqual(log('fun'), testCase.SnappedX0);
        end

        function cachedQueryAfterTheBudgetReadsNothing(testCase)
            [problem, log] = testCase.oracle();
            featured = FeaturedProblem(problem, testCase.quantized(true), 1, 1);
            featured.fun(testCase.X);
            TestQuantizedReads.clear(log);
            testCase.verifyEqual(featured.fun(0), testCase.SnappedX^2);
            testCase.verifyReads(log, [], [], []);
        end

        function composedViewKeepsOneReadPerChannel(testCase, groundTruth)
            [problem, log] = testCase.oracle();
            stage = struct('name', 'quantized', 'options', struct('mesh_size', testCase.Mesh, ...
                'mesh_type', 'absolute', 'ground_truth', groundTruth));
            noise = struct('name', 'noisy', 'options', struct('noise_level', 0));
            featured = FeaturedProblem(problem, Feature({stage, noise}), 10, 1);
            TestQuantizedReads.clear(log);
            testCase.verifyEqual(featured.fun(testCase.X), testCase.SnappedX^2);
            truth = testCase.truthPoint(groundTruth, testCase.X, testCase.SnappedX);
            % Each view reads its predecessor once: an observation and a reference.
            testCase.verifyReads(log, [testCase.SnappedX, truth], truth, truth);
            testCase.verifyEqual(featured.fun_hist, truth^2);
            TestQuantizedReads.clear(log);
            featured.cub(testCase.X);
            testCase.verifyEqual(log('cub'), [testCase.SnappedX, truth]);
        end
    end

    methods
        function [problem, log] = oracle(testCase)
            log = TestQuantizedReads.newLog();
            problem = Problem(struct('x0', testCase.X0, ...
                'fun', @(x) TestQuantizedReads.logged(log, 'fun', x(1)^2, x), ...
                'cub', @(x) TestQuantizedReads.logged(log, 'cub', x(1) - 2, x), ...
                'ceq', @(x) TestQuantizedReads.logged(log, 'ceq', x(1) - 3, x)));
            % Force the dimension probes before the log is cleared.
            problem.ptype; %#ok<VUNUS>
            TestQuantizedReads.clear(log);
        end

        function feature = quantized(testCase, groundTruth)
            feature = Feature('quantized', struct('mesh_size', testCase.Mesh, ...
                'mesh_type', 'absolute', 'ground_truth', groundTruth));
        end

        function verifyReads(testCase, log, fun, cub, ceq)
            testCase.verifyEqual(log('fun'), reshape(fun, 1, []), 'fun reads');
            testCase.verifyEqual(log('cub'), reshape(cub, 1, []), 'cub reads');
            testCase.verifyEqual(log('ceq'), reshape(ceq, 1, []), 'ceq reads');
        end
    end

    methods (Static)
        function point = truthPoint(groundTruth, unsnapped, snapped)
            point = unsnapped;
            if groundTruth, point = snapped; end
        end

        function value = violation(x)
            value = max([0, x - 2, abs(x - 3)]);
        end

        function log = newLog()
            log = containers.Map('KeyType', 'char', 'ValueType', 'any');
            TestQuantizedReads.clear(log);
        end

        function clear(log)
            for channel = {'fun', 'cub', 'ceq'}
                log(channel{1}) = zeros(1, 0);
            end
        end

        function value = logged(log, channel, value, x)
            log(channel) = [log(channel), x(1)];
        end
    end
end
