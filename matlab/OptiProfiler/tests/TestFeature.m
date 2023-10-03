classdef TestFeature < matlab.unittest.TestCase
    methods (Static)

        function f = rosen(x)
            f = sum(100*(x(2:end) - x(1:end-1).^2).^2 + (1 - x(1:end-1)).^2);
        end

    end

    methods (Test)

        function testPlain(testCase)
            % Test that the plain feature works correctly.
            nValues = [1, 10, 100];
            seedValues = [0, 1, 2];

            for n = nValues
                for seed = seedValues
                    % Generate random data
                    rng(seed);
                    x = randn(n, 1);
                    f = TestFeature.rosen(x);
                    
                    % Test the plain feature
                    feature = Feature("PLAIN");
                    testCase.verifyEqual(feature.name, "plain");
                    testCase.verifyEqual(feature.options, struct('n_runs', int32(1)));
                    testCase.verifyEqual(feature.modifier(x, f), f);
                end
            end
        end

        function testRegularized(testCase)
            % Test that the regularized feature works correctly.
            nValues = [1, 10, 100];
            seedValues = [0, 1, 2];

            for n = nValues
                for seed = seedValues
                    % Generate random data
                    rng(seed);
                    x = randn(n, 1);
                    f = TestFeature.rosen(x);
                    
                    % Test the regularized feature
                    feature = Feature("REGULARIZED");
                    testCase.verifyEqual(feature.name, "regularized");
                    expectedOptions = struct('n_runs', int32(1), 'parameter', 1.0, 'order', 2);
                    testCase.verifyEqual(feature.options, expectedOptions);
                    testCase.verifyEqual(feature.modifier(x, f), f + norm(x), 'AbsTol', 1e-9);
                    
                    % Add custom options
                    feature = Feature("REGULARIZED", 'parameter', 2.0, 'order', 3);
                    testCase.verifyEqual(feature.name, "regularized");
                    expectedOptions = struct('n_runs', int32(1), 'parameter', 2.0, 'order', 3);
                    testCase.verifyEqual(feature.options, expectedOptions);
                    testCase.verifyEqual(feature.modifier(x, f), f + 2.0 * norm(x, 3), 'AbsTol', 1e-9);
                end
            end
        end

        function testNoisy(testCase)
            % Test that the noisy feature works correctly.
            nValues = [1, 10, 100];
            seedValues = [0, 1, 2];
        
            for n = nValues
                for seed = seedValues
                    % Generate random data
                    rng(seed);
                    x = randn(n, 1);
                    f = TestFeature.rosen(x);
                    
                    % Test the noisy feature
                    feature = Feature("NOISY");
                    testCase.verifyEqual(feature.name, "noisy");
                    testCase.verifyTrue(isfield(feature.options, 'distribution') && strcmp(feature.options.type, 'relative'));
                    testCase.verifyEqual(feature.options.n_runs, int32(10));
                    feature.modifier(x, f);
                    
                    % Add custom options
                    feature = Feature("NOISY", 'distribution', @(rng) 1.0, 'type', 'absolute', 'n_runs', 5);
                    testCase.verifyEqual(feature.name, "noisy");
                    testCase.verifyTrue(isfield(feature.options, 'distribution') && strcmp(feature.options.type, 'absolute'));
                    testCase.verifyEqual(feature.options.n_runs, int32(5));
                    testCase.verifyEqual(feature.modifier(x, f), f + 1.0, 'AbsTol', 1e-9);
                end
            end
        end

        function testTruncated(testCase)
            % Test that the truncated feature works correctly.
            nValues = [1, 10, 100];
            seedValues = [0, 1, 2];
        
            for n = nValues
                for seed = seedValues
                    % Generate random data
                    rng(seed);
                    x = randn(n, 1);
                    f = TestFeature.rosen(x);
                    
                    % Test the truncated feature
                    feature = Feature("TRUNCATED");
                    testCase.verifyEqual(feature.name, "truncated");
                    expectedOptions = struct('n_runs', int32(10), 'significant_digits', int32(6));
                    testCase.verifyEqual(feature.options, expectedOptions);
                    testCase.verifyTrue(all(abs(feature.modifier(x, f) - f) <= (1e-5 + 1e-5 * abs(f))));
                    testCase.verifyTrue(all(abs(feature.modifier(x, -f) + f) <= (1e-5 + 1e-5 * abs(-f))));
                    
                    % Add custom options
                    feature = Feature("TRUNCATED", 'significant_digits', int32(4));
                    testCase.verifyEqual(feature.name, "truncated");
                    expectedOptions = struct('n_runs', int32(10), 'significant_digits', int32(4));
                    testCase.verifyEqual(feature.options, expectedOptions);
                    testCase.verifyTrue(all(abs(feature.modifier(x, f) - f) <= (1e-3 + 1e-3 * abs(f))));
                end
            end

            % Test that the truncated feature works correctly with float significant digits.
            Feature("TRUNCATED", 'significant_digits', 2.0)
        end

        function testTough(testCase)
            % Test that the tough feature works correctly.
            nValues = [1, 10, 100];
            seedValues = [0, 1, 2];
        
            for n = nValues
                for seed = seedValues
                    % Generate random data
                    rng(seed);
                    x = randn(n, 1);
                    f = TestFeature.rosen(x);
                    
                    % Test the tough feature
                    feature = Feature("TOUGH");
                    testCase.verifyEqual(feature.name, "tough");
                    expectedOptions = struct('n_runs', int32(10), 'rate_error', 0.0, 'rate_nan', 0.05);
                    testCase.verifyEqual(feature.options, expectedOptions);

                    f_tough = feature.modifier(x, f);
                    testCase.verifyTrue(isequal(f_tough, f) || isnan(f_tough));
                    
                    % Add custom options
                    feature = Feature("TOUGH", 'rate_error', 0.5, 'rate_nan', 0.5);
                    testCase.verifyEqual(feature.name, "tough");
                    expectedOptions = struct('n_runs', int32(10), 'rate_error', 0.5, 'rate_nan', 0.5);
                    testCase.verifyEqual(feature.options, expectedOptions);
                    
                    try
                        f_tough = feature.modifier(x, f);
                        testCase.verifyTrue(isequal(f_tough, f) || isnan(f_tough));
                    catch ME
                        if ~strcmp(ME.identifier,"MATLAB:Feature:Tough")
                            rethrow(ME);
                        end
                    end
                end
            end
        end

        function testCustom(testCase)
            % Test that the custom feature works correctly.
            nValues = [1, 10, 100];
            seedValues = [0, 1, 2];
        
            for n = nValues
                for seed = seedValues
                    % Generate random data
                    rng(seed);
                    x = randn(n, 1);
                    f = TestFeature.rosen(x);
                    
                    % Test the custom feature
                    feature = Feature("custom", 'modifier', @(x, f, seed) f + 1.0);
                    testCase.verifyEqual(feature.name, "custom");
                    testCase.verifyTrue(isfield(feature.options, 'modifier'));
                    testCase.verifyTrue(isfield(feature.options, 'n_runs'));
                    testCase.verifyEqual(feature.options.n_runs, int32(1));
                    testCase.verifyEqual(feature.modifier(x, f), f + 1.0);
                end
            end
        end

        function testExceptions(testCase)
            % Test that the exceptions are thrown correctly.

            testCase.verifyError(@() Feature(), "MATLAB:Feature:MissingArguments");
            testCase.verifyError(@() Feature("unknown"), "MATLAB:Feature:UnknownFeature");
            testCase.verifyError(@() Feature("REGULARIZED", 'unknown', 1), "MATLAB:Feature:UnknownOption");
            testCase.verifyError(@() Feature("PLAIN", 'parameter', 1.0), "MATLAB:Feature:InvalidOptionForFeature");
            testCase.verifyError(@() Feature("CUSTOM", 'modifier', 1.0), "MATLAB:Feature:modifier_NotFunctionHandle");
            testCase.verifyError(@() Feature("NOISY", 'distribution', 1.0), "MATLAB:Feature:distribution_NotFunctionHandle");
            testCase.verifyError(@() Feature("REGULARIZED", 'order', '1.0'), "MATLAB:Feature:order_NotNumber");
            testCase.verifyError(@() Feature("REGULARIZED", 'parameter', '1.0'), "MATLAB:Feature:parameter_NotNonnegativeNumber");
            testCase.verifyError(@() Feature("REGULARIZED", 'parameter', -1.0), "MATLAB:Feature:parameter_NotNonnegativeNumber");
            testCase.verifyError(@() Feature("TOUGH", 'rate_error', '1.0'), "MATLAB:Feature:rate_error_NotBetween_0_1");
            testCase.verifyError(@() Feature("TOUGH", 'rate_error', -1.0), "MATLAB:Feature:rate_error_NotBetween_0_1");
            testCase.verifyError(@() Feature("TOUGH", 'rate_error', 2.0), "MATLAB:Feature:rate_error_NotBetween_0_1");
            testCase.verifyError(@() Feature("TOUGH", 'rate_nan', '1.0'), "MATLAB:Feature:rate_nan_NotBetween_0_1");
            testCase.verifyError(@() Feature("TOUGH", 'rate_nan', -1.0), "MATLAB:Feature:rate_nan_NotBetween_0_1");
            testCase.verifyError(@() Feature("TOUGH", 'rate_nan', 2.0), "MATLAB:Feature:rate_nan_NotBetween_0_1");
            testCase.verifyError(@() Feature("TRUNCATED", 'significant_digits', 2.5), "MATLAB:Feature:significant_digits_NotPositiveInteger");
            testCase.verifyError(@() Feature("NOISY", 'type', '+'), "MATLAB:Feature:type_InvalidInput");
            testCase.verifyError(@() Feature("CUSTOM"), "MATLAB:Feature:MissingModifier");
        end

    end
end