classdef TestFeature < matlab.unittest.TestCase
    methods (Static)

        function f = rosen(x)
            f = sum(100*(x(2:end) - x(1:end-1).^2).^2 + (1 - x(1:end-1)).^2);
        end

    end

    methods (Test)

        function testStruct(testCase)
            % Test the case where the second input is a struct.

            % Generate random data
            rng(1);
            x = randn(10, 1);
            f = TestFeature.rosen(x);

            % Test the plain feature
            feature = Feature("PLAIN", struct('n_runs', 3));
            testCase.verifyEqual(feature.name, "plain");
            testCase.verifyEqual(feature.options, struct('n_runs', int32(3)));
            testCase.verifyEqual(feature.modifier(x, f), f);
        end

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

        function testPerturbed_x0(testCase)
            % Test that the perturbed_x0 feature works correctly.
            nValues = [1, 10, 100];
            seedValues = [0, 1, 2];

            for n = nValues
                for seed = seedValues
                    % Generate random data
                    rng(seed);
                    
                    % Test the perturbed_x0 feature
                    feature = Feature("PERTURBED_X0");
                    testCase.verifyEqual(feature.name, "perturbed_x0");
                    testCase.verifyTrue(isfield(feature.options, 'distribution'));
                    testCase.verifyEqual(feature.options.n_runs, int32(10));
                    
                    % Add custom options
                    feature = Feature("PERTURBED_X0", 'distribution', @(rng) 1.0, 'n_runs', 5);
                    testCase.verifyEqual(feature.name, "perturbed_x0");
                    testCase.verifyTrue(isfield(feature.options, 'distribution'));
                    testCase.verifyEqual(feature.options.n_runs, int32(5));
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

                    feature = Feature("TRUNCATED", 'significant_digits', int32(2), 'perturbed_trailing_zeros', false);
                    testCase.verifyEqual(feature.name, "truncated");
                    expectedOptions = struct('n_runs', int32(1), 'significant_digits', int32(2), 'perturbed_trailing_zeros', false);
                    testCase.verifyEqual(feature.options, expectedOptions);
                end
            end

            % Test that the truncated feature works correctly with float significant digits.
            Feature("TRUNCATED", 'significant_digits', 2.0)
        end

        function testRotated(testCase)
            % Test that the rotated feature works correctly.
            nValues = [1, 10, 100];
            seedValues = [0, 1, 2];
        
            for n = nValues
                for seed = seedValues
                    % Generate random data
                    rng(seed);
                    x = randn(n, 1);
                    f = TestFeature.rosen(x);
                    
                    % Test the rotated feature
                    feature = Feature("ROTATED");
                    testCase.verifyEqual(feature.name, "rotated");
                    expectedOptions = struct('n_runs', int32(10));
                    testCase.verifyEqual(feature.options, expectedOptions);
                end
            end
        end

        function testPermuted(testCase)
            % Test that the permuted feature works correctly.
            nValues = [1, 10, 100];
            seedValues = [0, 1, 2];
        
            for n = nValues
                for seed = seedValues
                    % Generate random data
                    rng(seed);
                    x = randn(n, 1);
                    f = TestFeature.rosen(x);
                    
                    % Test the permuted feature
                    feature = Feature("PERMUTED");
                    testCase.verifyEqual(feature.name, "permuted");
                    expectedOptions = struct('n_runs', int32(10));
                    testCase.verifyEqual(feature.options, expectedOptions);
                end
            end
        end

        function testRandom_nan(testCase)
            % Test that the random_nan feature works correctly.
            nValues = [1, 10, 100];
            seedValues = [0, 1, 2];
        
            for n = nValues
                for seed = seedValues
                    % Generate random data
                    rng(seed);
                    x = randn(n, 1);
                    f = TestFeature.rosen(x);
                    
                    % Test the random_nan feature
                    feature = Feature("RANDOM_NAN");
                    testCase.verifyEqual(feature.name, "random_nan");
                    expectedOptions = struct('n_runs', int32(10), 'rate_nan', 0.999);
                    testCase.verifyEqual(feature.options, expectedOptions);

                    f_tough = feature.modifier(x, f);
                    testCase.verifyTrue(isequal(f_tough, f) || isnan(f_tough));

                    % Add custom options
                    feature = Feature("RANDOM_NAN", 'rate_nan', 0.1);
                    testCase.verifyEqual(feature.name, "random_nan");
                    expectedOptions = struct('n_runs', int32(10), 'rate_nan', 0.1);
                    testCase.verifyEqual(feature.options, expectedOptions);
                end
            end
        end

        function testUnrelaxable_constraints(testCase)
            % Test that the unrelaxable_constraints feature works correctly.
            nValues = [1, 10, 100];
            seedValues = [0, 1, 2];

            for n = nValues
                for seed = seedValues
                    % Generate random data
                    rng(seed);
                    x = randn(n, 1);
                    f = TestFeature.rosen(x);
                    
                    % Test the unrelaxable_constraints feature
                    feature = Feature("UNRELAXABLE_CONSTRAINTS");
                    testCase.verifyEqual(feature.name, "unrelaxable_constraints");
                    expectedOptions = struct('n_runs', int32(1), 'unrelaxable_bounds', false, 'unrelaxable_linear_constraints', false, 'unrelaxable_nonlinear_constraints', false);
                    testCase.verifyEqual(feature.options, expectedOptions);

                    feature = Feature("UNRELAXABLE_CONSTRAINTS", 'unrelaxable_bounds', true, 'unrelaxable_linear_constraints', true, 'unrelaxable_nonlinear_constraints', true);
                    testCase.verifyEqual(feature.name, "unrelaxable_constraints");
                    expectedOptions = struct('n_runs', int32(1), 'unrelaxable_bounds', true, 'unrelaxable_linear_constraints', true, 'unrelaxable_nonlinear_constraints', true);
                    testCase.verifyEqual(feature.options, expectedOptions);

                    f_unrelaxable_bounds = feature.modifier(x, f, seed, 1);
                    testCase.verifyTrue(isinf(f_unrelaxable_bounds));
                    f_unrelaxable_linear_constraints = feature.modifier(x, f, seed, 0, 1);
                    testCase.verifyTrue(isinf(f_unrelaxable_linear_constraints));
                    f_unrelaxable_nonlinear_constraints = feature.modifier(x, f, seed, 0, 0, 1);
                    testCase.verifyTrue(isinf(f_unrelaxable_nonlinear_constraints));
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

        function testIsStochastic(testCase)
            % Test that the function IsStochastic works correctly.

            feature_plain = Feature("PLAIN");
            testCase.verifyFalse(feature_plain.isStochastic());
            feature_perturbed_x0 = Feature("PERTURBED_X0");
            testCase.verifyTrue(feature_perturbed_x0.isStochastic());
        end

        function testExceptions(testCase)
            % Test that the exceptions are thrown correctly.

            testCase.verifyError(@() Feature(1), "MATLAB:Feature:FeaturenameNotString");
            testCase.verifyError(@() Feature("CUSTOM", 'n_runs', 1, 'distribution'), "MATLAB:Feature:InvalidNumberOfArguments");
            testCase.verifyError(@() Feature(""), "MATLAB:Feature:MissingArguments");
            testCase.verifyError(@() Feature("unknown"), "MATLAB:Feature:UnknownFeature");
            testCase.verifyError(@() Feature("PLAIN", 'parameter', 1.0), "MATLAB:Feature:InvalidOptionForFeature");
            testCase.verifyError(@() Feature("CUSTOM", 'modifier', 1.0), "MATLAB:Feature:modifier_NotFunctionHandle");
            testCase.verifyError(@() Feature("NOISY", 'distribution', 1.0), "MATLAB:Feature:distribution_NotFunctionHandle");
            testCase.verifyError(@() Feature("RANDOM_NAN", 'rate_nan', '1.0'), "MATLAB:Feature:rate_nan_NotBetween_0_1");
            testCase.verifyError(@() Feature("RANDOM_NAN", 'rate_nan', -1.0), "MATLAB:Feature:rate_nan_NotBetween_0_1");
            testCase.verifyError(@() Feature("RANDOM_NAN", 'rate_nan', 2.0), "MATLAB:Feature:rate_nan_NotBetween_0_1");
            testCase.verifyError(@() Feature("TRUNCATED", 'significant_digits', 2.5), "MATLAB:Feature:significant_digits_NotPositiveInteger");
            testCase.verifyError(@() Feature("TRUNCATED", 'perturbed_trailing_zeros', 1), "MATLAB:Feature:perturbed_trailing_zeros_NotLogical");
            testCase.verifyError(@() Feature("NOISY", 'n_runs', 2.5), "MATLAB:Feature:n_runs_NotPositiveInteger");
            testCase.verifyError(@() Feature("NOISY", 'type', '+'), "MATLAB:Feature:type_InvalidInput");
            testCase.verifyError(@() Feature("UNRELAXABLE_CONSTRAINTS", 'unrelaxable_bounds', 1), "MATLAB:Feature:unrelaxable_bounds_NotLogical");
            testCase.verifyError(@() Feature("UNRELAXABLE_CONSTRAINTS", 'unrelaxable_linear_constraints', 1), "MATLAB:Feature:unrelaxable_linear_constraints_NotLogical");
            testCase.verifyError(@() Feature("UNRELAXABLE_CONSTRAINTS", 'unrelaxable_nonlinear_constraints', 1), "MATLAB:Feature:unrelaxable_nonlinear_constraints_NotLogical");
            testCase.verifyError(@() Feature("CUSTOM"), "MATLAB:Feature:MissingModifier");
        end

    end
end