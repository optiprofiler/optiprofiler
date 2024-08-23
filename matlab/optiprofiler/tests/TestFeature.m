classdef TestFeature < matlab.unittest.TestCase
    methods (Static)

        function f = rosen(x)
            f = sum(100*(x(2:end) - x(1:end-1).^2).^2 + (1 - x(1:end-1)).^2);
        end

        function value = custom_distribution(randstream, n)
            if nargin < 2
                n = 1;
            end
            value = randstream.rand(n, 1);
        end

    end

    methods (Test)

        function testStruct(testCase)
            % Test the case where the second input is a struct.

            % Test the plain feature with a struct as the second input
            ft = Feature('plain', struct('n_runs', 3));
            randstream = ft.default_rng(1);
            x = randstream.randn(10, 1);
            f = TestFeature.rosen(x);
            testCase.verifyEqual(ft.name, 'plain');
            testCase.verifyEqual(ft.options, struct('n_runs', 3));
            testCase.verifyEqual(ft.modifier(x, f), f);
        end

        function testPlain(testCase)
            % Test that the plain feature works correctly.
            nValues = [1, 10, 100];
            seedValues = [0, 1, 2, nan];

            for n = nValues
                for seed = seedValues
                    % Test the default plain feature
                    ft = Feature('plain');
                    randstream = ft.default_rng(seed);
                    x = randstream.randn(n, 1);
                    f = @(x) TestFeature.rosen(x);
                    testCase.verifyEqual(ft.name, 'plain');
                    expectedOptions = struct('n_runs', 1);
                    testCase.verifyEqual(ft.options, expectedOptions);
                    g = ft.modifier(x, f);
                    testCase.verifyEqual(g(x), f(x));

                    % Add custom options
                    ft = Feature('plain', 'n_runs', 5);
                    testCase.verifyEqual(ft.options.n_runs, 5);
                end
            end
        end

        function testPerturbed_x0(testCase)
            % Test that the perturbed_x0 feature works correctly.
            nValues = [1, 10, 100];
            seedValues = [0, 1, 2];

            for n = nValues
                for seed = seedValues               
                    % Test the default perturbed_x0 feature
                    ft = Feature('perturbed_x0');
                    randstream = ft.default_rng(seed);
                    x = randstream.randn(n, 1);
                    f = @(x) TestFeature.rosen(x);
                    testCase.verifyEqual(ft.name, 'perturbed_x0');
                    testCase.verifyEqual(ft.options.n_runs, 10);
                    testCase.verifyEqual(ft.options.noise_level, 1e-3);
                    testCase.verifyEqual(ft.options.noise_type, 'relative');
                    
                    p = Problem(struct('fun', @(x) f(x), 'x0', x));
                    ftp = FeaturedProblem(p, ft, 500, seed);
                    if n == 1  % To cover one line in the code
                        testCase.verifyEqual(ftp.x0, max(1e-8, abs(p.x0)) .* sign(p.x0) .* (1 + 1e-3 .* ft.default_distribution(ft.default_rng(seed))));
                    else
                        testCase.verifyEqual(ftp.x0, max(1e-8, abs(p.x0)) .* sign(p.x0) .* (1 + 1e-3 .* ft.default_distribution(ft.default_rng(seed), n)));
                    end
                    testCase.verifyEqual(ftp.fun(ones(n, 1)), p.fun(ones(n, 1)));

                    % Add custom options
                    ft = Feature('perturbed_x0', 'n_runs', 5, 'noise_level', 1e-2, 'noise_type', 'absolute', 'distribution', @(x, y) TestFeature.custom_distribution(x, y));
                    testCase.verifyEqual(ft.options.n_runs, 5);
                    testCase.verifyEqual(ft.options.noise_level, 1e-2);
                    ftp = FeaturedProblem(p, ft, 500, seed);
                    testCase.verifyEqual(ftp.x0, p.x0 + 1e-2 * TestFeature.custom_distribution(ft.default_rng(seed), n));
                end
            end
        end

        function testNoisy(testCase)
            % Test that the noisy feature works correctly.
            nValues = [1, 10, 100];
            seedValues = [0, 1, 2];
        
            for n = nValues
                for seed = seedValues
                    % Test the default noisy feature
                    ft = Feature('noisy');
                    randstream = ft.default_rng(seed);
                    x = randstream.randn(n, 1);
                    f = @(x) TestFeature.rosen(x);
                    testCase.verifyEqual(ft.name, 'noisy');
                    testCase.verifyEqual(ft.options.n_runs, 10);
                    testCase.verifyEqual(ft.options.noise_level, 1e-3);
                    testCase.verifyEqual(ft.options.noise_type, 'relative');

                    p = Problem(struct('fun', @(x) f(x), 'x0', x));
                    ftp = FeaturedProblem(p, ft, 500, seed);
                    ftp.fun(x);

                    % Add custom options
                    ft = Feature('noisy', 'n_runs', 5, 'noise_level', 1e-2, 'noise_type', 'absolute', 'distribution', @(x, y) TestFeature.custom_distribution(x, y));
                    testCase.verifyEqual(ft.options.n_runs, 5);
                    testCase.verifyEqual(ft.options.noise_level, 1e-2);
                    testCase.verifyEqual(ft.options.noise_type, 'absolute');
                    ftp = FeaturedProblem(p, ft, 500, seed);
                    ftp.fun(x);
                end
            end
        end

        function testTruncated(testCase)
            % Test that the truncated feature works correctly.
            nValues = [1, 10, 100];
            seedValues = [0, 1, 2];
        
            for n = nValues
                for seed = seedValues
                    % Test the default truncated feature
                    ft = Feature('truncated');
                    randstream = ft.default_rng(seed);
                    x = randstream.randn(n, 1);
                    f = @(x) TestFeature.rosen(x);
                    testCase.verifyEqual(ft.name, 'truncated');
                    g = @(x) -f(x);
                    testCase.verifyEqual(ft.options.n_runs, 10);
                    testCase.verifyEqual(ft.options.significant_digits, 6);
                    testCase.verifyEqual(ft.options.perturbed_trailing_zeros, true);

                    p1 = Problem(struct('fun', @(x) f(x), 'x0', x));
                    p2 = Problem(struct('fun', @(x) g(x), 'x0', x));
                    ftp1 = FeaturedProblem(p1, ft, 500, seed);
                    ftp1.fun(x);
                    ftp2 = FeaturedProblem(p2, ft, 500, seed);
                    ftp2.fun(x);
                    
                    % Add custom options
                    ft = Feature('truncated', 'n_runs', 5, 'significant_digits', 2, 'perturbed_trailing_zeros', false);
                    testCase.verifyEqual(ft.options.n_runs, 5);
                    testCase.verifyEqual(ft.options.significant_digits, 2);
                    testCase.verifyEqual(ft.options.perturbed_trailing_zeros, false);
                    ftp = FeaturedProblem(p1, ft, 500, seed);
                    ftp.fun(x);

                    ft = Feature('truncated', 'significant_digits', 2, 'perturbed_trailing_zeros', false);
                    testCase.verifyEqual(ft.options.n_runs, 1);
                end
            end
        end

        function testPermuted(testCase)
            % Test that the permuted feature works correctly.
            nValues = [1, 10, 100];
            seedValues = [0, 1, 2];
        
            for n = nValues
                for seed = seedValues
                    % Test the default permuted feature
                    ft = Feature('permuted');
                    randstream = ft.default_rng(seed);
                    x = randstream.randn(n, 1);
                    f = @(x) TestFeature.rosen(x);
                    testCase.verifyEqual(ft.name, 'permuted');
                    testCase.verifyEqual(ft.options.n_runs, 10);
                    
                    p = Problem(struct('fun', @(x) f(x), 'x0', x));
                    ftp = FeaturedProblem(p, ft, 500, seed);
                    testCase.verifyEqual(ftp.fun(ftp.x0), p.fun(p.x0));

                    % Add custom options
                    ft = Feature('permuted', 'n_runs', 5);
                    testCase.verifyEqual(ft.options.n_runs, 5);
                    ftp = FeaturedProblem(p, ft, 500, seed);
                    testCase.verifyEqual(ftp.fun(ftp.x0), p.fun(p.x0));
                end
            end
        end

        function testLinealy_transformed(testCase)
            % Test that the linearly_transformed feature works correctly.
            nValues = [1, 10, 100];
            seedValues = [0, 1, 2];
        
            for n = nValues
                for seed = seedValues
                    % Test the default linearly_transformed feature
                    ft = Feature('linearly_transformed');
                    randstream = ft.default_rng(seed);
                    x = randstream.randn(n, 1);
                    f = @(x) TestFeature.rosen(x);
                    testCase.verifyEqual(ft.name, 'linearly_transformed');
                    testCase.verifyEqual(ft.options.n_runs, 10);
                    testCase.verifyEqual(ft.options.rotated, true);
                    testCase.verifyEqual(ft.options.condition_number, 1);
                    
                    p = Problem(struct('fun', @(x) f(x), 'x0', x));
                    ftp = FeaturedProblem(p, ft, 500, seed);
                    ftp.fun(x);
                    
                    % Add custom options
                    ft = Feature('linearly_transformed', 'n_runs', 5, 'rotated', false, 'condition_number', 10);
                    testCase.verifyEqual(ft.options.n_runs, 5);
                    testCase.verifyEqual(ft.options.rotated, false);
                    testCase.verifyEqual(ft.options.condition_number, 10);
                    ftp = FeaturedProblem(p, ft, 500, seed);
                    ftp.fun(x);
                end
            end
        end

        function testRandom_nan(testCase)
            % Test that the random_nan feature works correctly.
            nValues = [1, 10, 100];
            seedValues = [0, 1, 2];
        
            for n = nValues
                for seed = seedValues
                    % Test the default random_nan feature
                    ft = Feature('random_nan');
                    randstream = ft.default_rng(seed);
                    x = randstream.randn(n, 1);
                    f = @(x) TestFeature.rosen(x);
                    testCase.verifyEqual(ft.name, 'random_nan');
                    testCase.verifyEqual(ft.options.n_runs, 10);
                    testCase.verifyEqual(ft.options.rate_nan, 0.05);
                    
                    p = Problem(struct('fun', @(x) f(x), 'x0', x));
                    ftp = FeaturedProblem(p, ft, 500, seed);
                    ftp.fun(x);
                    
                    % Add custom options
                    ft = Feature('random_nan', 'n_runs', 5, 'rate_nan', 0.2);
                    testCase.verifyEqual(ft.options.n_runs, 5);
                    testCase.verifyEqual(ft.options.rate_nan, 0.2);
                    ftp = FeaturedProblem(p, ft, 500, seed);
                    ftp.fun(x);
                end
            end
        end

        function testUnrelaxable_constraints(testCase)
            % Test that the unrelaxable_constraints feature works correctly.
            nValues = [1, 10, 100];
            seedValues = [0, 1, 2];

            for n = nValues
                for seed = seedValues
                    % Test the default unrelaxable_constraints feature
                    ft = Feature('unrelaxable_constraints');
                    randstream = ft.default_rng(seed);
                    x = randstream.randn(n, 1);
                    f = @(x) TestFeature.rosen(x);
                    testCase.verifyEqual(ft.name, 'unrelaxable_constraints');
                    testCase.verifyEqual(ft.options.n_runs, 1);
                    testCase.verifyEqual(ft.options.unrelaxable_bounds, false);
                    testCase.verifyEqual(ft.options.unrelaxable_linear_constraints, false);
                    testCase.verifyEqual(ft.options.unrelaxable_nonlinear_constraints, false);
                    
                    p = Problem(struct('fun', @(x) f(x), 'x0', x));
                    ftp = FeaturedProblem(p, ft, 500, seed);
                    ftp.fun(x);
                    
                    % Add custom options
                    ft = Feature('unrelaxable_constraints', 'n_runs', 5, 'unrelaxable_bounds', true, 'unrelaxable_linear_constraints', true, 'unrelaxable_nonlinear_constraints', 1);
                    testCase.verifyEqual(ft.options.n_runs, 5);
                    testCase.verifyEqual(ft.options.unrelaxable_bounds, true);
                    testCase.verifyEqual(ft.options.unrelaxable_linear_constraints, true);
                    testCase.verifyEqual(ft.options.unrelaxable_nonlinear_constraints, 1);
                    p1 = Problem(struct('fun', @(x) f(x), 'x0', x, 'xl', x + 1));
                    p2 = Problem(struct('fun', @(x) f(x), 'x0', x, 'aub', eye(n), 'bub', x - 2));
                    p3 = Problem(struct('fun', @(x) f(x), 'x0', x, 'cub', @(x) f(x) + 1e8));
                    ftp1 = FeaturedProblem(p1, ft, 500, seed);
                    ftp1.fun(x);
                    ftp2 = FeaturedProblem(p2, ft, 500, seed);
                    ftp2.fun(x);
                    ftp3 = FeaturedProblem(p3, ft, 500, seed);
                    ftp3.fun(x);
                end
            end
        end

        function testCustom(testCase)
            % Test that the custom feature works correctly.
            nValues = [1, 10, 100];
            seedValues = [0, 1, 2];

            % Test the simplest custom feature
            ft = Feature('custom', 'modifier', @(x, f) f + 1.0);
            testCase.verifyEqual(ft.options.n_runs, 1);
        
            for n = nValues
                for seed = seedValues
                    % Test the custom feature
                    ft = Feature('custom', 'n_runs', 5, 'modifier', @(x, f, seed) f + 1.0);
                    randstream = ft.default_rng(seed);
                    x = randstream.randn(n, 1);
                    f = @(x) TestFeature.rosen(x);
                    testCase.verifyEqual(ft.name, 'custom');
                    testCase.verifyEqual(ft.options.n_runs, 5);

                    p = Problem(struct('fun', @(x) f(x), 'x0', x));
                    ftp = FeaturedProblem(p, ft, 500, seed);
                    testCase.verifyEqual(ftp.fun(x), f(x) + 1.0);
                end
            end
        end

        function testIsStochastic(testCase)
            % Test that the function IsStochastic works correctly.

            feature_plain = Feature('plain');
            testCase.verifyFalse(feature_plain.isStochastic());
            feature_perturbed_x0 = Feature('noisy');
            testCase.verifyTrue(feature_perturbed_x0.isStochastic());
        end

        function testExceptions(testCase)
            % Test that the exceptions are thrown correctly.

            testCase.verifyError(@() Feature(1), "MATLAB:Feature:FeaturenameNotString");
            testCase.verifyError(@() Feature('custom', 'n_runs', 1, 'distribution'), "MATLAB:Feature:InvalidNumberOfArguments");
            testCase.verifyError(@() Feature(''), "MATLAB:Feature:UnknownFeature");
            testCase.verifyError(@() Feature('unknown'), "MATLAB:Feature:UnknownFeature");
            testCase.verifyError(@() Feature('plain', 'unknown_option', 1.0), "MATLAB:Feature:UnknownOption");
            testCase.verifyError(@() Feature('plain', 'noise_level', 1.0), "MATLAB:Feature:InvalidOptionForFeature");
            testCase.verifyError(@() Feature('custom', 'modifier', 1.0), "MATLAB:Feature:modifier_NotFunctionHandle");
            testCase.verifyError(@() Feature('noisy', 'distribution', 1.0), "MATLAB:Feature:distribution_NotFunctionHandle");
            testCase.verifyError(@() Feature('random_nan', 'rate_nan', '1.0'), "MATLAB:Feature:rate_nan_NotBetween_0_1");
            testCase.verifyError(@() Feature('random_nan', 'rate_nan', -1.0), "MATLAB:Feature:rate_nan_NotBetween_0_1");
            testCase.verifyError(@() Feature('random_nan', 'rate_nan', 2.0), "MATLAB:Feature:rate_nan_NotBetween_0_1");
            testCase.verifyError(@() Feature('truncated', 'significant_digits', 2.5), "MATLAB:Feature:significant_digits_NotPositiveInteger");
            testCase.verifyError(@() Feature('perturbed_x0', 'noise_level', -1), "MATLAB:Feature:noise_level_NotPositive");
            testCase.verifyError(@() Feature('truncated', 'perturbed_trailing_zeros', -1), "MATLAB:Feature:perturbed_trailing_zeros_NotLogical");
            testCase.verifyError(@() Feature('linearly_transformed', 'rotated', 2), "MATLAB:Feature:rotated_NotLogical");
            testCase.verifyError(@() Feature('linearly_transformed', 'condition_number', -1), "MATLAB:Feature:condition_number_InvalidInput");
            testCase.verifyError(@() Feature('noisy', 'n_runs', 2.5), "MATLAB:Feature:n_runs_NotPositiveInteger");
            testCase.verifyError(@() Feature('noisy', 'noise_type', '+'), "MATLAB:Feature:noise_type_InvalidInput");
            testCase.verifyError(@() Feature('unrelaxable_constraints', 'unrelaxable_bounds', -1), "MATLAB:Feature:unrelaxable_bounds_NotLogical");
            testCase.verifyError(@() Feature('unrelaxable_constraints', 'unrelaxable_linear_constraints', -1), "MATLAB:Feature:unrelaxable_linear_constraints_NotLogical");
            testCase.verifyError(@() Feature('unrelaxable_constraints', 'unrelaxable_nonlinear_constraints', -1), "MATLAB:Feature:unrelaxable_nonlinear_constraints_NotLogical");
            testCase.verifyError(@() Feature("CUSTOM"), "MATLAB:Feature:MissingModifier");
            ft = Feature('noisy');
            testCase.verifyError(@() ft.default_rng(1 + 1j), "MATLAB:Feature:SeedNotEvenReal");
        end

    end
end