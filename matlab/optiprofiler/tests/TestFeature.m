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

        function x0 = custom_mod_x0(rand_stream, p)
            x0 = p.x0 + rand_stream.rand(size(p.x0));
        end

        function [A, b, inv] = custom_mod_affine(rand_stream, p)
            n = size(p.x0, 1);
            A = diag(1 + rand_stream.rand(n, 1));
            b = rand_stream.rand(n, 1);
            inv = 0;
        end

        function [xl, xu] = custom_mod_bounds(rand_stream, p)
            n = size(p.x0, 1);
            xl = -rand_stream.rand(n, 1);
            xu = rand_stream.rand(n, 1);
        end

        function [aub, bub] = custom_mod_linear_ub(rand_stream, p)
            n = size(p.x0, 1);
            aub = rand_stream.rand(n, 1);
            bub = rand_stream.rand(1);
        end

        function [aeq, beq] = custom_mod_linear_eq(rand_stream, p)
            n = size(p.x0, 1);
            aeq = rand_stream.rand(n, 1);
            beq = rand_stream.rand(1);
        end

        function f = custom_mod_fun(x, rand_stream, p)
            f = p.fun(x) + rand_stream.rand(1);
        end

        function cub = custom_mod_cub(x, rand_stream, p)
            m = size(p.cub(x), 1);
            cub = p.cub(x) + rand_stream.rand(m, 1);
        end

        function ceq = custom_mod_ceq(x, rand_stream, p)
            m = size(p.ceq(x), 1);
            ceq = p.ceq(x) + rand_stream.rand(m, 1);
        end

    end

    methods (Test)

        function testStruct(testCase)
            % Test the case where the second input is a struct.

            ft = Feature('plain', struct('n_runs', 3));
            testCase.verifyEqual(ft.name, 'plain');
            testCase.verifyEqual(ft.options, struct('n_runs', 3));
        end

        function testPlain(testCase)
            % Test the 'plain' feature.

            ft = Feature('plain');
            testCase.verifyEqual(ft.name, 'plain');
            testCase.verifyEqual(ft.options, struct('n_runs', 1));

            options.n_runs = 5;
            ft = Feature('plain', options);
            testCase.verifyEqual(ft.name, 'plain');
            testCase.verifyEqual(ft.options, struct('n_runs', 5));
        end

        function testPerturbed_x0(testCase)
            % Test the 'perturbed_x0' feature.

            ft = Feature('perturbed_x0');
            testCase.verifyEqual(ft.name, 'perturbed_x0');
            testCase.verifyEqual(ft.options.n_runs, 10);
            testCase.verifyEqual(ft.options.noise_level, 1e-3);

            options.n_runs = 5;
            options.noise_level = 1e-2;
            ft = Feature('perturbed_x0', options);
            testCase.verifyEqual(ft.name, 'perturbed_x0');
            testCase.verifyEqual(ft.options.n_runs, 5);
            testCase.verifyEqual(ft.options.noise_level, 1e-2);
        end

        function testNoisy(testCase)
            % Test the 'noisy' feature.
        
            ft = Feature('noisy');
            testCase.verifyEqual(ft.name, 'noisy');
            testCase.verifyEqual(ft.options.n_runs, 10);
            testCase.verifyEqual(ft.options.noise_type, 'mixed');
            testCase.verifyEqual(ft.options.noise_level, 1e-3);

            options.n_runs = 5;
            options.noise_type = 'absolute';
            options.noise_level = 1e-2;
            ft = Feature('noisy', options);
            testCase.verifyEqual(ft.name, 'noisy');
            testCase.verifyEqual(ft.options.n_runs, 5);
            testCase.verifyEqual(ft.options.noise_type, 'absolute');
            testCase.verifyEqual(ft.options.noise_level, 1e-2);
        end

        function testTruncated(testCase)
            % Test the 'truncated' feature.

            ft = Feature('truncated');
            testCase.verifyEqual(ft.name, 'truncated');
            testCase.verifyEqual(ft.options.n_runs, 10);
            testCase.verifyEqual(ft.options.significant_digits, 6);
            testCase.verifyEqual(ft.options.perturbed_trailing_zeros, true);

            options.n_runs = 5;
            options.significant_digits = 2;
            options.perturbed_trailing_zeros = false;
            ft = Feature('truncated', options);
            testCase.verifyEqual(ft.name, 'truncated');
            testCase.verifyEqual(ft.options.n_runs, 5);
            testCase.verifyEqual(ft.options.significant_digits, 2);
            testCase.verifyEqual(ft.options.perturbed_trailing_zeros, false);
        end

        function testPermuted(testCase)
            % Test the 'permuted' feature.

            ft = Feature('permuted');
            testCase.verifyEqual(ft.name, 'permuted');
            testCase.verifyEqual(ft.options.n_runs, 10);
            
            options.n_runs = 5;
            ft = Feature('permuted', options);
            testCase.verifyEqual(ft.name, 'permuted');
            testCase.verifyEqual(ft.options.n_runs, 5);
        end

        function testLinealy_transformed(testCase)
            % Test the 'linearly_transformed' feature.

            ft = Feature('linearly_transformed');
            testCase.verifyEqual(ft.name, 'linearly_transformed');
            testCase.verifyEqual(ft.options.n_runs, 10);
            testCase.verifyEqual(ft.options.rotated, true);
            testCase.verifyEqual(ft.options.condition_factor, 0);

            options.n_runs = 5;
            options.rotated = false;
            options.condition_factor = 1e3;
            ft = Feature('linearly_transformed', options);
            testCase.verifyEqual(ft.name, 'linearly_transformed');
            testCase.verifyEqual(ft.options.n_runs, 5);
            testCase.verifyEqual(ft.options.rotated, false);
            testCase.verifyEqual(ft.options.condition_factor, 1e3);
        end

        function testRandom_nan(testCase)
            % Test the 'random_nan' feature.
            
            ft = Feature('random_nan');
            testCase.verifyEqual(ft.name, 'random_nan');
            testCase.verifyEqual(ft.options.n_runs, 10);
            testCase.verifyEqual(ft.options.rate_nan, 0.05);

            options.n_runs = 5;
            options.rate_nan = 0.1;
            ft = Feature('random_nan', options);
            testCase.verifyEqual(ft.name, 'random_nan');
            testCase.verifyEqual(ft.options.n_runs, 5);
            testCase.verifyEqual(ft.options.rate_nan, 0.1);
        end

        function testUnrelaxable_constraints(testCase)
            % Test the 'unrelaxable_constraints' feature.
            
            ft = Feature('unrelaxable_constraints');
            testCase.verifyEqual(ft.name, 'unrelaxable_constraints');
            testCase.verifyEqual(ft.options.n_runs, 1);
            testCase.verifyEqual(ft.options.unrelaxable_bounds, false);
            testCase.verifyEqual(ft.options.unrelaxable_linear_constraints, false);
            testCase.verifyEqual(ft.options.unrelaxable_nonlinear_constraints, false);

            options.n_runs = 5;
            options.unrelaxable_bounds = true;
            options.unrelaxable_linear_constraints = true;
            options.unrelaxable_nonlinear_constraints = true;
            ft = Feature('unrelaxable_constraints', options);
            testCase.verifyEqual(ft.name, 'unrelaxable_constraints');
            testCase.verifyEqual(ft.options.n_runs, 5);
            testCase.verifyEqual(ft.options.unrelaxable_bounds, true);
            testCase.verifyEqual(ft.options.unrelaxable_linear_constraints, true);
            testCase.verifyEqual(ft.options.unrelaxable_nonlinear_constraints, true);
        end

        function testNonquantifiable_constraints(testCase)
            % Test the 'nonquantifiable_constraints' feature.
            
            ft = Feature('nonquantifiable_constraints');
            testCase.verifyEqual(ft.name, 'nonquantifiable_constraints');
            testCase.verifyEqual(ft.options.n_runs, 1);
            
            options.n_runs = 5;
            ft = Feature('nonquantifiable_constraints', options);
            testCase.verifyEqual(ft.name, 'nonquantifiable_constraints');
            testCase.verifyEqual(ft.options.n_runs, 5);
        end

        function testQuantized(testCase)
            % Test the 'quantized' feature.
            
            ft = Feature('quantized');
            testCase.verifyEqual(ft.name, 'quantized');
            testCase.verifyEqual(ft.options.n_runs, 1);
            testCase.verifyEqual(ft.options.mesh_size, 1e-3);
            testCase.verifyEqual(ft.options.is_truth, true);
            
            options.n_runs = 5;
            options.mesh_size = 1e-2;
            options.is_truth = false;
            ft = Feature('quantized', options);
            testCase.verifyEqual(ft.name, 'quantized');
            testCase.verifyEqual(ft.options.n_runs, 5);
            testCase.verifyEqual(ft.options.mesh_size, 1e-2);
            testCase.verifyEqual(ft.options.is_truth, false);
        end

        function testCustom(testCase)
            % Test the 'custom' feature.
            
            ft = Feature('custom');
            testCase.verifyEqual(ft.name, 'custom');
            testCase.verifyEqual(ft.options.n_runs, 1);

            options.n_runs = 5;
            options.mod_x0 = @TestFeature.custom_mod_x0;
            options.mod_affine = @TestFeature.custom_mod_affine;
            options.mod_bounds = @TestFeature.custom_mod_bounds;
            options.mod_linear_ub = @TestFeature.custom_mod_linear_ub;
            options.mod_linear_eq = @TestFeature.custom_mod_linear_eq;
            options.mod_fun = @TestFeature.custom_mod_fun;
            options.mod_cub = @TestFeature.custom_mod_cub;
            options.mod_ceq = @TestFeature.custom_mod_ceq;
            ft = Feature('custom', options);
            testCase.verifyEqual(ft.name, 'custom');
            testCase.verifyEqual(ft.options.n_runs, 5);
            testCase.verifyEqual(ft.options.mod_x0, @TestFeature.custom_mod_x0);
            testCase.verifyEqual(ft.options.mod_affine, @TestFeature.custom_mod_affine);
            testCase.verifyEqual(ft.options.mod_bounds, @TestFeature.custom_mod_bounds);
            testCase.verifyEqual(ft.options.mod_linear_ub, @TestFeature.custom_mod_linear_ub);
            testCase.verifyEqual(ft.options.mod_linear_eq, @TestFeature.custom_mod_linear_eq);
            testCase.verifyEqual(ft.options.mod_fun, @TestFeature.custom_mod_fun);
            testCase.verifyEqual(ft.options.mod_cub, @TestFeature.custom_mod_cub);
            testCase.verifyEqual(ft.options.mod_ceq, @TestFeature.custom_mod_ceq);
        end
        
        function testErrors(testCase)
            % Test whether the function throws errors as expected.

            testCase.verifyError(@() Feature(1), "MATLAB:Feature:FeaturenameNotString")

            testCase.verifyError(@() Feature('plain', 'n_runs'), "MATLAB:Feature:InvalidNumberOfArguments")

            testCase.verifyError(@() Feature('unknown'), "MATLAB:Feature:UnknownFeature")

            testCase.verifyError(@() Feature('plain', struct('n_runs', 1, 'unknown', 1)), "MATLAB:Feature:UnknownOption")

            testCase.verifyError(@() Feature('plain', struct('n_runs', 1, 'noise_level', '1e-3')), "MATLAB:Feature:InvalidOptionForFeature")

            testCase.verifyError(@() Feature('plain', struct('n_runs', 1.1)), "MATLAB:Feature:n_runs_NotPositiveInteger")

            testCase.verifyError(@() Feature('noisy', struct('distribution', 'normal')), "MATLAB:Feature:distribution_NotFunctionHandle")

            testCase.verifyError(@() Feature('random_nan', struct('rate_nan', 1.1)), "MATLAB:Feature:rate_nan_NotBetween_0_1")

            testCase.verifyError(@() Feature('truncated', struct('significant_digits', 0)), "MATLAB:Feature:significant_digits_NotPositiveInteger")

            testCase.verifyError(@() Feature('noisy', struct('noise_level', 'unknown')), "MATLAB:Feature:noise_level_NotPositive")

            testCase.verifyError(@() Feature('noisy', struct('noise_type', 'unknown')), "MATLAB:Feature:noise_type_InvalidInput")

            testCase.verifyError(@() Feature('truncated', struct('perturbed_trailing_zeros', 'unknown')), "MATLAB:Feature:perturbed_trailing_zeros_NotLogical")

            testCase.verifyError(@() Feature('linearly_transformed', struct('rotated', 'unknown')), "MATLAB:Feature:rotated_NotLogical")

            testCase.verifyError(@() Feature('linearly_transformed', struct('condition_factor', 'unknown')), "MATLAB:Feature:condition_factor_InvalidInput")

            testCase.verifyError(@() Feature('unrelaxable_constraints', struct('unrelaxable_bounds', 'unknown')), "MATLAB:Feature:unrelaxable_bounds_NotLogical")

            testCase.verifyError(@() Feature('unrelaxable_constraints', struct('unrelaxable_linear_constraints', 'unknown')), "MATLAB:Feature:unrelaxable_linear_constraints_NotLogical")

            testCase.verifyError(@() Feature('unrelaxable_constraints', struct('unrelaxable_nonlinear_constraints', 'unknown')), "MATLAB:Feature:unrelaxable_nonlinear_constraints_NotLogical")

            testCase.verifyError(@() Feature('quantized', struct('mesh_size', 'unknown')), "MATLAB:Feature:mesh_size_NotPositive")

            testCase.verifyError(@() Feature('quantized', struct('is_truth', 'unknown')), "MATLAB:Feature:is_truth_NotLogical")

            testCase.verifyError(@() Feature('custom', struct('mod_x0', 'unknown')), "MATLAB:Feature:mod_x0_NotFunctionHandle")

            testCase.verifyError(@() Feature('custom', struct('mod_bounds', 'unknown')), "MATLAB:Feature:mod_bounds_NotFunctionHandle")

            testCase.verifyError(@() Feature('custom', struct('mod_linear_ub', 'unknown')), "MATLAB:Feature:mod_linear_ub_NotFunctionHandle")

            testCase.verifyError(@() Feature('custom', struct('mod_linear_eq', 'unknown')), "MATLAB:Feature:mod_linear_eq_NotFunctionHandle")

            testCase.verifyError(@() Feature('custom', struct('mod_affine', 'unknown')), "MATLAB:Feature:mod_affine_NotFunctionHandle")

            testCase.verifyError(@() Feature('custom', struct('mod_fun', 'unknown')), "MATLAB:Feature:mod_fun_NotFunctionHandle")

            testCase.verifyError(@() Feature('custom', struct('mod_cub', 'unknown')), "MATLAB:Feature:mod_cub_NotFunctionHandle")

            testCase.verifyError(@() Feature('custom', struct('mod_ceq', 'unknown')), "MATLAB:Feature:mod_ceq_NotFunctionHandle")

            options.mod_affine = @TestFeature.custom_mod_affine;
            ft = Feature('custom', options);
            p = Problem(struct('fun', @TestFeature.rosen, 'x0', [-1; -1; -1]));
            testCase.verifyError(@() FeaturedProblem(p, ft, 500, 1), "MATLAB:Feature:AffineTransformationNotInvertible")

            ft = Feature('plain');
            testCase.verifyError(@() ft.default_rng('a'), "MATLAB:Feature:SeedNotEvenReal")
        end

    end
end