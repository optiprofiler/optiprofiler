classdef TestBenchmark < matlab.unittest.TestCase
    methods (Test)
        
        function testWithValidInput(testCase)

            solvers = {@fmincon_test1, @fmincon_test2};
            options.benchmark_id = 'unit-test';
            options.solver_names = {'sqp', 'interior-point'};
            options.max_tol_order = 2;
            options.n_runs = 2;
            options.p_type = 'ubln';
            options.maxdim = 11;
            options.mindim = 11;

            % Test the case where problem set only contains one problem.
            options.problem = s_load('ROSENBR');
            benchmark(solvers, options);
            options = rmfield(options, 'problem');

            % Test plain feature.
            benchmark(solvers);

            % Test perturbed_x0 feature.
            benchmark(solvers, 'perturbed_x0');

            % Test noisy feature.
            options.feature_name = 'noisy';
            options.noise_level = 1e-4;
            options.distribution = 'uniform';
            benchmark(solvers, options);
            options = rmfield(options, 'noise_level');
            options = rmfield(options, 'distribution');

            % Test truncated feature.
            options.summarize_log_ratio_profiles = true;
            options.feature_name = 'truncated';
            options.significant_digits = 4;
            options.perturbed_trailing_zeros = true;
            benchmark(solvers, options);
            options = rmfield(options, 'significant_digits');
            options = rmfield(options, 'perturbed_trailing_zeros');

            % Test permuted feature.
            options.summarize_performance_profiles = false;
            options.feature_name = 'permuted';
            benchmark(solvers, options);

            % Test linearly_transformed feature.
            options.summarize_data_profiles = false;
            options.summarize_performance_profiles = true;
            options.feature_name = 'linearly_transformed';
            options.rotated = true;
            options.condition_factor = 1;
            benchmark(solvers, options);
            options = rmfield(options, 'rotated');
            options = rmfield(options, 'condition_factor');

            % Test random_nan feature.
            options.summarize_performance_profiles = false;
            options.feature_name = 'random_nan';
            options.nan_rate = 0.1;
            benchmark(solvers, options);
            options = rmfield(options, 'nan_rate');

            % Test unrelaxable_constraints feature.
            options.summarize_log_ratio_profiles = false;
            options.feature_name = 'unrelaxable_constraints';
            options.unrelaxable_bounds = true;
            options.unrelaxable_linear_constraints = true;
            options.unrelaxable_nonlinear_constraints = true;
            benchmark(solvers, options);
            options = rmfield(options, 'unrelaxable_bounds');
            options = rmfield(options, 'unrelaxable_linear_constraints');
            options = rmfield(options, 'unrelaxable_nonlinear_constraints');

            % Test nonquantifiable_constraints feature.
            options.summarize_data_profiles = true;
            options.summarize_performance_profiles = false;
            options.feature_name = 'nonquantifiable_constraints';
            benchmark(solvers, options);

            % Test quantized feature.
            options.summarize_log_ratio_profiles = true;
            options.summarize_data_profiles = false;
            options.feature_name = 'quantized';
            options.mesh_size = 0.01;
            options.ground_truth = false;
            benchmark(solvers, options);
            options = rmfield(options, 'mesh_size');
            options = rmfield(options, 'ground_truth');

            % Test custom feature.
            options.feature_name = 'custom';
            options.mod_x0 = @mod_x0;
            options.mod_fun = @mod_fun;
            options.mod_affine = @mod_affine;
            benchmark(solvers, options);

            % Test `load` option.
            options.load = 'latest';
            benchmark(solvers, options);
        end

        function testErrors(testCase)
            % Test whether the function throws errors as expected.
            solvers = {@fmincon_test1, @fmincon_test2};
            options = struct();

            testCase.verifyError(@() benchmark(), "MATLAB:benchmark:solverMustBeProvided")

            testCase.verifyError(@() benchmark(solvers, {1}), "MATLAB:benchmark:SecondArgumentWrongType")

            testCase.verifyError(@() benchmark(solvers, 'plain', options, 1), "MATLAB:benchmark:TooMuchInput")

            testCase.verifyError(@() benchmark({1, 2}, options), "MATLAB:benchmark:solversWrongType")

            testCase.verifyError(@() benchmark({@fminsearch_test}), "MATLAB:benchmark:solversAtLeastTwo")

            options.feature_name = 1;
            testCase.verifyError(@() benchmark(solvers, options), "MATLAB:benchmark:feature_nameNotcharstr")
            options = rmfield(options, 'feature_name');

            options.feature_name = 'a';
            testCase.verifyError(@() benchmark(solvers, options), "MATLAB:benchmark:feature_nameNotValid")
            options = rmfield(options, 'feature_name');

            options.a = {'a'};
            testCase.verifyError(@() benchmark(solvers, options), "MATLAB:benchmark:UnknownOptions")

        end
    end
end

function x = fmincon_test1(varargin)
    options = optimoptions('fmincon', 'SpecifyObjectiveGradient', false, 'Algorithm', 'sqp');
    if nargin == 2
        fun = varargin{1};
        x0 = varargin{2};
        x = fmincon(fun, x0, [], [], [], [], [], [], [], options);
    elseif nargin == 4
        fun = varargin{1};
        x0 = varargin{2};
        xl = varargin{3};
        xu = varargin{4};
        x = fmincon(fun, x0, [], [], [], [], xl, xu, [], options);
    elseif nargin == 8
        fun = varargin{1};
        x0 = varargin{2};
        xl = varargin{3};
        xu = varargin{4};
        aub = varargin{5};
        bub = varargin{6};
        aeq = varargin{7};
        beq = varargin{8};
        x = fmincon(fun, x0, aub, bub, aeq, beq, xl, xu, [], options);
    elseif nargin == 10
        fun = varargin{1};
        x0 = varargin{2};
        xl = varargin{3};
        xu = varargin{4};
        aub = varargin{5};
        bub = varargin{6};
        aeq = varargin{7};
        beq = varargin{8};
        cub = varargin{9};
        ceq = varargin{10};
        nonlcon = @(x) deal(cub(x), ceq(x));
        x = fmincon(fun, x0, aub, bub, aeq, beq, xl, xu, nonlcon, options);
    end
end

function x = fmincon_test2(varargin)

    options = optimoptions('fmincon', 'SpecifyObjectiveGradient', false, 'Algorithm', 'interior-point');
    if nargin == 2
        fun = varargin{1};
        x0 = varargin{2};
        x = fmincon(fun, x0, [], [], [], [], [], [], [], options);
    elseif nargin == 4
        fun = varargin{1};
        x0 = varargin{2};
        xl = varargin{3};
        xu = varargin{4};
        x = fmincon(fun, x0, [], [], [], [], xl, xu, [], options);
    elseif nargin == 8
        fun = varargin{1};
        x0 = varargin{2};
        xl = varargin{3};
        xu = varargin{4};
        aub = varargin{5};
        bub = varargin{6};
        aeq = varargin{7};
        beq = varargin{8};
        x = fmincon(fun, x0, aub, bub, aeq, beq, xl, xu, [], options);
    elseif nargin == 10
        fun = varargin{1};
        x0 = varargin{2};
        xl = varargin{3};
        xu = varargin{4};
        aub = varargin{5};
        bub = varargin{6};
        aeq = varargin{7};
        beq = varargin{8};
        cub = varargin{9};
        ceq = varargin{10};
        nonlcon = @(x) deal(cub(x), ceq(x));
        x = fmincon(fun, x0, aub, bub, aeq, beq, xl, xu, nonlcon, options);
    end
end

function x0 = mod_x0(rand_stream, problem)

    [Q, R] = qr(rand_stream.randn(problem.n));
    Q(:, diag(R) < 0) = -Q(:, diag(R) < 0);
    x0 = Q * problem.x0;
    x0 = x0 + 1e-3 * max(1, norm(x0)) * rand_stream.randn(problem.n, 1) / norm(rand_stream.randn(problem.n, 1));
end

function f = mod_fun(x, rand_stream, problem)

    f = problem.fun(x);
    f = f + max(1, abs(f)) * 1e-3 * rand_stream.randn(1);
end

function [A, b, inv] = mod_affine(rand_stream, problem)

    [Q, R] = qr(rand_stream.randn(problem.n));
    Q(:, diag(R) < 0) = -Q(:, diag(R) < 0);
    A = Q';
    b = zeros(problem.n, 1);
    inv = Q;
end