classdef Feature < handle
%FEATURE is a class that defines a mapping from an optimization problem to a
%   new one with specified features.
%
%   We are interested to test solvers on problems with different features. For
%   example, we want to test the performance of solvers under the case where
%   the objective function is noisy. For this purpose, we define FEATURE class.
%
%   Suppose we have an optimization problem
%
%       min fun(x)
%       s.t. xl <= x <= xu,
%            aub * x <= bub,
%            aeq * x = beq,
%            cub(x) <= 0,
%            ceq(x) = 0,
%       with initial point x0.
%
%   Then, FEATURE maps the above optimization problem to the following one.
%
%       min fun_mod(A * x + b)
%       s.t. xl_mod <= A * x + b <= xu_mod,
%            aub_mod * (A * x + b) <= bub_mod,
%            aeq_mod * (A * x + b) = beq_mod,
%            cub_mod(A * x + b) <= 0,
%            ceq_mod(A * x + b) = 0,
%       with initial point x0_mod.
%
%   FEATURE should be initialized by the following signature:
%
%       F = FEATURE(NAME) creates an instance of FEATURE with the name NAME;
%
%       F = FEATURE(NAME, OPTIONS) creates an instance of FEATURE with the
%       name NAME and options OPTIONS.
%
%   The input NAME should be one of the following char or string:
%
%       1. 'plain':
%           do nothing to the optimization problem.
%       2. 'perturbed_x0':
%           perturb the initial point.
%       3. 'noisy':
%           add noise to the objective function and nonlinear constraints.
%       4. 'truncated':
%           truncate values of the objective function and nonlinear constraints
%           to a given number of significant digits.
%       5. 'permuted':
%           randomly permute the variables. Note that x0, xl, xu, aub, bub,
%           aeq, beq will be modified since we want to keep the new problem
%           mathematically equivalent to the original one.
%       6. 'linearly_transformed':
%           generate a invertible linear transformation in the form D * Q' with
%           D being a diagonal matrix and Q being a random orthogonal matrix,
%           and apply the transformation to the variables. In this way, the
%           Hessian of the objective function becomes Q * D * H * D * Q' with H
%           being the original Hessian. Note that x0, xl, xu, aub, bub, aeq,
%           beq will be modified since we want to keep the new problem
%           mathematically equivalent to the original one.
%       7. 'random_nan':
%           randomly replace values of the objective function and nonlinear
%           constraints with NaN.
%       8. 'unrelaxable_constraints':
%           set the objective function to Inf outside the feasible region.
%       9. 'nonquantifiable_constraints':
%           replace values of nonlinear constraints with either 0 (if the
%           constraint is satisfied) or 1 (if the constraint is violated).
%       10. 'quantized':
%           quantize the objective function and nonlinear constraints.
%       11. 'custom':
%           user-defined FEATURE.
%
%   The optional input OPTIONS should be a struct which can contain the
%   following fields:
%
%       - n_runs: the number of runs of the experiments under the given
%         feature. Default is 5 for stochastic features and 1 for deterministic
%         features.
%       - distribution: the distribution of perturbation in 'perturbed_x0'
%         feature or noise in 'noisy' feature. It should be either a string
%         (or char), or a function handle
%               ``(random_stream, dimension) -> random vector``,
%         accepting a random_stream and the dimension of a problem and
%         returning a random vector with the given dimension. In 'perturbed_x0'
%         case, the char should be either 'spherical' or 'gaussian' (default is
%         'spherical'). In 'noisy' case, the char should be either 'gaussian'
%         or 'uniform' (default is 'gaussian').
%       - perturbation_level: the magnitude of the perturbation to the initial
%         guess in the 'perturbed_x0' feature. Default is 10^-3.
%       - noise_level: the magnitude of the noise in the 'noisy' feature.
%         Default is 10^-3.
%       - noise_type: the type of the noise in the 'noisy' features. It should
%         be either 'absolute', 'relative', or 'mixed'. Default is 'mixed'.
%       - significant_digits: the number of significant digits in the
%         'truncated' feature. Default is 6.
%       - perturbed_trailing_digits: whether we will randomize the trailing
%         zeros of the objective function value in the 'perturbed_x0' feature.
%         Default is false.
%       - rotated: whether to use a random or given rotation matrix to rotate
%         the coordinates of a problem in the 'linearly_transformed' feature.
%         Default is true.
%       - condition_factor: the scaling factor of the condition number of the
%         linear transformation in the 'linearly_transformed' feature. More
%         specifically, the condition number of the linear transformation will
%         2 ^ (condition_factor * n / 2), where `n` is the dimension of the
%         problem. Default is 0.
%       - nan_rate: the probability that the evaluation of the objective
%         function will return NaN in the 'random_nan' feature. Default is
%         0.05.
%       - unrelaxable_bounds: whether the bound constraints are unrelaxable or
%         not in the 'unrelaxable_constraints' feature. Default is true.
%       - unrelaxable_linear_constraints: whether the linear constraints are
%         unrelaxable or not in the 'unrelaxable_constraints' feature. Default
%         is false.
%       - unrelaxable_nonlinear_constraints: whether the nonlinear constraints
%         are unrelaxable or not in the 'unrelaxable_constraints' feature.
%         Default is false.
%       - mesh_size: the size of the mesh in the 'quantized' feature. Default
%         is 10^-3.
%       - mesh_type: the type of the mesh in the 'quantized' feature. It should
%         be either 'absolute' or 'relative'. Default is 'absolute'.
%       - ground_truth: whether the feature is the ground truth or not. Default
%         is true.
%       - mod_x0: the modifier function to modify the inital guess in the 
%         'custom' feature. It should be a function handle as follows:
%               ``(random_stream, problem) -> modified_x0``,
%         where `problem` is an instance of the class Problem, and
%         `modified_x0` is the modified initial guess. No default.
%       - mod_affine: the modifier function to generate the affine
%         transformation applied to the variables in the 'custom' feature. It
%         should be a function handle as follows:
%               ``(random_stream, problem) -> (A, b, inv)``,
%         where `problem` is an instance of the class Problem, `A` is the
%         matrix of the affine transformation, `b` is the vector of the affine
%         transformation, and `inv` is the inverse of matrix `A`. No default.
%       - mod_bounds: the modifier function to modify the bound constraints in
%         the 'custom' feature. It should be a function handle as follows:
%               ``(random_stream, problem) -> (modified_xl, modified_xu)``,
%         where `problem` is an instance of the class Problem, `modified_xl` is
%         the modified lower bound, and `modified_xu` is the modified upper
%         bound. No default.
%       - mod_linear_ub: the modifier function to modify the linear inequality
%         constraints in the 'custom' feature. It should be a function handle
%         as follows:
%               ``(random_stream, problem) -> (modified_aub, modified_bub)``,
%         where `problem` is an instance of the class Problem, `modified_aub`
%         is the modified matrix of the linear inequality constraints, and
%         `modified_bub` is the modified vector of the linear inequality
%         constraints. No default.
%       - mod_linear_eq: the modifier function to modify the linear equality
%         constraints in the 'custom' feature. It should be a function handle
%         as follows:
%               ``(random_stream, problem) -> (modified_aeq, modified_beq)``,
%         where `problem` is an instance of the class Problem, `modified_aeq`
%         is the modified matrix of the linear equality constraints, and
%         `modified_beq` is the modified vector of the linear equality
%         constraints. No default.
%       - mod_fun: the modifier function to modify the objective function in
%         the 'custom' feature. It should be a function handle as follows:
%               ``(x, random_stream, problem) -> modified_fun``,
%         where `x` is the evaluation point, `problem` is an instance of the
%         class Problem, and `modified_fun` is the modified objective function
%         value. No default.
%       - mod_cub: the modifier function to modify the nonlinear inequality
%         constraints in the 'custom' feature. It should be a function handle
%         as follows:
%               ``(x, random_stream, problem) -> modified_cub``,
%         where x is the evaluation point, `problem` is an instance of the
%         class Problem, and `modified_cub` is the modified vector of the
%         nonlinear inequality constraints. No default.
%       - mod_ceq: the modifier function to modify the nonlinear equality
%         constraints in the 'custom' feature. It should be a function handle
%         as follows:
%               ``(x, random_stream, problem) -> modified_ceq``,
%         where x is the evaluation point, `problem` is an instance of the
%         class Problem, and `modified_ceq` is the modified vector of the
%         nonlinear equality constraints. No default.
%
%   Different input NAME will have different valid fields of OPTIONS. We list
%   the valid fields for each input NAME as follows:
%
%       1. 'plain':
%           n_runs
%       2. 'perturbed_x0':
%           n_runs, distribution, perturbation_level
%       3. 'noisy':
%           n_runs, distribution, noise_level, noise_type
%       4. 'truncated':
%           n_runs, significant_digits, perturbed_trailing_digits
%       5. 'permuted':
%           n_runs
%       6. 'linearly_transformed':
%           n_runs, rotated, condition_factor
%       7. 'random_nan':
%           n_runs, nan_rate
%       8. 'unrelaxable_constraints':
%           n_runs, unrelaxable_bounds, unrelaxable_linear_constraints,
%           unrelaxable_nonlinear_constraints
%       9. 'nonquantifiable_constraints':
%           n_runs
%       10. 'quantized':
%           n_runs, mesh_size, mesh_type, ground_truth
%       11. 'custom':
%           n_runs, mod_x0, mod_affine, mod_bounds, mod_linear_ub,
%           mod_linear_eq, mod_fun, mod_cub, mod_ceq
%
%   The output F contains the following methods:
%
%       1. modifier_x0:
%           a function handle to modify the initial point.
%       2. modifier_affine:
%           a function handle to generate an invertible matrix A (and its
%           inverse) and a vector b for the affine transformation applied to
%           the variables.
%       3. modifier_bounds:
%           a function handle to modify the lower and upper bounds.
%       4. modifier_linear_ub:
%           a function handle to modify the linear inequality constraints.
%       5. modifier_linear_eq:
%           a function handle to modify the linear equality constraints.
%       6. modifier_fun:
%           a function handle to modify the objective function value.
%       7. modifier_cub:
%           a function handle to modify the values of the nonlinear inequality
%           constraints.
%       8. modifier_ceq:
%           a function handle to modify the values of the nonlinear equality
%           constraints.
%
%   All the methods of F will be used later to modify the optimization
%   problem.
%

    properties (GetAccess = public, SetAccess = public)

        name
        options
        
    end

    methods

        function obj = Feature(name, varargin)
            %{
            Initialize a FEATURE.

            Parameters
            ----------
            name : char or str
                Name of the FEATURE. The legal values are 'plain',
                'perturbed_x0', 'noisy', 'truncated', 'permuted', 
                'linearly_transformed', 'random_nan',
                'unrelaxable_constraints', 'nonquantifiable_constraints',
                'custom'.
            varargin : cell
                Options of the FEATURE.

            Returns
            -------
            obj : Feature
                A FEATURE object.
            %}

            % Preprocess the feature name.
            if ~ischarstr(name)
                error("MATLAB:Feature:FeaturenameNotString", "The first input argument for `Feature` must be a char or string.")
            end
            % Convert the feature name to lowercase characters.
            obj.name = char(lower(name));
            % Check whether the feature is valid.
            % validFeatureNames = {enumeration('FeatureName').value};
            % Only for MATLAB R2021b and later.
            validFeatureNames = cellfun(@(x) x.value, num2cell(enumeration('FeatureName')), 'UniformOutput', false);
            if ~ismember(obj.name, validFeatureNames)
                error("MATLAB:Feature:UnknownFeature", "Unknown feature: " + obj.name + ".")
            end

            % Preprocess the feature options.
            obj.options = struct();
            if nargin > 1
                if isstruct(varargin{1})
                    obj.options = varargin{1};
                elseif mod(length(varargin), 2) ~= 0
                    error("MATLAB:Feature:InvalidNumberOfArguments", "The number of input arguments for `Feature` must be odd if the second argument is not a struct.")
                else
                    for i = 1:2:length(varargin)
                        obj.options.(lower(varargin{i})) = varargin{i+1};
                    end
                end
            end
            optionKeys = fieldnames(obj.options);
            validOptionKeys = cellfun(@(x) x.value, num2cell(enumeration('FeatureOptionKey')), 'UniformOutput', false);
            for i = 1:numel(optionKeys)
                key = optionKeys{i};
                if ~ismember(key, validOptionKeys)
                    error("MATLAB:Feature:UnknownOption", "Unknown option for feature: " + key + ".")
                end

                % Check whether the option is valid for the feature.
                known_options = {FeatureOptionKey.N_RUNS.value};
                switch obj.name
                    case FeatureName.CUSTOM.value
                        known_options = [known_options, {FeatureOptionKey.MOD_X0.value, FeatureOptionKey.MOD_BOUNDS.value, FeatureOptionKey.MOD_LINEAR_UB.value, FeatureOptionKey.MOD_LINEAR_EQ.value, FeatureOptionKey.MOD_AFFINE.value, FeatureOptionKey.MOD_FUN.value, FeatureOptionKey.MOD_CUB.value, FeatureOptionKey.MOD_CEQ.value}];
                    case FeatureName.NOISY.value
                        known_options = [known_options, {FeatureOptionKey.DISTRIBUTION.value, FeatureOptionKey.NOISE_LEVEL.value, FeatureOptionKey.NOISE_TYPE.value}];
                    case FeatureName.PERTURBED_X0.value
                        known_options = [known_options, {FeatureOptionKey.DISTRIBUTION.value, FeatureOptionKey.PERTURBATION_LEVEL.value}];
                    case FeatureName.RANDOM_NAN.value
                        known_options = [known_options, {FeatureOptionKey.NAN_RATE.value}];
                    case FeatureName.TRUNCATED.value
                        known_options = [known_options, {FeatureOptionKey.PERTURBED_TRAILING_DIGITS.value, FeatureOptionKey.SIGNIFICANT_DIGITS.value}];
                    case FeatureName.UNRELAXABLE_CONSTRAINTS.value
                        known_options = [known_options, {FeatureOptionKey.UNRELAXABLE_BOUNDS.value, FeatureOptionKey.UNRELAXABLE_LINEAR_CONSTRAINTS.value, FeatureOptionKey.UNRELAXABLE_NONLINEAR_CONSTRAINTS.value}];
                    case FeatureName.LINEARLY_TRANSFORMED.value
                        known_options = [known_options, {FeatureOptionKey.ROTATED.value, FeatureOptionKey.CONDITION_FACTOR.value}];
                    case FeatureName.QUANTIZED.value
                        known_options = [known_options, {FeatureOptionKey.MESH_SIZE.value, FeatureOptionKey.MESH_TYPE.value, FeatureOptionKey.GROUND_TRUTH.value}];
                    case FeatureName.PERMUTED.value
                        % Do nothing
                    case FeatureName.NONQUANTIFIABLE_CONSTRAINTS.value
                        % Do nothing
                    case FeatureName.PLAIN.value
                        % Do nothing
                    otherwise
                        error("MATLAB:Feature:UnknownFeature", "Unknown feature: " + obj.name + ".")
                end
                if ~ismember(key, known_options)
                    error("MATLAB:Feature:InvalidOptionForFeature", "Option `" + key + "` is not valid for feature '" + obj.name + "'.")
                end

                % Check whether the option type is valid.
                switch key
                    case FeatureOptionKey.N_RUNS.value
                        if ~isintegerscalar(obj.options.(key)) || obj.options.(key) <= 0
                            error("MATLAB:Feature:n_runs_NotPositiveInteger", "Option `" + key + "` must be a positive integer.")
                        end
                    case FeatureOptionKey.DISTRIBUTION.value
                        if ~isa(obj.options.(key), 'function_handle') && ~ischarstr(obj.options.(key))
                            error("MATLAB:Feature:distribution_NotFunctionHandle", "Option `" + key + "` must be the char (or string), or a function handle.")
                        end
                        if strcmp(obj.name, FeatureName.NOISY.value) && ~ismember(char(obj.options.(key)), {'gaussian', 'uniform'})
                            error("MATLAB:Feature:distribution_NotFunctionHandle_noisy", "Option `" + key + "` for feature 'noisy' must be the char (or string) 'gaussian' or 'uniform', or a function handle.")
                        elseif strcmp(obj.name, FeatureName.PERTURBED_X0.value) && ~ismember(char(obj.options.(key)), {'gaussian', 'spherical'})
                            error("MATLAB:Feature:distribution_NotFunctionHandle_perturbed_x0", "Option `" + key + "` for feature 'perturbed_x0' must be the char (or string) 'gaussian' or 'spherical', or a function handle.")
                        end
                        if ischarstr(obj.options.(key))
                            obj.options.(key) = char(obj.options.(key));
                        end
                    case FeatureOptionKey.NAN_RATE.value
                        if ~isrealscalar(obj.options.(key)) || obj.options.(key) < 0.0 || obj.options.(key) > 1.0
                            error("MATLAB:Feature:nan_rate_NotBetween_0_1", "Option `" + key + "` must be a real number between 0 and 1.")
                        end
                    case FeatureOptionKey.SIGNIFICANT_DIGITS.value
                        if ~isintegerscalar(obj.options.(key)) || obj.options.(key) <= 0
                            error("MATLAB:Feature:significant_digits_NotPositiveInteger", "Option `" + key + "` must be a positive integer.")
                        end
                    case FeatureOptionKey.NOISE_LEVEL.value
                        if ~isrealscalar(obj.options.(key)) || obj.options.(key) < 0.0
                            error("MATLAB:Feature:noise_level_NotPositive", "Option `" + key + "` must be a nonnegative real number.")
                        end
                    case FeatureOptionKey.NOISE_TYPE.value
                        if ~ischarstr(obj.options.(key)) || ~ismember(obj.options.(key), {'absolute', 'relative', 'mixed'})
                            error("MATLAB:Feature:noise_type_InvalidInput", "Option `" + key + "` must be 'absolute', 'relative', or 'mixed'.")
                        end
                        obj.options.(key) = char(obj.options.(key));
                    case FeatureOptionKey.PERTURBED_TRAILING_DIGITS.value
                        if ~islogicalscalar(obj.options.(key))
                            error("MATLAB:Feature:perturbed_trailing_digits_NotLogical", "Option `" + key + "` must be a logical value.")
                        end
                    case FeatureOptionKey.ROTATED.value
                        if ~islogicalscalar(obj.options.(key))
                            error("MATLAB:Feature:rotated_NotLogical", "Option `" + key + "` must be a logical value.")
                        end
                    case FeatureOptionKey.CONDITION_FACTOR.value
                        if ~(isrealscalar(obj.options.(key)) && obj.options.(key) >= 0)
                            error("MATLAB:Feature:condition_factor_InvalidInput", "Option `" + key + "` must be a nonnegative real number.")
                        end
                    case FeatureOptionKey.UNRELAXABLE_BOUNDS.value
                        if ~islogicalscalar(obj.options.(key))
                            error("MATLAB:Feature:unrelaxable_bounds_NotLogical", "Option `" + key + "` must be a logical value.")
                        end
                    case FeatureOptionKey.UNRELAXABLE_LINEAR_CONSTRAINTS.value
                        if ~islogicalscalar(obj.options.(key))
                            error("MATLAB:Feature:unrelaxable_linear_constraints_NotLogical", "Option `" + key + "` must be a logical value.")
                        end
                    case FeatureOptionKey.UNRELAXABLE_NONLINEAR_CONSTRAINTS.value
                        if ~islogicalscalar(obj.options.(key))
                            error("MATLAB:Feature:unrelaxable_nonlinear_constraints_NotLogical", "Option `" + key + "` must be a logical value.")
                        end
                    case FeatureOptionKey.MESH_SIZE.value
                        if ~isrealscalar(obj.options.(key)) || obj.options.(key) <= 0.0
                            error("MATLAB:Feature:mesh_size_NotPositive", "Option `" + key + "` must be a positive real number.")
                        end
                    case FeatureOptionKey.MESH_TYPE.value
                        if ~ischarstr(obj.options.(key)) || ~ismember(obj.options.(key), {'absolute', 'relative'})
                            error("MATLAB:Feature:mesh_type_InvalidInput", "Option `" + key + "` must be 'absolute' or 'relative'.")
                        end
                        obj.options.(key) = char(obj.options.(key));
                    case FeatureOptionKey.GROUND_TRUTH.value
                        if ~islogicalscalar(obj.options.(key))
                            error("MATLAB:Feature:ground_truth_NotLogical", "Option `" + key + "` must be a logical value.")
                        end
                    case FeatureOptionKey.MOD_X0.value
                        if ~isa(obj.options.(key), 'function_handle')
                            error("MATLAB:Feature:mod_x0_NotFunctionHandle", "Option `" + key + "` must be a function handle.")
                        end
                    case FeatureOptionKey.MOD_BOUNDS.value
                        if ~isa(obj.options.(key), 'function_handle')
                            error("MATLAB:Feature:mod_bounds_NotFunctionHandle", "Option `" + key + "` must be a function handle.")
                        end
                    case FeatureOptionKey.MOD_LINEAR_UB.value
                        if ~isa(obj.options.(key), 'function_handle')
                            error("MATLAB:Feature:mod_linear_ub_NotFunctionHandle", "Option `" + key + "` must be a function handle.")
                        end
                    case FeatureOptionKey.MOD_LINEAR_EQ.value
                        if ~isa(obj.options.(key), 'function_handle')
                            error("MATLAB:Feature:mod_linear_eq_NotFunctionHandle", "Option `" + key + "` must be a function handle.")
                        end
                    case FeatureOptionKey.MOD_AFFINE.value
                        if ~isa(obj.options.(key), 'function_handle')
                            error("MATLAB:Feature:mod_affine_NotFunctionHandle", "Option `" + key + "` must be a function handle.")
                        end
                    case FeatureOptionKey.MOD_FUN.value
                        if ~isa(obj.options.(key), 'function_handle')
                            error("MATLAB:Feature:mod_fun_NotFunctionHandle", "Option `" + key + "` must be a function handle.")
                        end
                    case FeatureOptionKey.MOD_CUB.value
                        if ~isa(obj.options.(key), 'function_handle')
                            error("MATLAB:Feature:mod_cub_NotFunctionHandle", "Option `" + key + "` must be a function handle.")
                        end
                    case FeatureOptionKey.MOD_CEQ.value
                        if ~isa(obj.options.(key), 'function_handle')
                            error("MATLAB:Feature:mod_ceq_NotFunctionHandle", "Option `" + key + "` must be a function handle.")
                        end
                end
            end

            % Set default options for the unspecified options.
            obj = set_default_options(obj);
        end

        function is_stochastic = is_stochastic(obj)
            %{
            Determine whether the feature is stochastic.

            Returns
            -------
            is_stochastic : bool
            %}
            switch obj.name
                case {FeatureName.PERTURBED_X0.value, FeatureName.NOISY.value, FeatureName.PERMUTED.value, FeatureName.RANDOM_NAN.value, FeatureName.CUSTOM.value}
                    is_stochastic = true;
                case FeatureName.TRUNCATED.value
                    is_stochastic = obj.options.(FeatureOptionKey.PERTURBED_TRAILING_DIGITS.value);
                case FeatureName.LINEARLY_TRANSFORMED.value
                    is_stochastic = obj.options.(FeatureOptionKey.ROTATED.value);
                otherwise
                    is_stochastic = false;
            end
        end

        function x0 = modifier_x0(obj, seed, problem)
            %{
            Modify the initial point.

            Parameters
            ----------
            seed : int
                Seed used to generate random numbers.
            problem : Problem
                Problem for which the initial point is modified.
            
            Returns
            -------
            x0 : double, size (n,)
                Modified initial point.
            %}

            switch obj.name
                case FeatureName.CUSTOM.value
                    % If the user specifies a custom modifier for the initial point, use it.
                    if isfield(obj.options, FeatureOptionKey.MOD_X0.value)
                        rand_stream_custom = obj.default_rng(seed);
                        x0 = obj.options.(FeatureOptionKey.MOD_X0.value)(rand_stream_custom, problem);
                        return;
                    end
                    % If the user does not specify a custom modifier for the initial point but specifies
                    % a custom affine transformation, we need to apply the inverse of the affine
                    % transformation to the initial point.
                    if isfield(obj.options, FeatureOptionKey.MOD_AFFINE.value)
                        [~, b, inv] = obj.modifier_affine(seed, problem);
                        x0 = inv * (problem.x0 - b);
                    else
                        x0 = problem.x0;
                    end
                case FeatureName.PERTURBED_X0.value
                    % Use max(1, norm(x0)) to avoid no perturbation when x0 is zero.
                    rand_stream_perturbed_x0 = obj.default_rng(seed);
                    perturbation_level = obj.options.(FeatureOptionKey.PERTURBATION_LEVEL.value) * max(1, norm(problem.x0));
                    if strcmp(obj.options.(FeatureOptionKey.DISTRIBUTION.value), 'gaussian')
                        x0 = problem.x0 + perturbation_level * rand_stream_perturbed_x0.randn(problem.n, 1);
                    elseif strcmp(obj.options.(FeatureOptionKey.DISTRIBUTION.value), 'spherical')
                        perturbation = rand_stream_perturbed_x0.randn(problem.n, 1);
                        x0 = problem.x0 + perturbation_level * perturbation / norm(perturbation);
                    else
                        x0 = problem.x0 + perturbation_level * obj.options.(FeatureOptionKey.DISTRIBUTION.value)(rand_stream_perturbed_x0, problem.n);
                    end
                case FeatureName.PERMUTED.value
                    % Note that we need to apply the reverse permutation to the initial point so that
                    % the new problem is mathematically equivalent to the original one.
                    rand_stream_permuted = obj.default_rng(seed);
                    permutation = rand_stream_permuted.randperm(problem.n);
                    [~, reverse_permutation] = sort(permutation);
                    x0 = problem.x0(reverse_permutation);
                case FeatureName.LINEARLY_TRANSFORMED.value
                    % Apply the inverse of the affine transformation to the initial point.
                    [~, ~, inv] = obj.modifier_affine(seed, problem);
                    x0 = inv * problem.x0;
                otherwise
                    x0 = problem.x0;
            end
        end

        function [A, b, inv] = modifier_affine(obj, seed, problem)
            %{
            Generate an invertible matrix A and a vector b for the affine transformation applied to the variables.

            Parameters
            ----------
            seed : int
                Seed used to generate random numbers.
            problem : Problem
                Problem for which the affine transformation is generated.
            Returns
            -------
            A : double, size (n, n)
                Matrix of the affine transformation.
            b : double, size (n, 1)
                Vector of the affine transformation.
            inv : double, size (n, n)
                Inverse of the matrix A.
            %}

            % Default values
            A = eye(problem.n);
            b = zeros(problem.n, 1);
            inv = eye(problem.n);

            switch obj.name
                case FeatureName.CUSTOM.value
                    if isfield(obj.options, FeatureOptionKey.MOD_AFFINE.value)
                        rand_stream_custom = obj.default_rng(seed);
                        [A, b, inv] = obj.options.(FeatureOptionKey.MOD_AFFINE.value)(rand_stream_custom, problem);
                    end
                    % Check whether A * inv is an identity matrix.
                    if norm(A * inv - eye(problem.n)) > 1e-8 * problem.n
                        error("MATLAB:Feature:AffineTransformationNotInvertible", "The multiplication of the affine transformation matrix and its inverse is not an identity matrix.")
                    end
                case FeatureName.PERMUTED.value
                    % Generate a random permutation matrix.
                    rand_stream_permuted = obj.default_rng(seed);
                    permutation = rand_stream_permuted.randperm(problem.n);
                    A = eye(problem.n);
                    A = A(permutation, :);
                    inv = A';
                    b = zeros(problem.n, 1);
                case FeatureName.LINEARLY_TRANSFORMED.value
                    % Generate A in the form D * Q' with D being a diagonal matrix and Q being an
                    % orthogonal matrix.

                    % Generate a random rotation matrix Q if the 'rotated' option is set to true.
                    if obj.options.(FeatureOptionKey.ROTATED.value)
                        %{
                        We generate a random orthogonal matrix Q following the uniform distribution
                        on O(n). The method refers to a note written by the late Professor Nicholas
                        Higham, see:
                        https://nhigham.com/2020/04/22/what-is-a-random-orthogonal-matrix/
                        and an answer from MathStackExchange, see:
                        https://math.stackexchange.com/a/4891933/1088047
                        %}
                        rand_stream_linearly_transformed = obj.default_rng(seed);
                        [Q, R] = qr(rand_stream_linearly_transformed.randn(problem.n));
                        Q(:, diag(R) < 0) = -Q(:, diag(R) < 0);
                    else
                        Q = eye(problem.n);
                    end
                    % Generate a positive definite diagonal matrix D with condition number equal to
                    % 2^(condition_factor * n / 2), where n is the dimension of the problem. In this
                    % way, the condition number of Q * D^2 * Q' is 2^(condition_factor * n).
                    log_condition_number = sqrt(obj.options.(FeatureOptionKey.CONDITION_FACTOR.value) * problem.n / 2);
                    power = linspace(-log_condition_number/2, log_condition_number/2, problem.n);
                    A = diag(2.^power) * Q';
                    inv = Q * diag(2.^-power);
                otherwise
                    % Do nothing
            end
        end

        function [xl, xu] = modifier_bounds(obj, seed, problem)
            %{
            Modify the bounds.

            Parameters
            ----------
            seed : int
                Seed used to generate random numbers.
            problem : Problem
                Problem for which the bounds are modified.

            Returns
            -------
            xl : double, size (n,)
                Modified lower bounds.
            xu : double, size (n,)
                Modified upper bounds.
            %}

            switch obj.name
                case FeatureName.CUSTOM.value
                    % If the user specifies a custom modifier for the bounds, use it.
                    if isfield(obj.options, FeatureOptionKey.MOD_BOUNDS.value)
                        rand_stream_custom = obj.default_rng(seed);
                        [xl, xu] = obj.options.(FeatureOptionKey.MOD_BOUNDS.value)(rand_stream_custom, problem);
                        return;
                    end
                    if ~isfield(obj.options, FeatureOptionKey.MOD_AFFINE.value)
                        xl = problem.xl;
                        xu = problem.xu;
                        return;
                    end
                    % If the user does not specify a custom modifier for the bounds but specifies a
                    % custom affine transformation, we need to specially handle the bounds.
                    [~, b, inv] = obj.modifier_affine(seed, problem);
                    if ~isdiag(inv)
                        xl = -Inf(problem.n, 1);
                        xu = Inf(problem.n, 1);
                        return;
                    end
                    % If the inverse of the affine transformation is diagonal, we can apply it to get
                    % the modified bounds.
                    xl_tmp = diag(inv) .* (problem.xl - b);
                    xu_tmp = diag(inv) .* (problem.xu - b);
                    xl = min(xl_tmp, xu_tmp);
                    xu = max(xl_tmp, xu_tmp);
                case FeatureName.PERMUTED.value
                    % Note that we need to apply the reverse permutation to the bounds so that the new
                    % problem is mathematically equivalent to the original one.
                    rand_stream_permuted = obj.default_rng(seed);
                    permutation = rand_stream_permuted.randperm(problem.n);
                    [~, reverse_permutation] = sort(permutation);
                    xl = problem.xl(reverse_permutation);
                    xu = problem.xu(reverse_permutation);
                case FeatureName.LINEARLY_TRANSFORMED.value
                    % Apply the inverse of the affine transformation to the bounds.
                    [~, ~, inv] = obj.modifier_affine(seed, problem);
                    if ~isdiag(inv)
                        xl = -Inf(problem.n, 1);
                        xu = Inf(problem.n, 1);
                        return;
                    end
                    xl_tmp = diag(inv) .* problem.xl;
                    xu_tmp = diag(inv) .* problem.xu;
                    xl = min(xl_tmp, xu_tmp);
                    xu = max(xl_tmp, xu_tmp);
                otherwise
                    xl = problem.xl;
                    xu = problem.xu;    
            end
        end

        function [aub, bub] = modifier_linear_ub(obj, seed, problem)
            %{
            Modify the linear inequality constraints.

            Parameters
            ----------
            seed : int
                Seed used to generate random numbers.
            problem : Problem
                Problem for which the linear inequality constraints are modified.

            Returns
            -------
            aub : double, size (m_linear_ub, n)
                Modified matrix of linear inequality constraints.
            bub : double, size (m_linear_ub,)
                Modified right-hand side vector of linear inequality
                constraints.
            %}
            
            switch obj.name
                case FeatureName.CUSTOM.value
                    % If the user specifies a custom modifier for the linear inequality constraints,
                    % use it.
                    if isfield(obj.options, FeatureOptionKey.MOD_LINEAR_UB.value)
                        rand_stream_custom = obj.default_rng(seed);
                        [aub, bub] = obj.options.(FeatureOptionKey.MOD_LINEAR_UB.value)(rand_stream_custom, problem);
                        return;
                    end
                    if ~isfield(obj.options, FeatureOptionKey.MOD_AFFINE.value)
                        aub = problem.aub;
                        bub = problem.bub;
                        return;
                    end
                    % If the user does not specify a custom modifier for the linear inequality
                    % constraints but specifies a custom affine transformation, we need to specially
                    % handle the linear inequality constraints.
                    [A, b] = obj.modifier_affine(seed, problem);
                    if isdiag(A)
                        aub = problem.aub * A;
                        bub = problem.bub - problem.aub * b;
                        return;
                    end
                    %{
                    We need to specially handle bound constraints and linear inequality constraints.

                    Bound constraints
                    xl <= A * x + b <= xu
                    should be modified to
                    A * x <= xu - b
                    and
                    -A * x <= -xl + b
                    when A is not diagonal.

                    Linear inequality constraints
                    turn out to be
                    (aub * A) * x <= bub - aub * b
                    %}
                    
                    % Pick out the indices of lower bounds who are not -Inf and upper bounds who are
                    % not Inf since later we will not transform them into linear inequality constraints.
                    idx_lb = ~isinf(problem.xl);
                    idx_ub = ~isinf(problem.xu);
                    % Remove the indices, of which the lower and upper bound are equal, since later
                    % we will put them into the linear equality constraints.
                    idx_eq = find(problem.xl == problem.xu);
                    idx_lb(idx_eq) = false;
                    idx_ub(idx_eq) = false;
                    if isempty(problem.aub)
                        aub = [A(idx_ub, :); -A(idx_lb, :)];
                        bub = [problem.xu(idx_ub) - b(idx_ub); -(problem.xl(idx_lb) - b(idx_lb))];
                        return;
                    end
                    aub = [A(idx_ub, :); -A(idx_lb, :); problem.aub * A];
                    bub = [problem.xu(idx_ub) - b(idx_ub); -(problem.xl(idx_lb) - b(idx_lb)); problem.bub - problem.aub * b];
                case FeatureName.PERMUTED.value
                    rand_stream_permuted = obj.default_rng(seed);
                    permutation = rand_stream_permuted.randperm(problem.n);
                    [~, reverse_permutation] = sort(permutation);
                    aub = problem.aub(:, reverse_permutation);
                    bub = problem.bub;
                case FeatureName.LINEARLY_TRANSFORMED.value
                    % Similar to the case in the custom feature where a custom affine transformation
                    % is specified.
                    A = obj.modifier_affine(seed, problem);
                    if isdiag(A)
                        aub = problem.aub * A;
                        bub = problem.bub;
                        return;
                    end
                    idx_lb = ~isinf(problem.xl);
                    idx_ub = ~isinf(problem.xu);
                    idx_eq = find(problem.xl == problem.xu);
                    idx_lb(idx_eq) = false;
                    idx_ub(idx_eq) = false;
                    if isempty(problem.aub)
                        aub = [A(idx_ub, :); -A(idx_lb, :)];
                        bub = [problem.xu(idx_ub); -problem.xl(idx_lb)];
                        return;
                    end
                    aub = [A(idx_ub, :); -A(idx_lb, :); problem.aub * A];
                    bub = [problem.xu(idx_ub); -problem.xl(idx_lb); problem.bub];
                otherwise
                    aub = problem.aub;
                    bub = problem.bub;
            end
        end

        function [aeq, beq] = modifier_linear_eq(obj, seed, problem)
            %{
            Modify the linear equality constraints.

            Parameters
            ----------
            seed : int
                Seed used to generate random numbers.
            problem : Problem
                Problem for which the linear equality constraints are modified.

            Returns
            -------
            aeq : double, size (m_linear_eq, n)
                Modified matrix of linear equality constraints.
            beq : double, size (m_linear_eq,)
                Modified right-hand side vector of linear equality constraints.
            %}
            
            switch obj.name
                case FeatureName.CUSTOM.value
                    % If the user specifies a custom modifier for the linear equality constraints, use it.
                    if isfield(obj.options, FeatureOptionKey.MOD_LINEAR_EQ.value)
                        rand_stream_custom = obj.default_rng(seed);
                        [aeq, beq] = obj.options.(FeatureOptionKey.MOD_LINEAR_EQ.value)(rand_stream_custom, problem);
                        return;
                    end
                    if ~isfield(obj.options, FeatureOptionKey.MOD_AFFINE.value)
                        aeq = problem.aeq;
                        beq = problem.beq;
                        return;
                    end
                    % If the user does not specify a custom modifier for the linear equality constraints but specifies a custom affine transformation, we need to specially handle the linear equality constraints.
                    [A, b] = obj.modifier_affine(seed, problem);
                    if isdiag(A)
                        aeq = problem.aeq * A;
                        beq = problem.beq - problem.aeq * b;
                        return;
                    end
                    %{
                    We need to specially handle bound constraints and linear equality constraints.

                    Bound constraints
                    xl <= A * x + b <= xu
                    with xl = xu should be modified to
                    A * x = xu - b
                    when A is not diagonal.

                    Linear equality constraints
                    turn out to be
                    (aeq * A) * x = beq - aeq * b
                    %}

                    % Pick out the indices, of which the lower and upper bound are equal.
                    idx_eq = find(problem.xl == problem.xu);
                    if isempty(problem.aeq)
                        aeq = A(idx_eq, :);
                        beq = problem.xu(idx_eq) - b(idx_eq);
                        return;
                    end
                    aeq = [A(idx_eq, :); problem.aeq * A];
                    beq = [problem.xu(idx_eq) - b(idx_eq); problem.beq - problem.aeq * b];
                case FeatureName.PERMUTED.value
                    rand_stream_permuted = obj.default_rng(seed);
                    permutation = rand_stream_permuted.randperm(problem.n);
                    [~, reverse_permutation] = sort(permutation);
                    aeq = problem.aeq(:, reverse_permutation);
                    beq = problem.beq;
                case FeatureName.LINEARLY_TRANSFORMED.value
                    % Similar to the case in the custom feature where a custom affine transformation is specified.
                    A = obj.modifier_affine(seed, problem);
                    if isdiag(A)
                        aeq = problem.aeq * A;
                        beq = problem.beq;
                        return;
                    end
                    idx_eq = find(problem.xl == problem.xu);
                    if isempty(problem.aeq)
                        aeq = A(idx_eq, :);
                        beq = problem.xu(idx_eq);
                        return;
                    end
                    aeq = [A(idx_eq, :); problem.aeq * A];
                    beq = [problem.xu(idx_eq); problem.beq];
                otherwise
                    aeq = problem.aeq;
                    beq = problem.beq;
            end
        end

        function f = modifier_fun(obj, x, seed, problem, n_eval)
            %{
            Modify the objective function value.

            Parameters
            ----------
            x : double, size (n,)
                Point at which the objective function is evaluated.
            seed : int
                Seed used to generate random numbers.
            problem : Problem
                Optimization problem to be modified.
            n_eval : int
                Number of evaluations of the objective function.
                (We will use it to generate random streams so that evaluating
                the same point multiple times will not lead to the same
                random numbers.)

            Returns
            -------
            f : double
                Modified objective function value.
            %}

            % Convert x into a cell array. We will later use it to generate 
            % random streams so that randomness of each point is independent.
            xCell = num2cell(x);

            f = problem.fun(x);

            switch obj.name
                case FeatureName.CUSTOM.value
                    if isfield(obj.options, FeatureOptionKey.MOD_FUN.value)
                        rand_stream_custom = obj.default_rng(seed, f, xCell{:}, n_eval);
                        f = obj.options.(FeatureOptionKey.MOD_FUN.value)(x, rand_stream_custom, problem);
                        return;
                    end
                case FeatureName.NOISY.value
                    rand_stream_noisy = obj.default_rng(seed, f, xCell{:}, n_eval);
                    if strcmp(obj.options.(FeatureOptionKey.DISTRIBUTION.value), 'gaussian')
                        noise = rand_stream_noisy.randn();
                    elseif strcmp(obj.options.(FeatureOptionKey.DISTRIBUTION.value), 'uniform')
                        noise = 2 * rand_stream_noisy.rand() - 1;
                    else
                        noise = obj.options.(FeatureOptionKey.DISTRIBUTION.value)(rand_stream_noisy, 1);
                    end
                    if strcmp(obj.options.(FeatureOptionKey.NOISE_TYPE.value), 'absolute')
                        f = f + obj.options.(FeatureOptionKey.NOISE_LEVEL.value) * noise;
                    elseif strcmp(obj.options.(FeatureOptionKey.NOISE_TYPE.value), 'relative')
                        f = f * (1.0 + obj.options.(FeatureOptionKey.NOISE_LEVEL.value) * noise);
                    else
                        % We need the distribution of the noise to be symmetric with respect to 0.
                        f = f + max(1, abs(f)) * obj.options.(FeatureOptionKey.NOISE_LEVEL.value) * noise;
                    end
                case FeatureName.RANDOM_NAN.value
                    rand_stream_random_nan = obj.default_rng(seed, f, xCell{:}, n_eval);
                    if rand_stream_random_nan.rand() < obj.options.(FeatureOptionKey.NAN_RATE.value)
                        f = NaN;
                    end
                case FeatureName.TRUNCATED.value
                    if isnan(f) || isinf(f)
                        % If f is NaN or Inf, we do not need to truncate it.
                        % Note that if f is NaN or Inf, digits will be set to NaN or Inf respectively, which will lead
                        % to an error when calling 'round(f, digits)'.
                        return;
                    end
                    rand_stream_truncated = obj.default_rng(seed, f, xCell{:}, n_eval);
                    if f == 0
                        digits = obj.options.(FeatureOptionKey.SIGNIFICANT_DIGITS.value) - 1;
                    else
                        digits = obj.options.(FeatureOptionKey.SIGNIFICANT_DIGITS.value) - floor(log10(abs(f))) - 1;
                    end
                    f = round(f, digits);
                    % Round f to the desired number of significant digits. (We can also use 'fix(x)'
                    % to round towards zero.) We can shorten the above code by using 'round(f, digits, "significant")'
                    % directly, but we want to keep the same format as the Python code. The code before
                    % 'round(f, digits)' cannot be removed since the default choice of 'round' in MATLAB
                    % will consider digits in relation to the decimal point.
                    if obj.options.(FeatureOptionKey.PERTURBED_TRAILING_DIGITS.value)
                        if f >= 0
                            f = f + rand_stream_truncated.rand() * 10 ^ (-digits);
                        else
                            f = f - rand_stream_truncated.rand() * 10 ^ (-digits);
                        end
                    end
                case FeatureName.UNRELAXABLE_CONSTRAINTS.value
                    [~, maxcv_bounds, maxcv_linear, maxcv_nonlinear] = problem.maxcv(x, true);
                    if obj.options.(FeatureOptionKey.UNRELAXABLE_BOUNDS.value) && maxcv_bounds > 0
                        f = Inf;
                    elseif obj.options.(FeatureOptionKey.UNRELAXABLE_LINEAR_CONSTRAINTS.value) && maxcv_linear > 0
                        f = Inf;
                    elseif obj.options.(FeatureOptionKey.UNRELAXABLE_NONLINEAR_CONSTRAINTS.value) && maxcv_nonlinear > 0
                        f = Inf;
                    end
                case FeatureName.QUANTIZED.value
                    mesh_size = obj.options.(FeatureOptionKey.MESH_SIZE.value);
                    if strcmp(obj.options.(FeatureOptionKey.MESH_TYPE.value), 'relative')
                        mesh_size = mesh_size .* max(1, abs(x));
                    end
                    x = mesh_size .* round(x ./ mesh_size);  % round(x) rounds x to the closest integer.
                    f = problem.fun(x);
                otherwise
                    % Do nothing
            end
        end

        function cub_ = modifier_cub(obj, x, seed, problem, n_eval_cub)
            %{
            Modify the values of the nonlinear inequality constraints.

            Parameters
            ----------
            x : double, size (n,)
                Point at which the nonlinear inequality constraints are evaluated.
            seed : int
                Seed used to generate random numbers.
            problem : Problem
                Problem for which the nonlinear inequality constraints are modified.
            n_eval_cub : int
                Number of evaluations of the nonlinear inequality constraints.
                (We will use it to generate random streams so that evaluating
                the same point multiple times will not lead to the same
                random numbers.)
            
            Returns
            -------
            cub_ : double, size (m_nonlinear_ub,)
                Modified values of the nonlinear inequality constraints.
            %}
            
            % Convert x into a cell array. We will later use it to generate 
            % random streams so that randomness of each point is independent.
            xCell = num2cell(x);

            cub_ = problem.cub(x);
            % If cub_ is empty, return directly!
            if isempty(cub_)
                return;
            end
            cubCell = num2cell(cub_);

            switch obj.name
                case FeatureName.CUSTOM.value
                    if isfield(obj.options, FeatureOptionKey.MOD_CUB.value)
                        rand_stream_custom = obj.default_rng(seed, cubCell{:}, xCell{:}, n_eval_cub);
                        cub_ = obj.options.(FeatureOptionKey.MOD_CUB.value)(x, rand_stream_custom, problem);
                        return;
                    end
                case FeatureName.NOISY.value
                    % Similar to the case in the modifier_fun method.
                    rand_stream_noisy = obj.default_rng(seed, cubCell{:}, xCell{:}, n_eval_cub);
                    if strcmp(obj.options.(FeatureOptionKey.DISTRIBUTION.value), 'gaussian')
                        noise = rand_stream_noisy.randn(size(cub_));
                    elseif strcmp(obj.options.(FeatureOptionKey.DISTRIBUTION.value), 'uniform')
                        noise = 2 * rand_stream_noisy.rand(size(cub_)) - 1;
                    else
                        noise = obj.options.(FeatureOptionKey.DISTRIBUTION.value)(rand_stream_noisy, size(cub_));
                    end
                    if strcmp(obj.options.(FeatureOptionKey.NOISE_TYPE.value), 'absolute')
                        cub_ = cub_ + obj.options.(FeatureOptionKey.NOISE_LEVEL.value) * noise;
                    elseif strcmp(obj.options.(FeatureOptionKey.NOISE_TYPE.value), 'relative')
                        cub_ = cub_ .* (1.0 + obj.options.(FeatureOptionKey.NOISE_LEVEL.value) * noise);
                    else
                        cub_ = cub_ + max(1, abs(cub_)) .* (obj.options.(FeatureOptionKey.NOISE_LEVEL.value) * noise);
                    end
                case FeatureName.RANDOM_NAN.value
                    % Similar to the case in the modifier_fun method.
                    rand_stream_random_nan = obj.default_rng(seed, cubCell{:}, xCell{:}, n_eval_cub);
                    cub_(rand_stream_random_nan.rand(size(cub_)) < obj.options.(FeatureOptionKey.NAN_RATE.value)) = NaN;
                case FeatureName.TRUNCATED.value
                    % Similar to the case in the modifier_fun method.
                    rand_stream_truncated = obj.default_rng(seed, cubCell{:}, xCell{:}, n_eval_cub);
                    digits = zeros(size(cub_));
                    digits(cub_ == 0) = obj.options.(FeatureOptionKey.SIGNIFICANT_DIGITS.value) - 1;
                    digits(cub_ ~= 0) = obj.options.(FeatureOptionKey.SIGNIFICANT_DIGITS.value) - floor(log10(abs(cub_(cub_ ~= 0)))) - 1;
                    % Note that `floor(log10(abs(NaN))) = NaN` and `floor(log10(abs(Inf))) = Inf` in MATLAB.
                    for i_cub = 1:numel(cub_)
                        if ~isnan(cub_(i_cub)) && ~isinf(cub_(i_cub))
                            cub_(i_cub) = round(cub_(i_cub), digits(i_cub));
                        end
                    end
                    if obj.options.(FeatureOptionKey.PERTURBED_TRAILING_DIGITS.value)
                        cub_(cub_ >= 0) = cub_(cub_ >= 0) + rand_stream_truncated.rand(size(cub_(cub_ >= 0))) .* (10 .^ (-digits(cub_ >= 0)));
                        cub_(cub_ < 0) = cub_(cub_ < 0) - rand_stream_truncated.rand(size(cub_(cub_ < 0))) .* (10 .^ (-digits(cub_ < 0)));
                    end
                case FeatureName.NONQUANTIFIABLE_CONSTRAINTS.value
                    % Set the elements whose value are less than or equal to 0 to 0.
                    cub_(cub_ <= 0) = 0;
                    % Set the rest to 1.
                    cub_(~(cub_ <= 0)) = 1;
                case FeatureName.QUANTIZED.value
                    % Similar to the case in the modifier_fun method.
                    mesh_size = obj.options.(FeatureOptionKey.MESH_SIZE.value);
                    if strcmp(obj.options.(FeatureOptionKey.MESH_TYPE.value), 'relative')
                        mesh_size = mesh_size .* max(1, abs(x));
                    end
                    x = mesh_size .* round(x ./ mesh_size);
                    cub_ = problem.cub(x);
                otherwise
                    % Do nothing
            end
        end

        function ceq_ = modifier_ceq(obj, x, seed, problem, n_eval_ceq)
            %{
            Modify the values of the nonlinear equality constraints.

            Parameters
            ----------
            x : double, size (n,)
                Point at which the nonlinear equality constraints are evaluated.
            seed : int
                Seed used to generate random numbers.
            problem : Problem
                Problem for which the nonlinear equality constraints are modified.
            n_eval_ceq : int
                Number of evaluations of the nonlinear equality constraints.
                (We will use it to generate random streams so that evaluating
                the same point multiple times will not lead to the same
                random numbers.)
            
            Returns
            -------
            ceq_ : double, size (m_nonlinear_eq,)
                Modified values of the nonlinear equality constraints.
            %}

            % Convert x into a cell array. We will later use it to generate
            % random streams so that randomness of each point is independent.
            xCell = num2cell(x);

            ceq_ = problem.ceq(x);
            % If ceq_ is empty, return directly!
            if isempty(ceq_)
                return;
            end
            ceqCell = num2cell(ceq_);

            switch obj.name
                case FeatureName.CUSTOM.value
                    if isfield(obj.options, FeatureOptionKey.MOD_CEQ.value)
                        rand_stream_custom = obj.default_rng(seed, ceqCell{:}, xCell{:}, n_eval_ceq);
                        ceq_ = obj.options.(FeatureOptionKey.MOD_CEQ.value)(x, rand_stream_custom, problem);
                        return;
                    end
                case FeatureName.NOISY.value
                    % Similar to the case in the modifier_fun method.
                    rand_stream_noisy = obj.default_rng(seed, ceqCell{:}, xCell{:}, n_eval_ceq);
                    if strcmp(obj.options.(FeatureOptionKey.DISTRIBUTION.value), 'gaussian')
                        noise = rand_stream_noisy.randn(size(ceq_));
                    elseif strcmp(obj.options.(FeatureOptionKey.DISTRIBUTION.value), 'uniform')
                        noise = 2 * rand_stream_noisy.rand(size(ceq_)) - 1;
                    else
                        noise = obj.options.(FeatureOptionKey.DISTRIBUTION.value)(rand_stream_noisy, size(ceq_));
                    end
                    if strcmp(obj.options.(FeatureOptionKey.NOISE_TYPE.value), 'absolute')
                        ceq_ = ceq_ + obj.options.(FeatureOptionKey.NOISE_LEVEL.value) * noise;
                    elseif strcmp(obj.options.(FeatureOptionKey.NOISE_TYPE.value), 'relative')
                        ceq_ = ceq_ .* (1.0 + obj.options.(FeatureOptionKey.NOISE_LEVEL.value) * noise);
                    else
                        ceq_ = ceq_ + max(1, abs(ceq_)) .* (obj.options.(FeatureOptionKey.NOISE_LEVEL.value) * noise);
                    end
                case FeatureName.RANDOM_NAN.value
                    % Similar to the case in the modifier_fun method.
                    rand_stream_random_nan = obj.default_rng(seed, ceqCell{:}, xCell{:}, n_eval_ceq);
                    ceq_(rand_stream_random_nan.rand(size(ceq_)) < obj.options.(FeatureOptionKey.NAN_RATE.value)) = NaN;
                case FeatureName.TRUNCATED.value
                    % Similar to the case in the modifier_fun method.
                    rand_stream_truncated = obj.default_rng(seed, ceqCell{:}, xCell{:}, n_eval_ceq);
                    digits = zeros(size(ceq_));
                    digits(ceq_ == 0) = obj.options.(FeatureOptionKey.SIGNIFICANT_DIGITS.value) - 1;
                    digits(ceq_ ~= 0) = obj.options.(FeatureOptionKey.SIGNIFICANT_DIGITS.value) - floor(log10(abs(ceq_(ceq_ ~= 0)))) - 1;
                    % Note that `floor(log10(abs(NaN))) = NaN` and `floor(log10(abs(Inf))) = Inf` in MATLAB.
                    for i_ceq = 1:numel(ceq_)
                        if ~isnan(ceq_(i_ceq)) && ~isinf(ceq_(i_ceq))
                            ceq_(i_ceq) = round(ceq_(i_ceq), digits(i_ceq));
                        end
                    end
                    if obj.options.(FeatureOptionKey.PERTURBED_TRAILING_DIGITS.value)
                        ceq_(ceq_ >= 0) = ceq_(ceq_ >= 0) + rand_stream_truncated.rand(size(ceq_(ceq_ >= 0))) .* (10 .^ (-digits(ceq_ >= 0)));
                        ceq_(ceq_ < 0) = ceq_(ceq_ < 0) - rand_stream_truncated.rand(size(ceq_(ceq_ < 0))) .* (10 .^ (-digits(ceq_ < 0)));
                    end
                case FeatureName.NONQUANTIFIABLE_CONSTRAINTS.value
                    % Set the elements whose absolute value are less than or equal to 10^(-6) to 0.
                    ceq_(abs(ceq_) <= 1e-6) = 0;
                    % Set the rest to 1.
                    ceq_(~(abs(ceq_) <= 1e-6)) = 1;
                case FeatureName.QUANTIZED.value
                    % Similar to the case in the modifier_fun method.
                    mesh_size = obj.options.(FeatureOptionKey.MESH_SIZE.value);
                    if strcmp(obj.options.(FeatureOptionKey.MESH_TYPE.value), 'relative')
                        mesh_size = mesh_size .* max(1, abs(x));
                    end
                    x = mesh_size .* round(x ./ mesh_size);
                    ceq_ = problem.ceq(x);
                otherwise
                    % Do nothing
            end
        end

        function obj = set_default_options(obj)
            % Set default options.
            switch obj.name
                case FeatureName.PLAIN.value
                    if ~isfield(obj.options, FeatureOptionKey.N_RUNS.value)
                        obj.options.(FeatureOptionKey.N_RUNS.value) = 1;
                    end
                case FeatureName.CUSTOM.value
                    if ~isfield(obj.options, FeatureOptionKey.N_RUNS.value)
                        obj.options.(FeatureOptionKey.N_RUNS.value) = 1;
                    end
                case FeatureName.NOISY.value
                    if ~isfield(obj.options, FeatureOptionKey.DISTRIBUTION.value)
                        obj.options.(FeatureOptionKey.DISTRIBUTION.value) = 'gaussian';
                    end
                    if ~isfield(obj.options, FeatureOptionKey.N_RUNS.value)
                        obj.options.(FeatureOptionKey.N_RUNS.value) = 5;
                    end
                    if ~isfield(obj.options, FeatureOptionKey.NOISE_LEVEL.value)
                        obj.options.(FeatureOptionKey.NOISE_LEVEL.value) = 1e-3;
                    end
                    if ~isfield(obj.options, FeatureOptionKey.NOISE_TYPE.value)
                        obj.options.(FeatureOptionKey.NOISE_TYPE.value) = 'mixed';
                    end
                case FeatureName.PERMUTED.value
                    if ~isfield(obj.options, FeatureOptionKey.N_RUNS.value)
                        obj.options.(FeatureOptionKey.N_RUNS.value) = 5;
                    end
                case FeatureName.LINEARLY_TRANSFORMED.value
                    if ~isfield(obj.options, FeatureOptionKey.N_RUNS.value)
                        obj.options.(FeatureOptionKey.N_RUNS.value) = 5;
                    end
                    if ~isfield(obj.options, FeatureOptionKey.ROTATED.value)
                        obj.options.(FeatureOptionKey.ROTATED.value) = true;
                    end
                    if ~isfield(obj.options, FeatureOptionKey.CONDITION_FACTOR.value)
                        obj.options.(FeatureOptionKey.CONDITION_FACTOR.value) = 0;
                    end
                case FeatureName.PERTURBED_X0.value
                    if ~isfield(obj.options, FeatureOptionKey.DISTRIBUTION.value)
                        obj.options.(FeatureOptionKey.DISTRIBUTION.value) = 'spherical';
                    end
                    if ~isfield(obj.options, FeatureOptionKey.PERTURBATION_LEVEL.value)
                        obj.options.(FeatureOptionKey.PERTURBATION_LEVEL.value) = 1e-3;
                    end
                    if ~isfield(obj.options, FeatureOptionKey.N_RUNS.value)
                        obj.options.(FeatureOptionKey.N_RUNS.value) = 5;
                    end
                case FeatureName.RANDOM_NAN.value
                    if ~isfield(obj.options, FeatureOptionKey.N_RUNS.value)
                        obj.options.(FeatureOptionKey.N_RUNS.value) = 5;
                    end
                    if ~isfield(obj.options, FeatureOptionKey.NAN_RATE.value)
                        obj.options.(FeatureOptionKey.NAN_RATE.value) = 0.05;
                    end
                case FeatureName.TRUNCATED.value
                    if ~isfield(obj.options, FeatureOptionKey.PERTURBED_TRAILING_DIGITS.value)
                        obj.options.(FeatureOptionKey.PERTURBED_TRAILING_DIGITS.value) = false;
                    end
                    if ~isfield(obj.options, FeatureOptionKey.N_RUNS.value)
                        if obj.options.(FeatureOptionKey.PERTURBED_TRAILING_DIGITS.value)
                            obj.options.(FeatureOptionKey.N_RUNS.value) = 5;
                        else
                            obj.options.(FeatureOptionKey.N_RUNS.value) = 1;
                        end
                    end
                    if ~isfield(obj.options, FeatureOptionKey.SIGNIFICANT_DIGITS.value)
                        obj.options.(FeatureOptionKey.SIGNIFICANT_DIGITS.value) = 6;
                    end
                case FeatureName.UNRELAXABLE_CONSTRAINTS.value
                    if ~isfield(obj.options, FeatureOptionKey.N_RUNS.value)
                        obj.options.(FeatureOptionKey.N_RUNS.value) = 1;
                    end
                    if ~isfield(obj.options, FeatureOptionKey.UNRELAXABLE_BOUNDS.value)
                        obj.options.(FeatureOptionKey.UNRELAXABLE_BOUNDS.value) = true;
                    end
                    if ~isfield(obj.options, FeatureOptionKey.UNRELAXABLE_LINEAR_CONSTRAINTS.value)
                        obj.options.(FeatureOptionKey.UNRELAXABLE_LINEAR_CONSTRAINTS.value) = false;
                    end
                    if ~isfield(obj.options, FeatureOptionKey.UNRELAXABLE_NONLINEAR_CONSTRAINTS.value)
                        obj.options.(FeatureOptionKey.UNRELAXABLE_NONLINEAR_CONSTRAINTS.value) = false;
                    end
                case FeatureName.NONQUANTIFIABLE_CONSTRAINTS.value
                    if ~isfield(obj.options, FeatureOptionKey.N_RUNS.value)
                        obj.options.(FeatureOptionKey.N_RUNS.value) = 1;
                    end
                case FeatureName.QUANTIZED.value
                    if ~isfield(obj.options, FeatureOptionKey.N_RUNS.value)
                        obj.options.(FeatureOptionKey.N_RUNS.value) = 1;
                    end
                    if ~isfield(obj.options, FeatureOptionKey.MESH_SIZE.value)
                        obj.options.(FeatureOptionKey.MESH_SIZE.value) = 1e-3;
                    end
                    if ~isfield(obj.options, FeatureOptionKey.MESH_TYPE.value)
                        obj.options.(FeatureOptionKey.MESH_TYPE.value) = 'absolute';
                    end
                    if ~isfield(obj.options, FeatureOptionKey.GROUND_TRUTH.value)
                        obj.options.(FeatureOptionKey.GROUND_TRUTH.value) = true;
                    end
                otherwise
                    error("MATLAB:Feature:UnknownFeature", "Unknown feature: " + obj.name + ".")
            end
        end

    end

    methods (Static)
        function rand_stream = default_rng(seed, varargin)
            % Generate a random number generator.
            %
            % Parameters
            % ----------
            % seed : double, but default one is 'shuffle'
            %     Seed used to generate an initial random number generator.
            % varargin : array of double
            %     Arguments used to generate the returned random number generator.
            %
            % Returns
            % -------
            % rand_stream : RandStream
            %     Random number generator.

            % Create an initial rand_stream with the given seed
            if nargin < 1 || isempty(seed)
                seed = 'shuffle';  % Default behavior is to shuffle the seed.
            end
            if ~strcmp(seed, 'shuffle')
                if isnan(seed) || isinf(seed)
                    seed = 0;
                end
                if ~isrealscalar(seed)
                    error("MATLAB:Feature:SeedNotEvenReal", "The input `seed` for method `default_rng` in `Feature` should be a real number.")
                end
                seed = mod(floor(seed), 2^32);
            end

            % Convert all elements in varargin to double
            varargin = cellfun(@double, varargin, 'UniformOutput', false);

            % Create a RandStream object with the specified seed
            rand_stream = RandStream('mt19937ar', 'Seed', seed);

            % Generate a new seed based on the initial rand_stream and the additional arguments
            newSeed = abs(sin(1e5 * randn(rand_stream, 1)) + sum(sin(1e5 * prod(cellfun(@(x) x, varargin))))) * 1e9;
            newSeed = mod(floor(newSeed), 2^32);
            if isnan(newSeed) || isinf(newSeed)
                newSeed = 0;
            end
            rand_stream = RandStream('mt19937ar', 'Seed', floor(newSeed));
        end
    end
end