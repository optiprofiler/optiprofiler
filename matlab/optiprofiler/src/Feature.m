classdef Feature < handle
%FEATURE maps an optimization problem to a new one with special properties.
%
%   We are interested to test solvers under different oracles. For example, we
%   want to test the performance of solvers under the case where the objective
%   function is noisy. For this purpose, we define FEATURE class. ``FEATURE''
%   means the oracle or the testing environment. Mathematically speaking, a
%   FEATURE is a mapping from an optimization problem to a new one.
%
%   Suppose we have an optimization problem
%    min f(x)
%    s.t. xl <= x <= xu
%         aub @ x <= bub
%         aeq @ x = beq
%         cub(x) <= 0
%         ceq(x) = 0
%   with initial point x0.
%   Then, FEATURE maps the optimization problem to the following one.
%    min f_mod(A @ x + b)
%    s.t. xl_mod <= A @ x + b <= xu_mod
%         aub_mod @ (A @ x + b) <= bub_mod
%         aeq_mod @ (A @ x + b) = beq_mod
%         cub_mod(A @ x + b) <= 0
%         ceq_mod(A @ x + b) = 0
%   with initial point x0_mod.
%
%Usage:
%   FT = FEATURE(NAME) creates a FEATURE with the name NAME.
%
%   FT = FEATURE(NAME, OPTIONS) creates a FEATURE with the name NAME and the
%   options OPTIONS.
%
%FEATURE properties:
%   NAME: str
%       Name of the FEATURE. The legal values are 'plain', 'perturbed_x0',
%       'noisy', 'truncated', 'permuted', 'linearly_transformed', 'random_nan',
%       'unrelaxable_constraints', 'custom'. The default value is 'plain'. We
%       will introduce the meaning of each FEATURE later in Details part.
%   OPTIONS: struct
%       Options of the FEATURE. With a given NAME, we can specify different
%       legal OPTIONS. We will introduce the OPTIONS of each FEATURE later in
%       Details part.
%
%FEATURE methods: we give 8 freedom to define a FEATURE by following methods.
%   1. modifier_x0:
%       a function handle to modify the initial point.
%   2. modifier_affine:
%       a function handle to generate an invertible matrix A (and its inverse)
%       and a vector b for the affine transformation applied to the variables.
%   3. modifier_bounds:
%       a function handle to modify the lower and upper bounds.
%   4. modifier_linear_ub:
%       a function handle to modify the linear inequality constraints.
%   5. modifier_linear_eq:
%       a function handle to modify the linear equality constraints.
%   6. modifier_fun:
%       a function handle to modify the objective function value.
%   7. modifier_cub:
%       a function handle to modify the values of the nonlinear inequality
%       constraints.
%   8. modifier_ceq:
%       a function handle to modify the values of the nonlinear equality
%       constraints.
%
%Details: we provide some predefined features as follows.
%   1. plain:
%       do nothing.
%   2. perturbed_x0:
%       perturb the initial point.
%   3. noisy:
%       add noise to the objective function value.
%   4. truncated:
%       truncate values of the objective function and nonlinearity constraints
%       to a given number of significant digits.
%   5. permuted:
%       randomly permute the variables. Note that x0, xl, xu, aub, bub, aeq,
%       beq should be modified since we need to keep the new problem
%       mathematically equivalent to the original one.
%   6. linearly_transformed:
%       generate a invertible linear transformation in the form D * Q' with D
%       being a diagonal matrix and Q being an orthogonal matrix, and apply the
%       transformation to the variables. In this way, the Hessian of the
%       objective function becomes Q * D * H * D * Q' with H being the original
%       Hessian. Note that x0, xl, xu, aub, bub, aeq, beq should be modified.
%   7. random_nan:
%       replace values of the objective function and nonlinear constraints with
%       NaN with a given rate.
%   8. unrelaxable_constraints:
%       set the objective function value to Inf outside the feasible region.
%   9. nonquantifiable_constraints:
%       replace values of nonlinear constraints with either 0 (if the
%       constraint is satisfied) or 1 (if the constraint is violated).
%   10. custom:
%       a user-defined FEATURE. Users can define their own FEATURE by
%       specifying any modifier listed in the above.



    properties (GetAccess = public, SetAccess = private)

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

            if ~ischarstr(name)
                error("MATLAB:Feature:FeaturenameNotString", "The feature name must be a char or a string.")
            end
            % Convert the feature name to lowercase characters.
            obj.name = char(lower(name));

            obj.options = struct();
            if nargin > 1
                if isstruct(varargin{1})
                    obj.options = varargin{1};
                elseif mod(length(varargin), 2) ~= 0
                    error("MATLAB:Feature:InvalidNumberOfArguments", "The input of the feature options must be in pairs.")
                else
                    for i = 1:2:length(varargin)
                        obj.options.(lower(varargin{i})) = varargin{i+1};
                    end
                end
            end

            % Check whether the feature is valid.
            % validFeatureNames = {enumeration('FeatureName').value};
            % Only for MATLAB R2021b and later.
            validFeatureNames = cellfun(@(x) x.value, num2cell(enumeration('FeatureName')), 'UniformOutput', false);
            if ~ismember(obj.name, validFeatureNames)
                error("MATLAB:Feature:UnknownFeature", "Unknown feature: " + obj.name + ".")
            end
            optionKeys = fieldnames(obj.options);
            validOptionKeys = cellfun(@(x) x.value, num2cell(enumeration('FeatureOptionKey')), 'UniformOutput', false);
            for i = 1:numel(optionKeys)
                key = optionKeys{i};
                if ~ismember(key, validOptionKeys)
                    error("MATLAB:Feature:UnknownOption", "Unknown option: " + key + ".")
                end

                % Check whether the option is valid for the feature.
                known_options = {FeatureOptionKey.N_RUNS.value};
                switch obj.name
                    case FeatureName.CUSTOM.value
                        known_options = [known_options, {FeatureOptionKey.MOD_X0.value, FeatureOptionKey.MOD_BOUNDS.value, FeatureOptionKey.MOD_LINEAR_UB.value, FeatureOptionKey.MOD_LINEAR_EQ.value, FeatureOptionKey.MOD_AFFINE.value, FeatureOptionKey.MOD_FUN.value, FeatureOptionKey.MOD_CUB.value, FeatureOptionKey.MOD_CEQ.value}];
                    case FeatureName.NOISY.value
                        known_options = [known_options, {FeatureOptionKey.DISTRIBUTION.value, FeatureOptionKey.NOISE_LEVEL.value, FeatureOptionKey.NOISE_TYPE.value}];
                    case FeatureName.PERTURBED_X0.value
                        known_options = [known_options, {FeatureOptionKey.DISTRIBUTION.value, FeatureOptionKey.NOISE_LEVEL.value}];
                    case FeatureName.RANDOM_NAN.value
                        known_options = [known_options, {FeatureOptionKey.RATE_NAN.value}];
                    case FeatureName.TRUNCATED.value
                        known_options = [known_options, {FeatureOptionKey.PERTURBED_TRAILING_ZEROS.value, FeatureOptionKey.SIGNIFICANT_DIGITS.value}];
                    case FeatureName.UNRELAXABLE_CONSTRAINTS.value
                        known_options = [known_options, {FeatureOptionKey.UNRELAXABLE_BOUNDS.value, FeatureOptionKey.UNRELAXABLE_LINEAR_CONSTRAINTS.value, FeatureOptionKey.UNRELAXABLE_NONLINEAR_CONSTRAINTS.value}];
                    case FeatureName.LINEARLY_TRANSFORMED.value
                        known_options = [known_options, {FeatureOptionKey.ROTATED.value, FeatureOptionKey.CONDITION_FACTOR.value}];
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
                    error("MATLAB:Feature:InvalidOptionForFeature", "Option " + key + " is not valid for feature " + obj.name + ".")
                end

                % Check whether the option type is valid.
                switch key
                    case FeatureOptionKey.N_RUNS.value
                        if ~isintegerscalar(obj.options.(key)) || obj.options.(key) <= 0
                            error("MATLAB:Feature:n_runs_NotPositiveInteger", "Option " + key + " must be a positive integer.")
                        end
                    case FeatureOptionKey.DISTRIBUTION.value
                        if ~isa(obj.options.(key), 'function_handle')
                            error("MATLAB:Feature:distribution_NotFunctionHandle", "Option " + key + " must be a function handle.")
                        end
                    case FeatureOptionKey.RATE_NAN.value
                        if ~isrealscalar(obj.options.(key)) || obj.options.(key) < 0.0 || obj.options.(key) > 1.0
                            error("MATLAB:Feature:rate_nan_NotBetween_0_1", "Option " + key + " must be a real number between 0 and 1.")
                        end
                    case FeatureOptionKey.SIGNIFICANT_DIGITS.value
                        if ~isintegerscalar(obj.options.(key)) || obj.options.(key) <= 0
                            error("MATLAB:Feature:significant_digits_NotPositiveInteger", "Option " + key + " must be a positive integer.")
                        end
                    case FeatureOptionKey.NOISE_LEVEL.value
                        if ~isrealscalar(obj.options.(key)) || obj.options.(key) <= 0.0
                            error("MATLAB:Feature:noise_level_NotPositive", "Option " + key + " must be a positive real number.")
                        end
                    case FeatureOptionKey.NOISE_TYPE.value
                        validNoiseTypes = cellfun(@(x) x.value, num2cell(enumeration('NoiseType')), 'UniformOutput', false);
                        if ~ischarstr(obj.options.(key)) || ~ismember(obj.options.(key), validNoiseTypes)
                            error("MATLAB:Feature:noise_type_InvalidInput", "Option " + key + " must be '" + NoiseType.ABSOLUTE.value + "' or '" + NoiseType.RELATIVE.value + "' or '" + NoiseType.MIXED.value + "'.")
                        end
                    case FeatureOptionKey.PERTURBED_TRAILING_ZEROS.value
                        if ~islogicalscalar(obj.options.(key))
                            error("MATLAB:Feature:perturbed_trailing_zeros_NotLogical", "Option " + key + " must be a logical.")
                        end
                    case FeatureOptionKey.ROTATED.value
                        if ~islogicalscalar(obj.options.(key))
                            error("MATLAB:Feature:rotated_NotLogical", "Option " + key + " must be a logical.")
                        end
                    case FeatureOptionKey.CONDITION_FACTOR.value
                        if ~isa(obj.options.(key), 'function_handle')
                            error("MATLAB:Feature:condition_factor_InvalidInput", "Option " + key + " must be a function handle.")
                        end
                    case FeatureOptionKey.UNRELAXABLE_BOUNDS.value
                        if ~islogicalscalar(obj.options.(key))
                            error("MATLAB:Feature:unrelaxable_bounds_NotLogical", "Option " + key + " must be a logical.")
                        end
                    case FeatureOptionKey.UNRELAXABLE_LINEAR_CONSTRAINTS.value
                        if ~islogicalscalar(obj.options.(key))
                            error("MATLAB:Feature:unrelaxable_linear_constraints_NotLogical", "Option " + key + " must be a logical.")
                        end
                    case FeatureOptionKey.UNRELAXABLE_NONLINEAR_CONSTRAINTS.value
                        if ~islogicalscalar(obj.options.(key))
                            error("MATLAB:Feature:unrelaxable_nonlinear_constraints_NotLogical", "Option " + key + " must be a logical.")
                        end
                    case FeatureOptionKey.MOD_X0.value
                        if ~isa(obj.options.(key), 'function_handle')
                            error("MATLAB:Feature:mod_x0_NotFunctionHandle", "Option " + key + " must be a function handle.")
                        end
                    case FeatureOptionKey.MOD_BOUNDS.value
                        if ~isa(obj.options.(key), 'function_handle')
                            error("MATLAB:Feature:mod_bounds_NotFunctionHandle", "Option " + key + " must be a function handle.")
                        end
                    case FeatureOptionKey.MOD_LINEAR_UB.value
                        if ~isa(obj.options.(key), 'function_handle')
                            error("MATLAB:Feature:mod_linear_ub_NotFunctionHandle", "Option " + key + " must be a function handle.")
                        end
                    case FeatureOptionKey.MOD_LINEAR_EQ.value
                        if ~isa(obj.options.(key), 'function_handle')
                            error("MATLAB:Feature:mod_linear_eq_NotFunctionHandle", "Option " + key + " must be a function handle.")
                        end
                    case FeatureOptionKey.MOD_AFFINE.value
                        if ~isa(obj.options.(key), 'function_handle')
                            error("MATLAB:Feature:mod_affine_NotFunctionHandle", "Option " + key + " must be a function handle.")
                        end
                    case FeatureOptionKey.MOD_FUN.value
                        if ~isa(obj.options.(key), 'function_handle')
                            error("MATLAB:Feature:mod_fun_NotFunctionHandle", "Option " + key + " must be a function handle.")
                        end
                    case FeatureOptionKey.MOD_CUB.value
                        if ~isa(obj.options.(key), 'function_handle')
                            error("MATLAB:Feature:mod_cub_NotFunctionHandle", "Option " + key + " must be a function handle.")
                        end
                    case FeatureOptionKey.MOD_CEQ.value
                        if ~isa(obj.options.(key), 'function_handle')
                            error("MATLAB:Feature:mod_ceq_NotFunctionHandle", "Option " + key + " must be a function handle.")
                        end
                end
            end

            % Set default options for the unspecified options.
            obj = set_default_options(obj);
        end

        function is_stochastic = isStochastic(obj)
            stochasticFeatures = {FeatureName.CUSTOM.value, FeatureName.NOISY.value, FeatureName.PERMUTED.value, FeatureName.LINEARLY_TRANSFORMED.value, FeatureName.PERTURBED_X0.value, FeatureName.RANDOM_NAN.value, FeatureName.TRUNCATED.value};
            is_stochastic = ismember(obj.name, stochasticFeatures);
        end

        function x0 = modifier_x0(obj, seed, problem)
            %{
            Modify the initial point.

            Parameters
            ----------
            seed : int
                Seed used to generate random numbers.
            problem : Problem
                Optimization problem to be modified.
            
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
                    rand_stream_perturbed_x0 = obj.default_rng(seed, problem.x0);
                    x0 = problem.x0 + obj.options.(FeatureOptionKey.NOISE_LEVEL.value) * max(1, norm(problem.x0)) * obj.options.(FeatureOptionKey.DISTRIBUTION.value)(rand_stream_perturbed_x0, problem.n);
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
            % Define the affine transformation A @ x + b and the inverse inv of A. The transformation
            % is applied to the variables later.

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
                        error("MATLAB:Feature:AffineTransformationNotInvertible", "The affine transformation A * inv is not an identity matrix.")
                    end
                case FeatureName.PERMUTED.value
                    % Generate a random permutation matrix.
                    rand_stream_permuted = obj.default_rng(seed);
                    permutation = rand_stream_permuted.randperm(problem.n);
                    [~, reverse_permutation] = sort(permutation);
                    A = eye(problem.n);
                    A = A(reverse_permutation, :);
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
            Modify the lower bounds.

            Parameters
            ----------
            seed : int
                Seed used to generate random numbers.
            problem : Problem
                Optimization problem to be modified.

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
                    % If the inverse of the affine transformation is diagonal, we can aplly it to get
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
                Optimization problem to be modified.

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
                    xl <= A @ x + b <= xu
                    should be modified to
                    A @ x <= xu - b
                    and
                    -A @ x <= -xl + b
                    when A is not diagonal.

                    Linear inequality constraints
                    turn out to be
                    (aub @ A) @ x <= bub - aub @ b
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
                Optimization problem to be modified.

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
                    xl <= A @ x + b <= xu
                    with xl = xu should be modified to
                    A @ x = xu - b
                    when A is not diagonal.

                    Linear equality constraints
                    turn out to be
                    (aeq @ A) @ x = beq - aeq @ b
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

        function f = modifier_fun(obj, x, f, seed, problem)
            %{
            Modify the objective function value.

            Parameters
            ----------
            x : double, size (n,)
                Decision variables.
            f : double
                Objective function value.
            seed : int
                Seed used to generate random numbers.
            problem : Problem
                Optimization problem to be modified.

            Returns
            -------
            f : double
                Modified objective function value.
            %}

            % Convert x into a cell array. We will later use it to generate 
            % random streams so that randomness of each point is independent.
            xCell = num2cell(x);

            switch obj.name
                case FeatureName.CUSTOM.value
                    if isfield(obj.options, FeatureOptionKey.MOD_FUN.value)
                        rand_stream_custom = obj.default_rng(seed, f, xCell{:});
                        f = obj.options.(FeatureOptionKey.MOD_FUN.value)(x, f, rand_stream_custom, problem);
                        return;
                    end
                case FeatureName.NOISY.value
                    rand_stream_noisy = obj.default_rng(seed, f, obj.options.(FeatureOptionKey.NOISE_LEVEL.value), sum(double(obj.options.(FeatureOptionKey.NOISE_TYPE.value))), xCell{:});
                    if strcmp(obj.options.(FeatureOptionKey.NOISE_TYPE.value), NoiseType.ABSOLUTE.value)
                        f = f + obj.options.(FeatureOptionKey.NOISE_LEVEL.value) * obj.options.(FeatureOptionKey.DISTRIBUTION.value)(rand_stream_noisy, 1);
                    elseif strcmp(obj.options.(FeatureOptionKey.NOISE_TYPE.value), NoiseType.RELATIVE.value)
                        f = f * (1.0 + obj.options.(FeatureOptionKey.NOISE_LEVEL.value) * obj.options.(FeatureOptionKey.DISTRIBUTION.value)(rand_stream_noisy, 1));
                    else
                        % We need the distribution of the noise to be symmetric with respect to 0.
                        f = f + max(1, abs(f)) * obj.options.(FeatureOptionKey.NOISE_LEVEL.value) * obj.options.(FeatureOptionKey.DISTRIBUTION.value)(rand_stream_noisy, 1);
                    end
                case FeatureName.RANDOM_NAN.value
                    rand_stream_random_nan = obj.default_rng(seed, f, obj.options.(FeatureOptionKey.RATE_NAN.value), xCell{:});
                    if rand_stream_random_nan.rand() < obj.options.(FeatureOptionKey.RATE_NAN.value)
                        f = NaN;
                    end
                case FeatureName.TRUNCATED.value
                    rand_stream_truncated = obj.default_rng(seed, f, obj.options.(FeatureOptionKey.SIGNIFICANT_DIGITS.value), xCell{:});
                    if f == 0
                        digits = obj.options.(FeatureOptionKey.SIGNIFICANT_DIGITS.value) - 1;
                    else
                        digits = obj.options.(FeatureOptionKey.SIGNIFICANT_DIGITS.value) - floor(log10(abs(f))) - 1;
                    end
                    digits = double(digits);
                    f = round(f, digits);
                    % Round f to the desired number of significant digits. (We can also use 'fix(x)'
                    % to round towards zero.) We can shorten the above code by using 'round(f, digits, "significant")'
                    % directly, but we want to keep the same format as the Python code. The code before
                    % 'round(f, digits)' cannot be removed since the default choice of 'round' in MATLAB
                    % will consider digits in relation to the decimal point.
                    if obj.options.(FeatureOptionKey.PERTURBED_TRAILING_ZEROS.value)
                        if f >= 0
                            f = f + rand_stream_truncated.rand() * 10 ^ (-digits);
                        else
                            f = f - rand_stream_truncated.rand() * 10 ^ (-digits);
                        end
                    end
                case FeatureName.UNRELAXABLE_CONSTRAINTS.value
                    [~, maxcv_bounds, maxcv_linear, maxcv_nonlinear] = problem.maxcv(x);
                    if obj.options.(FeatureOptionKey.UNRELAXABLE_BOUNDS.value) && maxcv_bounds > 0
                        f = Inf;
                    elseif obj.options.(FeatureOptionKey.UNRELAXABLE_LINEAR_CONSTRAINTS.value) && maxcv_linear > 0
                        f = Inf;
                    elseif obj.options.(FeatureOptionKey.UNRELAXABLE_NONLINEAR_CONSTRAINTS.value) && maxcv_nonlinear > 0
                        f = Inf;
                    end
                otherwise
                    % Do nothing
            end
        end

        function cub_ = modifier_cub(obj, x, cub_, seed, problem)
            %{
            Modify the values of the nonlinear inequality constraints.

            Parameters
            ----------
            x : double, size (n,)
                Decision variables.
            cub_ : double, size (m_nonlinear_ub,)
                Values of the upper bounds of the nonlinear inequality
                constraints.
            seed : int
                Seed used to generate random numbers.
            problem : Problem
                Optimization problem to be modified.
            
            Returns
            -------
            cub_ : double, size (m_nonlinear_ub,)
                Modified values of the upper bounds of the nonlinear inequality
                constraints.
            %}
            
            % Convert x into a cell array. We will later use it to generate 
            % random streams so that randomness of each point is independent.
            xCell = num2cell(x);

            switch obj.name
                case FeatureName.CUSTOM.value
                    if isfield(obj.options, FeatureOptionKey.MOD_CUB.value)
                        rand_stream_custom = obj.default_rng(seed, cub_, xCell{:});
                        cub_ = obj.options.(FeatureOptionKey.MOD_CUB.value)(x, cub_, rand_stream_custom, problem);
                        return;
                    end
                case FeatureName.NOISY.value
                    % Similar to the case in the modifier_fun method.
                    rand_stream_noisy = obj.default_rng(seed, cub_, obj.options.(FeatureOptionKey.NOISE_LEVEL.value), sum(double(obj.options.(FeatureOptionKey.NOISE_TYPE.value))), xCell{:});
                    if strcmp(obj.options.(FeatureOptionKey.NOISE_TYPE.value), NoiseType.ABSOLUTE.value)
                        cub_ = cub_ + obj.options.(FeatureOptionKey.NOISE_LEVEL.value) * obj.options.(FeatureOptionKey.DISTRIBUTION.value)(rand_stream_noisy, 1);
                    elseif strcmp(obj.options.(FeatureOptionKey.NOISE_TYPE.value), NoiseType.RELATIVE.value)
                        cub_ = cub_ * (1.0 + obj.options.(FeatureOptionKey.NOISE_LEVEL.value) * obj.options.(FeatureOptionKey.DISTRIBUTION.value)(rand_stream_noisy, 1));
                    else
                        cub_ = cub_ + max(1, abs(f)) * obj.options.(FeatureOptionKey.NOISE_LEVEL.value) * obj.options.(FeatureOptionKey.DISTRIBUTION.value)(rand_stream_noisy, 1);
                    end
                case FeatureName.RANDOM_NAN.value
                    % Similar to the case in the modifier_fun method.
                    rand_stream_random_nan = obj.default_rng(seed, cub_, obj.options.(FeatureOptionKey.RATE_NAN.value), xCell{:});
                    if rand_stream_random_nan.rand() < obj.options.(FeatureOptionKey.RATE_NAN.value)
                        cub_ = NaN;
                    end
                case FeatureName.TRUNCATED.value
                    % Similar to the case in the modifier_fun method.
                    rand_stream_truncated = obj.default_rng(seed, cub_, obj.options.(FeatureOptionKey.SIGNIFICANT_DIGITS.value), xCell{:});
                    if cub_ == 0
                        digits = obj.options.(FeatureOptionKey.SIGNIFICANT_DIGITS.value) - 1;
                    else
                        digits = obj.options.(FeatureOptionKey.SIGNIFICANT_DIGITS.value) - floor(log10(abs(f))) - 1;
                    end
                    digits = double(digits);
                    cub_ = round(cub_, digits);
                    if obj.options.(FeatureOptionKey.PERTURBED_TRAILING_ZEROS.value)
                        if cub_ >= 0
                            cub_ = cub_ + rand_stream_truncated.rand() * 10 ^ (-digits);
                        else
                            cub_ = cub_ - rand_stream_truncated.rand() * 10 ^ (-digits);
                        end
                    end
                case FeatureName.NONQUANTIFIABLE_CONSTRAINTS.value
                    % Set the elements whose value are less than or equal to 0 to 0.
                    cub_(cub_ <= 0) = 0;
                    % Set the rest to 1.
                    cub_(~(cub_ <= 0)) = 1;
                otherwise
                    % Do nothing
            end
        end

        function ceq_ = modifier_ceq(obj, x, ceq_, seed, problem)
            %{
            Modify the values of the nonlinear equality constraints.

            Parameters
            ----------
            x : double, size (n,)
                Decision variables.
            ceq_ : double, size (m_nonlinear_eq,)
                Values of the upper bounds of the nonlinear equality
                constraints.
            seed : int
                Seed used to generate random numbers.
            problem : Problem
                Optimization problem to be modified.
            
            Returns
            -------
            ceq_ : double, size (m_nonlinear_eq,)
                Modified values of the upper bounds of the nonlinear equality
                constraints.
            %}

            % Convert x into a cell array. We will later use it to generate
            % random streams so that randomness of each point is independent.
            xCell = num2cell(x);

            switch obj.name
                case FeatureName.CUSTOM.value
                    if isfield(obj.options, FeatureOptionKey.MOD_CEQ.value)
                        rand_stream_custom = obj.default_rng(seed, ceq_, xCell{:});
                        ceq_ = obj.options.(FeatureOptionKey.MOD_CEQ.value)(x, ceq_, rand_stream_custom, problem);
                        return;
                    end
                case FeatureName.NOISY.value
                    % Similar to the case in the modifier_fun method.
                    rand_stream_noisy = obj.default_rng(seed, ceq_, obj.options.(FeatureOptionKey.NOISE_LEVEL.value), sum(double(obj.options.(FeatureOptionKey.NOISE_TYPE.value))), xCell{:});
                    if strcmp(obj.options.(FeatureOptionKey.NOISE_TYPE.value), NoiseType.ABSOLUTE.value)
                        ceq_ = ceq_ + obj.options.(FeatureOptionKey.NOISE_LEVEL.value) * obj.options.(FeatureOptionKey.DISTRIBUTION.value)(rand_stream_noisy, 1);
                    elseif strcmp(obj.options.(FeatureOptionKey.NOISE_TYPE.value), NoiseType.RELATIVE.value)
                        ceq_ = ceq_ * (1.0 + obj.options.(FeatureOptionKey.NOISE_LEVEL.value) * obj.options.(FeatureOptionKey.DISTRIBUTION.value)(rand_stream_noisy, 1));
                    else
                        ceq_ = ceq_ + max(1, abs(f)) * obj.options.(FeatureOptionKey.NOISE_LEVEL.value) * obj.options.(FeatureOptionKey.DISTRIBUTION.value)(rand_stream_noisy, 1);
                    end
                case FeatureName.RANDOM_NAN.value
                    % Similar to the case in the modifier_fun method.
                    rand_stream_random_nan = obj.default_rng(seed, ceq_, obj.options.(FeatureOptionKey.RATE_NAN.value), xCell{:});
                    if rand_stream_random_nan.rand() < obj.options.(FeatureOptionKey.RATE_NAN.value)
                        ceq_ = NaN;
                    end
                case FeatureName.TRUNCATED.value
                    % Similar to the case in the modifier_fun method.
                    rand_stream_truncated = obj.default_rng(seed, ceq_, obj.options.(FeatureOptionKey.SIGNIFICANT_DIGITS.value), xCell{:});
                    if ceq_ == 0
                        digits = obj.options.(FeatureOptionKey.SIGNIFICANT_DIGITS.value) - 1;
                    else
                        digits = obj.options.(FeatureOptionKey.SIGNIFICANT_DIGITS.value) - floor(log10(abs(f))) - 1;
                    end
                    digits = double(digits);
                    ceq_ = round(ceq_, digits);
                    if obj.options.(FeatureOptionKey.PERTURBED_TRAILING_ZEROS.value)
                        if ceq_ >= 0
                            ceq_ = ceq_ + rand_stream_truncated.rand() * 10 ^ (-digits);
                        else
                            ceq_ = ceq_ - rand_stream_truncated.rand() * 10 ^ (-digits);
                        end
                    end
                case FeatureName.NONQUANTIFIABLE_CONSTRAINTS.value
                    % Set the elements whose absolute value are less than or equal to 10^(-6) to 0.
                    ceq_(abs(ceq_) <= 1e-6) = 0;
                    % Set the rest to 1.
                    ceq_(~(abs(ceq_) <= 1e-6)) = 1;
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
                        obj.options.(FeatureOptionKey.DISTRIBUTION.value) = @obj.default_distribution;
                    end
                    if ~isfield(obj.options, FeatureOptionKey.N_RUNS.value)
                        obj.options.(FeatureOptionKey.N_RUNS.value) = 10;
                    end
                    if ~isfield(obj.options, FeatureOptionKey.NOISE_LEVEL.value)
                        obj.options.(FeatureOptionKey.NOISE_LEVEL.value) = 1e-3;
                    end
                    if ~isfield(obj.options, FeatureOptionKey.NOISE_TYPE.value)
                        obj.options.(FeatureOptionKey.NOISE_TYPE.value) = NoiseType.MIXED.value;
                    end
                case FeatureName.PERMUTED.value
                    if ~isfield(obj.options, FeatureOptionKey.N_RUNS.value)
                        obj.options.(FeatureOptionKey.N_RUNS.value) = 10;
                    end
                case FeatureName.LINEARLY_TRANSFORMED.value
                    if ~isfield(obj.options, FeatureOptionKey.N_RUNS.value)
                        obj.options.(FeatureOptionKey.N_RUNS.value) = 10;
                    end
                    if ~isfield(obj.options, FeatureOptionKey.ROTATED.value)
                        obj.options.(FeatureOptionKey.ROTATED.value) = true;
                    end
                    if ~isfield(obj.options, FeatureOptionKey.CONDITION_FACTOR.value)
                        obj.options.(FeatureOptionKey.CONDITION_FACTOR.value) = 0;
                    end
                case FeatureName.PERTURBED_X0.value
                    if ~isfield(obj.options, FeatureOptionKey.DISTRIBUTION.value)
                        obj.options.(FeatureOptionKey.DISTRIBUTION.value) = @obj.default_distribution;
                    end
                    if ~isfield(obj.options, FeatureOptionKey.NOISE_LEVEL.value)
                        obj.options.(FeatureOptionKey.NOISE_LEVEL.value) = 1e-3;
                    end
                    if ~isfield(obj.options, FeatureOptionKey.N_RUNS.value)
                        obj.options.(FeatureOptionKey.N_RUNS.value) = 10;
                    end
                case FeatureName.RANDOM_NAN.value
                    if ~isfield(obj.options, FeatureOptionKey.N_RUNS.value)
                        obj.options.(FeatureOptionKey.N_RUNS.value) = 10;
                    end
                    if ~isfield(obj.options, FeatureOptionKey.RATE_NAN.value)
                        obj.options.(FeatureOptionKey.RATE_NAN.value) = 0.05;
                    end
                case FeatureName.TRUNCATED.value
                    if ~isfield(obj.options, FeatureOptionKey.PERTURBED_TRAILING_ZEROS.value)
                        obj.options.(FeatureOptionKey.PERTURBED_TRAILING_ZEROS.value) = true;
                    end
                    if ~isfield(obj.options, FeatureOptionKey.N_RUNS.value)
                        if obj.options.(FeatureOptionKey.PERTURBED_TRAILING_ZEROS.value)
                            obj.options.(FeatureOptionKey.N_RUNS.value) = 10;
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
                        obj.options.(FeatureOptionKey.UNRELAXABLE_BOUNDS.value) = false;
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
                otherwise
                    error("MATLAB:Feature:UnknownFeature", "Unknown feature: " + obj.name + ".")
            end
        end

    end

    methods (Static)
        function value = default_distribution(rand_stream, n)
            if nargin < 2
                n = 1;
            end
            value = rand_stream.randn(n, 1);
        end

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
            if ~strcmp(seed, 'shuffle')
                if isnan(seed)
                    seed = 0;
                end
                if ~isrealscalar(seed)
                    error("MATLAB:Feature:SeedNotEvenReal", "The input for seed should at least be a real number.")
                end
                seed = mod(floor(seed), 2^32);
            end
            rand_stream = RandStream('mt19937ar', 'Seed', seed);

            % Convert all elements in varargin to double
            varargin = cellfun(@double, varargin, 'UniformOutput', false);

            % Generate a new seed based on the initial rand_stream and the additional arguments
            newSeed = abs(sin(1e5 * randn(rand_stream, 1)) + sum(sin(1e5 * prod(cellfun(@(x) x, varargin))))) * 1e9;
            newSeed = mod(floor(newSeed), 2^32);

            % Create a new rand_stream with the new seed
            if isnan(newSeed)
                newSeed = 0;
            end
            rand_stream = RandStream('mt19937ar', 'Seed', floor(newSeed));
        end
    end
end