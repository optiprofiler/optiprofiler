classdef FeaturedProblem < Problem
%FEATUREDPROBLEM modifies the original optimization problem by applying a feature.

    properties (GetAccess = public, SetAccess = private)

        feature
        max_eval
        seed
        permutation
        rotation
        scaler
        fun_hist
        maxcv_hist
        last_cub
        last_ceq

    end

    properties (Dependent)

        n_eval

    end

    methods

        function obj = FeaturedProblem(problem, feature, max_eval, seed)
            %{
            Initialize an optimization problem.

            Parameters
            ----------
            problem : `optiprofiler.problems.Problem`
                Problem to be used in the benchmarking.
            feature : `optiprofiler.features.Feature`
                Feature to be used in the benchmarking.
            max_eval : int
                Maximum number of function evaluations.
            seed : int, optional
                Seed for the random number generator.
            %}

            % Set the default value of the seed.
            if nargin < 4
                seed = mod(floor(posixtime(datetime('now'))), 2^32);
            end

            % Preprocess the feature.
            if ~isa(feature, 'Feature')
                error("MATLAB:FeaturedProblem:NotFeatureClass", "The argument feature must be an instance of the class Feature.");
            end

            % Preprocess the maximum number of function evaluations.
            if ~(isintegerscalar(max_eval) && max_eval > 0)
                error("MATLAB:FeaturedProblem:max_evalNotPositiveInteger", "The argument max_eval must be a positive integer.");
            end

            % Preprocess the seed.
            if ~(isintegerscalar(seed) && seed >= 0 && seed < 2^32)
                error("MATLAB:FeaturedProblem:seedNotNonnegativeInteger", "The argument seed must be a nonnegative integer seed less than 2^32");
            end

            rand_stream = feature.default_rng(seed);
            % Generate a random permutation.
            permutation = [];
            if strcmp(feature.name, FeatureName.PERMUTED.value)
                permutation = rand_stream.randperm(problem.n);
                [~, reverse_permutation] = sort(permutation);
            end
            % Generate a random rotation.
            rotation = [];
            if strcmp(feature.name, FeatureName.LINEARLY_TRANSFORMED.value)
                if feature.options.(FeatureOptionKey.ROTATED.value)
                    [Q, R] = qr(rand_stream.randn(problem.n));
                    rotation = Q * diag(sign(diag(R)));
                else
                    rotation = eye(problem.n);
                end
            end
            % Generate a random scaling matrix.
            scaler = [];
            if strcmp(feature.name, FeatureName.LINEARLY_TRANSFORMED.value)
                if strcmp(feature.options.(FeatureOptionKey.CONDITION_NUMBER.value), 'dimension_dependent')
                    condition_number = min(1e8, 2^(problem.n));
                else
                    condition_number = feature.options.(FeatureOptionKey.CONDITION_NUMBER.value);
                end
                power = log2(condition_number) / 2;
                scaler = -power + 2 * power * rand_stream.rand(problem.n - 2, 1);
                scaler = 2 .^ [-power; scaler; power];
                scaler = scaler(rand_stream.randperm(problem.n));
            end
            
            % Copy the problem.
            if ~isa(problem, 'Problem')
                error("MATLAB:FeaturedProblem:NotProblemClass", "The argument problem must be an instance of the class Problem.");
            end
            pb_struct = struct('fun', problem.fun_, 'x0', problem.x0, 'xl', problem.xl, ...
                'xu', problem.xu, 'aub', problem.aub, 'bub', problem.bub, 'aeq', problem.aeq, ...
                'beq', problem.beq, 'cub', problem.cub_, 'ceq', problem.ceq_, 'm_nonlinear_ub', problem.m_nonlinear_ub_, ...
                'm_nonlinear_eq', problem.m_nonlinear_eq_);

            % Randomize the initial point if feature is 'perturbed_x0', and permute the initial point if feature is 'permuted'.
            if strcmp(feature.name, FeatureName.PERTURBED_X0.value)
                if feature.options.(FeatureOptionKey.NOISE_TYPE.value) == NoiseType.ABSOLUTE.value
                    pb_struct.x0 = pb_struct.x0 + feature.options.(FeatureOptionKey.NOISE_LEVEL.value) * feature.options.(FeatureOptionKey.DISTRIBUTION.value)(rand_stream, length(pb_struct.x0));
                else
                    % Use max(1, norm(x0)) to avoid no perturbation when x0 is zero.
                    pb_struct.x0 = pb_struct.x0 + feature.options.(FeatureOptionKey.NOISE_LEVEL.value) * max(1, norm(pb_struct.x0)) * feature.options.(FeatureOptionKey.DISTRIBUTION.value)(rand_stream, length(pb_struct.x0));
                end
            elseif strcmp(feature.name, FeatureName.PERMUTED.value)
                pb_struct.x0 = pb_struct.x0(reverse_permutation);
            elseif strcmp(feature.name, FeatureName.LINEARLY_TRANSFORMED.value)
                pb_struct.x0 = (1 ./ scaler) .* (rotation' * pb_struct.x0);
            end

            % Permute xl, xu, aub, aeq if feature is 'permuted'.
            if strcmp(feature.name, FeatureName.PERMUTED.value)
                pb_struct.xl = pb_struct.xl(reverse_permutation);
                pb_struct.xu = pb_struct.xu(reverse_permutation);
                if ~isempty(pb_struct.aub)
                    pb_struct.aub = pb_struct.aub(:, reverse_permutation);
                end
                if ~isempty(pb_struct.aeq)
                    pb_struct.aeq = pb_struct.aeq(:, reverse_permutation);
                end
            elseif strcmp(feature.name, FeatureName.LINEARLY_TRANSFORMED.value)
                if ~isempty(pb_struct.aub)
                    pb_struct.aub = [rotation * diag(scaler); -rotation * diag(scaler); pb_struct.aub * rotation * diag(scaler)];
                    pb_struct.bub = [pb_struct.xu; -pb_struct.xl; pb_struct.bub];
                else
                    pb_struct.aub = [rotation * diag(scaler); -rotation * diag(scaler)];
                    pb_struct.bub = [pb_struct.xu; -pb_struct.xl];
                end
                if ~isempty(pb_struct.aeq)
                    pb_struct.aeq = pb_struct.aeq * rotation * diag(scaler);
                end
                pb_struct.xl = -Inf * ones(problem.n, 1);
                pb_struct.xu = Inf * ones(problem.n, 1);
            end

            % Initialize the FeaturedProblem object.
            obj@Problem(pb_struct);

            obj.feature = feature;
            obj.max_eval = max_eval;
            obj.seed = seed;
            obj.permutation = permutation;
            obj.rotation = rotation;
            obj.scaler = scaler;

            % Initialize the history of the objective function and the maximum constraint violation.
            obj.fun_hist = [];
            obj.maxcv_hist = [];

            % Initialize the last evaluated nonlinear inequality and equality constraints.
            obj.last_cub = [];
            obj.last_ceq = [];
        end

        function value = get.n_eval(obj)
            % Return number of objective function evaluations.

            value = length(obj.fun_hist);
        end

        function value = get.fun_hist(obj)
            value = obj.fun_hist;
        end

        function value = get.maxcv_hist(obj)
            value = obj.maxcv_hist;
        end

        function f = fun(obj, x)
            %{
            Evaluate the objective function.

            Parameters
            ----------
            x : array_like, shape (n,)
                Point at which to evaluate the objective function.

            Returns
            -------
            float
                Value of the objective function at `x`.

            Raises
            ------
            ValueError
                If the argument `x` has an invalid shape.
            %}

            if obj.n_eval >= obj.max_eval
                % If the maximum number of function evaluations has been reached, return the value of the objective function at the last point.
                f = obj.fun_hist(obj.max_eval);
                return
            end

            % Permute the variables if necessary.
            if strcmp(obj.feature.name, FeatureName.PERMUTED.value)
                x = x(obj.permutation);
            end

            % Evaluate the objective function and store the results.
            if strcmp(obj.feature.name, FeatureName.LINEARLY_TRANSFORMED.value)
                f_true = fun@Problem(obj, obj.rotation * (obj.scaler .* x));
            else
                f_true = fun@Problem(obj, x);
            end
            obj.fun_hist = [obj.fun_hist, f_true];
            [maxcv, maxcv_bounds, maxcv_linear, maxcv_nonlinear] = obj.maxcv(x, true);
            obj.maxcv_hist = [obj.maxcv_hist, maxcv];

            % Modified the objective function value according to the feature and return the 
            % modified value. We should not store the modified value because the performance 
            % of an optimization solver should be measured using the original objective function.
            f = obj.feature.modifier(x, f_true, obj.seed, maxcv_bounds, maxcv_linear, maxcv_nonlinear);
        end

        function f = cub(obj, x)
            %{
            Evaluate the nonlinear constraints ``cub(x) <= 0``.

            Parameters
            ----------
            x : array_like, shape (n,)
                Point at which to evaluate the nonlinear inequality constraints.

            Returns
            -------
            `vector`, shape (m_nonlinear_ub,)
                Values of the nonlinear inequality constraints at `x`.

            Raises
            ------
            ValueError
                If the argument `x` has an invalid shape or if the return value of
                the argument `cub` has an invalid shape.
            %}

            % Permute the variables if necessary.
            if strcmp(obj.feature.name, FeatureName.PERMUTED.value)
                x = x(obj.permutation);
            end

            % Evaluate the nonlinear inequality constraints.
            if obj.n_eval >= obj.max_eval
                % If the maximum number of function evaluations has been reached, return the value of the nonlinear inequality constraints at the last point.
                f = obj.last_cub;
                return
            else
                if strcmp(obj.feature.name, FeatureName.LINEARLY_TRANSFORMED.value)
                    f = cub@Problem(obj, obj.rotation * (obj.scaler .* x));
                else
                    f = cub@Problem(obj, x);
                end
                obj.last_cub = f;
            end
        end

        function f = ceq(obj, x)
            %{
            Evaluate the nonlinear constraints ``ceq(x) = 0``.

            Parameters
            ----------
            x : array_like, shape (n,)
                Point at which to evaluate the nonlinear equality constraints.

            Returns
            -------
            `vector`, shape (m_nonlinear_eq,)
                Values of the nonlinear equality constraints at `x`.

            Raises
            ------
            ValueError
                If the argument `x` has an invalid shape or if the return value of
                the argument `ceq` has an invalid shape.
            %}

            % Permute the variables if necessary.
            if strcmp(obj.feature.name, FeatureName.PERMUTED.value)
                x = x(obj.permutation);
            end

            % Evaluate the nonlinear equality constraints.
            if obj.n_eval >= obj.max_eval
                % If the maximum number of function evaluations has been reached, return the value of the nonlinear equality constraints at the last point.
                f = obj.last_ceq;
                return
            else
                if strcmp(obj.feature.name, FeatureName.LINEARLY_TRANSFORMED.value)
                    f = ceq@Problem(obj, obj.rotation * (obj.scaler .* x));
                else
                    f = ceq@Problem(obj, x);
                end
                obj.last_ceq = f;
            end
        end
        
    end
end