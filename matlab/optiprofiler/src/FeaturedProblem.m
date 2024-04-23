classdef FeaturedProblem < Problem
%FEATUREDPROBLEM modifies the original optimization problem by applying a feature.

    properties (GetAccess = public, SetAccess = private)

        feature
        max_eval
        seed
        permutation
        rotation
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
                error("The argument feature must be an instance of the class Feature.");
            end

            % Preprocess the maximum number of function evaluations.
            if ~isreal(max_eval)
                error("The argument feature max_eval must be a real number.");
            end
            if max_eval < 1 || max_eval ~= floor(max_eval)
                error("The argument max_eval must be a positive integer.");
            end

            % Preprocess the seed.
            if ~isreal(seed)
                error("The argument feature seed must be a real number.");
            end
            if seed < 0 || seed ~= floor(seed)
                error("The argument seed must be a nonnegative integer.");
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
            if strcmp(feature.name, FeatureName.ROTATED.value)
                [Q, R] = qr(rand_stream.randn(problem.n));
                rotation = Q * diag(sign(diag(R)));
            end
            
            % Copy the problem.
            if ~isa(problem, 'Problem')
                error("The argument problem must be an instance of the class Problem.");
            end
            pb_struct = struct('fun', problem.fun_, 'x0', problem.x0, 'xl', problem.xl, ...
                'xu', problem.xu, 'aub', problem.aub, 'bub', problem.bub, 'aeq', problem.aeq, ...
                'beq', problem.beq, 'cub', problem.cub_, 'ceq', problem.ceq_, 'm_nonlinear_ub', problem.m_nonlinear_ub_, ...
                'm_nonlinear_eq', problem.m_nonlinear_eq_);

            % Randomize the initial point if feature is 'perturbed_x0', and permute the initial point if feature is 'permuted'.
            if strcmp(feature.name, FeatureName.PERTURBED_X0.value)
                x0_Cell = num2cell(pb_struct.x0);
                rand_stream = feature.default_rng(seed, x0_Cell{:});
                pb_struct.x0 = pb_struct.x0 + feature.options.(FeatureOptionKey.DISTRIBUTION.value)(rand_stream, length(pb_struct.x0));
            elseif strcmp(feature.name, FeatureName.PERMUTED.value)
                pb_struct.x0 = pb_struct.x0(reverse_permutation);
            elseif strcmp(feature.name, FeatureName.ROTATED.value)
                pb_struct.x0 = rotation' * pb_struct.x0;
            end

            % Permute xl, xu, aub, aeq if feature is 'permuted'.
            if strcmp(feature.name, FeatureName.PERMUTED.value)
                pb_struct.xl = pb_struct.xl(reverse_permutation);
                pb_struct.xu = pb_struct.xu(reverse_permutation);
                if ~empty(pb_struct.aub)
                    pb_struct.aub = pb_struct.aub(:, reverse_permutation);
                end
                if ~empty(pb_struct.aeq)
                    pb_struct.aeq = pb_struct.aeq(:, reverse_permutation);
                end
            elseif strcmp(feature.name, FeatureName.ROTATED.value)
                if ~empty(pb_struct.aub)
                    pb_struct.aub = [rotation; -rotation; pb_struct.aub * rotation];
                    pb_struct.bub = [pb_struct.xu; -pb_struct.xl; pb_struct.bub];
                else
                    pb_struct.aub = [rotation; -rotation];
                    pb_struct.bub = [pb_struct.xu; -pb_struct.xl];
                end
                if ~empty(pb_struct.aeq)
                    pb_struct.aeq = pb_struct.aeq * rotation;
                end
                pb_struct.xl = -Inf * ones(problem.n, 1);
                pb_struct.xu = Inf * ones(problem.n, 1);
            end

            % Initialize the FeaturedProblem object.
            obj@Problem(pb_struct);
            % try
            %     obj.m_nonlinear_ub = problem.m_nonlinear_ub;
            % catch
            %     % pass
            % end
            % try
            %     obj.m_nonlinear_eq = problem.m_nonlinear_eq;
            % catch
            %     % pass
            % end

            obj.feature = feature;
            obj.max_eval = max_eval;
            obj.seed = seed;
            obj.permutation = permutation;
            obj.rotation = rotation;

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
            if strcmp(obj.feature.name, FeatureName.ROTATED.value)
                f_true = fun@Problem(obj, obj.rotation * x);
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
            elseif strcmp(obj.feature.name, FeatureName.ROTATED.value)
                f = cub@Problem(obj, obj.rotation * x);
                obj.last_cub = f;
            else
                f = cub@Problem(obj, x);
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
            elseif strcmp(obj.feature.name, FeatureName.ROTATED.value)
                f = ceq@Problem(obj, obj.rotation * x);
                obj.last_ceq = f;
            else
                f = ceq@Problem(obj, x);
                obj.last_ceq = f;
            end
        end
        
    end
end