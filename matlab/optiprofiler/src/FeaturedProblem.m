classdef FeaturedProblem < Problem
%FEATUREDPROBLEM modifies the original optimization problem by applying a feature.

    properties (GetAccess = public, SetAccess = private)

        problem
        feature
        max_eval
        seed
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
            Initialize an optimization problem with a specific feature.

            Parameters
            ----------
            problem : Problem
                The original optimization problem.
            feature : Feature
                The feature to apply to the optimization problem.
            max_eval : int
                The maximum number of function evaluations.
            seed : int
                The seed for the random number generator.
            %}

            % Preprocess the problem and the feature.
            if ~isa(problem, 'Problem')
                error("MATLAB:FeaturedProblem:NotProblemClass", "The first argument of FeaturedProblem must be an instance of the class Problem.");
            end
            if ~isa(feature, 'Feature')
                error("MATLAB:FeaturedProblem:NotFeatureClass", "The second argument of FeaturedProblem must be an instance of the class Feature.");
            end

            % Preprocess the maximum number of function evaluations.
            if ~(isintegerscalar(max_eval) && max_eval > 0)
                error("MATLAB:FeaturedProblem:max_evalNotPositiveInteger", "The argument max_eval must be a positive integer.");
            end

            % Preprocess the seed.
            if ~(isintegerscalar(seed) && seed >= 0 && seed < 2^32)
                error("MATLAB:FeaturedProblem:seedNotNonnegativeInteger", "The argument seed must be a nonnegative integer seed less than 2^32");
            end

            pb_struct = struct();
            pb_struct.name = problem.name;
            pb_struct.x_type = problem.x_type;
            % Modify the initial point.
            pb_struct.x0 = feature.modifier_x0(seed, problem);
            % Modify the bounds.
            [pb_struct.xl, pb_struct.xu] = feature.modifier_bounds(seed, problem);
            % Modify the linear inequality constraints.
            [pb_struct.aub, pb_struct.bub] = feature.modifier_linear_ub(seed, problem);
            % Modify the linear equality constraints.
            [pb_struct.aeq, pb_struct.beq] = feature.modifier_linear_eq(seed, problem);
            % First inherit some properties from the original problem.
            pb_struct.fun = problem.fun_;
            pb_struct.cub = problem.cub_;
            pb_struct.ceq = problem.ceq_;
            pb_struct.m_nonlinear_ub = problem.m_nonlinear_ub_;
            pb_struct.m_nonlinear_eq = problem.m_nonlinear_eq_;

            % Initialize the FeaturedProblem object.
            obj@Problem(pb_struct);
            obj.problem = problem;
            obj.feature = feature;
            obj.max_eval = max_eval;
            obj.seed = seed;

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
            x : double, size (n,)
                Decision variables.

            Returns
            -------
            f : double
                Modified objective function value.
            %}

            if obj.n_eval >= obj.max_eval
                % If the maximum number of function evaluations has been reached, return the value
                % of the objective function at the last point.
                f = obj.fun_hist(obj.max_eval);
                return
            end

            % Generate the affine transformation.
            [A, b] = obj.feature.modifier_affine(obj.seed, obj.problem);

            % Modified the objective function value according to the feature and return the 
            % modified value. We should not store the modified value because the performance 
            % of an optimization solver should be measured using the original objective function.
            f = obj.feature.modifier_fun(A * x + b, obj.seed, obj.problem);

            % Evaluate the objective function and store the results.
            f_true = obj.problem.fun(A * x + b);

            % If the Feature is ``quantized'' and the option ``is_truth'' is set to true, we should
            % set f_true to f.
            if strcmp(obj.feature.name, FeatureName.QUANTIZED.value) && obj.feature.options.(FeatureOptionKey.IS_TRUTH.value)
                f_true = f;
            end
            obj.fun_hist = [obj.fun_hist, f_true];
            obj.maxcv_hist = [obj.maxcv_hist, obj.maxcv(x)];
        end

        function cub_ = cub(obj, x)
            %{
            Evaluate the nonlinear inequality constraints.

            Parameters
            ----------
            x : double, size (n,)
                Decision variables.

            Returns
            -------
            cub_ : double, size (m_nonlinear_ub,)
                Modified nonlinear inequality constraints.    
            %}

            if obj.n_eval >= obj.max_eval
                % If the maximum number of function evaluations has been reached, return the value
                % of the nonlinear inequality constraints at the last point.
                cub_ = obj.last_cub;
                return
            end

            % Generate the affine transformation.
            [A, b] = obj.feature.modifier_affine(obj.seed, obj.problem);

            % Evaluate the nonlinear inequality constraints and store the results.
            cub_ = obj.feature.modifier_cub(A * x + b, obj.seed, obj.problem);
            obj.last_cub = cub_;
        end

        function ceq_ = ceq(obj, x)
            %{
            Evaluate the nonlinear equality constraints.

            Parameters
            ----------
            x : double, size (n,)
                Decision variables.

            Returns
            -------
            ceq_ : double, size (m_nonlinear_eq,)
                Modified nonlinear equality constraints.
            %}

            if obj.n_eval >= obj.max_eval
                % If the maximum number of function evaluations has been reached, return the value
                % of the nonlinear equality constraints at the last point.
                ceq_ = obj.last_ceq;
                return
            end

            % Generate the affine transformation.
            [A, b] = obj.feature.modifier_affine(obj.seed, obj.problem);

            % Evaluate the nonlinear equality constraints and store the results.
            ceq_ = obj.feature.modifier_ceq(A * x + b, obj.seed, obj.problem);
            obj.last_ceq = ceq_;
        end

        function cv = maxcv(obj, x)
            %{
            Evaluate the maximum constraint violation.

            Parameters
            ----------
            x : double, size (n,)
                Decision variables.

            Returns
            -------
            cv : double
                Maximum constraint violation.
            %}

            % Generate the affine transformation.
            [A, b] = obj.feature.modifier_affine(obj.seed, obj.problem);

            cv = obj.problem.maxcv(A * x + b);

            % If the Feature is ``quantized'' and the option ``is_truth'' is set to true, we should
            % use the modified constraint violation.
            if strcmp(obj.feature.name, FeatureName.QUANTIZED.value) && obj.feature.options.(FeatureOptionKey.IS_TRUTH.value)
                cv = maxcv@Problem(obj, x);
            end

        end
        
    end
end