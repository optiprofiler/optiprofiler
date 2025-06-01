classdef FeaturedProblem < Problem
%FEATUREDPROBLEM is a subclass of Problem class and defines an optimization problem
%   with a specific feature.
%
%   Problem and its subclass FEATUREDPROBLEM describe the following
%   optimization problem:
%
%       min fun(x)
%       s.t. xl <= x <= xu,
%            aub * x <= bub,
%            aeq * x = beq,
%            cub(x) <= 0,
%            ceq(x) = 0,
%       with initial point x0.
%
%   FEATUREDPROBLEM should be initialized by the following signature:
%
%       FP = FEATUREDPROBLEM(P, F, MAX_EVAL, SEED);
%
%   where the return FP is an instance of FEATUREDPROBLEM, the input P is an
%   instance of Problem, the input F is an instance of Feature, the input
%   MAX_EVAL is a positive integer, and the input SEED is a nonnegative
%   integer seed less than 2^32.
%
%   The output FP contains the following properties:
%
%       - problem: the original optimization problem.
%       - feature: the feature applied to the optimization problem.
%       - max_eval: the maximum number of function evaluations.
%       - seed: the seed for the random number generator.
%       - fun_hist: the history of the evaluated objective function values.
%       - cub_hist: the history of the evaluated nonlinear inequality
%         constraints.
%       - ceq_hist: the history of the evaluated nonlinear equality
%         constraints.
%       - maxcv_hist: the history of the maximum constraint violation.
%       - n_eval: the number of objective function evaluations.
%
%   The output FP contains all the methods of Problem, but the methods `fun`,
%   `cub`, `ceq`, and `maxcv` are modified by the input Feature.
%
%   Pay attention that when the number of function evaluations reaches the
%   input MAX_EVAL, the methods `fun`, `cub`, and `ceq` will return the
%   values of the objective function and constraints at the point where the
%   maximum number of function evaluations is reached, respectively.
%

    properties (GetAccess = public, SetAccess = private)

        problem
        feature
        max_eval
        seed
        fun_hist
        cub_hist
        ceq_hist
        maxcv_hist

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
                error("MATLAB:FeaturedProblem:NotProblemClass", "The first argument of `FeaturedProblem` must be an instance of the class Problem.");
            end
            if ~isa(feature, 'Feature')
                error("MATLAB:FeaturedProblem:NotFeatureClass", "The second argument of `FeaturedProblem` must be an instance of the class Feature.");
            end

            % Preprocess the maximum number of function evaluations.
            if ~(isintegerscalar(max_eval) && max_eval > 0)
                error("MATLAB:FeaturedProblem:max_evalNotPositiveInteger", "The argument `max_eval` of `FeaturedProblem` must be a positive integer.");
            end

            % Preprocess the seed.
            if ~(isintegerscalar(seed) && seed >= 0 && seed < 2^32)
                error("MATLAB:FeaturedProblem:seedNotNonnegativeInteger", "The argument `seed` of `FeaturedProblem` must be a nonnegative integer seed less than 2^32");
            end

            pb_struct = struct();
            pb_struct.name = problem.name;
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
            pb_struct.grad = problem.grad_;
            pb_struct.hess = problem.hess_;
            pb_struct.cub = problem.cub_;
            pb_struct.ceq = problem.ceq_;
            pb_struct.jcub = problem.jcub_;
            pb_struct.jceq = problem.jceq_;
            pb_struct.hcub = problem.hcub_;
            pb_struct.hceq = problem.hceq_;

            % Initialize the FeaturedProblem object.
            obj@Problem(pb_struct);
            obj.problem = problem;
            obj.feature = feature;
            obj.max_eval = max_eval;
            obj.seed = seed;

            % Initialize the history of the objective function, nonlinear inequality and equality constraints, and
            % the maximum constraint violation.
            % Note: maxcv_hist records the maximum constraint violation only at the points where the objective function
            % is evaluated.
            obj.fun_hist = [];
            obj.cub_hist = [];
            obj.ceq_hist = [];
            obj.maxcv_hist = [];

        end

        function value = get.n_eval(obj)
            % Return number of objective function evaluations.

            value = length(obj.fun_hist);
        end

        function value = get.fun_hist(obj)
            value = obj.fun_hist;
        end

        function value = get.cub_hist(obj)
            value = obj.cub_hist;
        end

        function value = get.ceq_hist(obj)
            value = obj.ceq_hist;
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

            % If the Feature is ``quantized'' and the option ``ground_truth'' is set to true, we should
            % set f_true to f.
            if strcmp(obj.feature.name, FeatureName.QUANTIZED.value) && obj.feature.options.(FeatureOptionKey.GROUND_TRUTH.value)
                f_true = f;
            end
            obj.fun_hist = [obj.fun_hist, f_true];
            try
                obj.maxcv_hist = [obj.maxcv_hist, obj.maxcv(x)];
            catch
                obj.maxcv_hist = [obj.maxcv_hist, NaN];
            end
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

            if length(obj.cub_hist) >= obj.max_eval
                % If the maximum number of function evaluations has been reached, return the value
                % of the nonlinear inequality constraints at the last point.
                cub_ = obj.cub_hist(obj.max_eval);
                return
            end

            % Generate the affine transformation.
            [A, b] = obj.feature.modifier_affine(obj.seed, obj.problem);

            % Evaluate the nonlinear inequality constraints and store the results.
            cub_ = obj.feature.modifier_cub(A * x + b, obj.seed, obj.problem);
            
            % Evaluate the nonlinear inequality constraints and store the results.
            cub_true = obj.problem.cub(A * x + b);

            % If the Feature is ``quantized'' and the option ``ground_truth'' is set to true, we should
            % use the modified constraint violation.
            if strcmp(obj.feature.name, FeatureName.QUANTIZED.value) && obj.feature.options.(FeatureOptionKey.GROUND_TRUTH.value)
                cub_true = cub_;
            end
            obj.cub_hist = [obj.cub_hist, cub_true];
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

            if length(obj.ceq_hist) >= obj.max_eval
                % If the maximum number of function evaluations has been reached, return the value
                % of the nonlinear equality constraints at the last point.
                ceq_ = obj.ceq_hist(obj.max_eval);
                return
            end

            % Generate the affine transformation.
            [A, b] = obj.feature.modifier_affine(obj.seed, obj.problem);

            % Evaluate the nonlinear equality constraints and store the results.
            ceq_ = obj.feature.modifier_ceq(A * x + b, obj.seed, obj.problem);

            % Evaluate the nonlinear equality constraints and store the results.
            ceq_true = obj.problem.ceq(A * x + b);

            % If the Feature is ``quantized'' and the option ``ground_truth'' is set to true, we should
            % use the modified constraint violation.
            if strcmp(obj.feature.name, FeatureName.QUANTIZED.value) && obj.feature.options.(FeatureOptionKey.GROUND_TRUTH.value)
                ceq_true = ceq_;
            end
            obj.ceq_hist = [obj.ceq_hist, ceq_true];
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

            % If the Feature is ``quantized'' and the option ``ground_truth'' is set to true, we should
            % use the modified constraint violation.
            if strcmp(obj.feature.name, FeatureName.QUANTIZED.value) && obj.feature.options.(FeatureOptionKey.GROUND_TRUTH.value)
                % `maxcv@Problem` is a method of the class `Problem`. By using this, maxcv will use obj.cub_ and
                % obj.ceq_ instead of obj.problem.cub_ and obj.problem.ceq_.
                cv = maxcv@Problem(obj, x);
            end

        end
        
    end
end