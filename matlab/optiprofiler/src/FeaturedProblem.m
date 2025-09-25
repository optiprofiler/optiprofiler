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
%   MAX_EVAL is a positive integer, and the input SEED is a nonnegative integer
%   seed less than 2^32.
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
%       - n_eval_fun: the minimum between the number of objective function
%         evaluations and max_eval.
%       - n_eval_cub: the minimum between the number of nonlinear inequality
%         constraint evaluations and max_eval.
%       - n_eval_ceq: the minimum between the number of nonlinear equality
%         constraint evaluations and max_eval.
%       - fun_init: the objective function value at the initial point.
%       - maxcv_init: the maximum constraint violation at the initial point.
%
%   The output FP contains all the methods of Problem, but the methods `fun`,
%   `cub`, `ceq`, and `maxcv` are modified by the input Feature.
%
%   Note the following two points.
%
%   1. When the number of function evaluations reaches the input MAX_EVAL, the
%   methods `fun`, `cub`, and `ceq` will return the values of the objective
%   function and constraints at the point where the maximum number of function
%   evaluations is reached, respectively.
%
%   2. When the number of function evaluations reaches two times the input
%   MAX_EVAL, the methods `fun`, `cub`, and `ceq` will raise an error
%   to terminate the optimization process.
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
        fun_init
        maxcv_init
    end

    properties (GetAccess = private, SetAccess = private)

        real_n_eval_fun
        real_n_eval_cub
        real_n_eval_ceq
        last_fun
        last_cub
        last_ceq
    end

    properties (Dependent)

        n_eval_fun
        n_eval_cub
        n_eval_ceq
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
            obj.real_n_eval_fun = 0;
            obj.real_n_eval_cub = 0;
            obj.real_n_eval_ceq = 0;

            % Initialize the history of the objective function, nonlinear inequality and equality constraints, and
            % the maximum constraint violation.
            % Note: maxcv_hist records the maximum constraint violation only at the points where the objective function
            % is evaluated.
            obj.fun_hist = [];
            obj.cub_hist = [];
            obj.ceq_hist = [];
            obj.maxcv_hist = [];
            obj.last_fun = NaN;
            obj.last_cub = NaN;
            obj.last_ceq = NaN;

            % Evaluate the objective function and the maximum constraint violation at the initial point.
            % Pay attention to the case when the feature is 'quantized' and the option ``ground_truth'' is set to true.
            [A, b] = obj.feature.modifier_affine(obj.seed, obj.problem);
            if strcmp(obj.feature.name, FeatureName.QUANTIZED.value) && obj.feature.options.(FeatureOptionKey.GROUND_TRUTH.value)
                obj.fun_init = obj.feature.modifier_fun(A * obj.x0 + b, obj.seed, obj.problem, obj.n_eval_fun);
            else
                obj.fun_init = obj.problem.fun(A * obj.x0 + b);
            end
            obj.maxcv_init = obj.problem.maxcv(A * obj.x0 + b);

        end

        function value = get.n_eval_fun(obj)
            % Return number of objective function evaluations.

            value = length(obj.fun_hist);
        end

        function value = get.n_eval_cub(obj)
            % Return number of nonlinear inequality constraint evaluations.

            value = length(obj.cub_hist);
        end

        function value = get.n_eval_ceq(obj)
            % Return number of nonlinear equality constraint evaluations.

            value = length(obj.ceq_hist);
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

            if obj.real_n_eval_fun >= 2 * obj.max_eval
                error("MATLAB:FeaturedProblem:funExceedTerminationEval", "The number of the objective function evaluations has reached %d (two times the maximum function evaluations).", 2 * obj.max_eval);
            end
            obj.real_n_eval_fun = obj.real_n_eval_fun + 1;

            if obj.n_eval_fun >= obj.max_eval
                % If the maximum number of function evaluations has been reached,
                % return the last evaluated objective function value.
                f = obj.last_fun;
                return
            end

            % Generate the affine transformation.
            [A, b] = obj.feature.modifier_affine(obj.seed, obj.problem);

            % Evaluate the modified the objective function value according to the feature and return the
            % modified value. We should not store the modified value because the performance 
            % of an optimization solver should be measured using the original objective function.
            f = obj.feature.modifier_fun(A * x + b, obj.seed, obj.problem, obj.n_eval_fun);
            obj.last_fun = f;

            % Evaluate the objective function and store the results.
            f_true = obj.problem.fun(A * x + b);

            % If the feature is 'quantized' and the option ``ground_truth'' is set to true, we should
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

        function cub_ = cub(obj, x, record_hist)
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

            if obj.real_n_eval_cub >= 2 * obj.max_eval
                error("MATLAB:FeaturedProblem:cubExceedTerminationEval", "The number of the nonlinear inequality constraint evaluations has reached %d (two times the maximum function evaluations).", 2 * obj.max_eval);
            end
            obj.real_n_eval_cub = obj.real_n_eval_cub + 1;

            if obj.n_eval_cub >= obj.max_eval
                % If the maximum number of function evaluations has been reached,
                % return the last evaluated nonlinear inequality constraints.
                cub_ = obj.last_cub;
                return
            end

            % Generate the affine transformation.
            [A, b] = obj.feature.modifier_affine(obj.seed, obj.problem);

            % Evaluate the nonlinear inequality constraints and store the results.
            cub_ = obj.feature.modifier_cub(A * x + b, obj.seed, obj.problem, obj.n_eval_cub);
            obj.last_cub = cub_;
            
            % Evaluate the nonlinear inequality constraints and store the results.
            cub_true = obj.problem.cub(A * x + b);

            % If the feature is 'quantized' and the option ``ground_truth'' is set to true, we should
            % use the modified constraint violation.
            if strcmp(obj.feature.name, FeatureName.QUANTIZED.value) && obj.feature.options.(FeatureOptionKey.GROUND_TRUTH.value)
                cub_true = cub_;
            end

            % Record the history of the nonlinear inequality constraints only when `record_hist` is true.
            if nargin < 3 || record_hist
                obj.cub_hist = [obj.cub_hist, cub_true];
            end
        end

        function ceq_ = ceq(obj, x, record_hist)
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

            if obj.real_n_eval_ceq >= 2 * obj.max_eval
                error("MATLAB:FeaturedProblem:ceqExceedTerminationEval", "The number of the nonlinear equality constraint evaluations has reached %d (two times the maximum function evaluations).", 2 * obj.max_eval);
            end
            obj.real_n_eval_ceq = obj.real_n_eval_ceq + 1;

            if obj.n_eval_ceq >= obj.max_eval
                % If the maximum number of function evaluations has been reached,
                % return the last evaluated nonlinear equality constraints.
                ceq_ = obj.last_ceq;
                return
            end

            % Generate the affine transformation.
            [A, b] = obj.feature.modifier_affine(obj.seed, obj.problem);

            % Evaluate the nonlinear equality constraints and store the results.
            ceq_ = obj.feature.modifier_ceq(A * x + b, obj.seed, obj.problem, obj.n_eval_ceq);
            obj.last_ceq = ceq_;

            % Evaluate the nonlinear equality constraints and store the results.
            ceq_true = obj.problem.ceq(A * x + b);

            % If the Feature is ``quantized'' and the option ``ground_truth'' is set to true, we should
            % use the modified constraint violation.
            if strcmp(obj.feature.name, FeatureName.QUANTIZED.value) && obj.feature.options.(FeatureOptionKey.GROUND_TRUTH.value)
                ceq_true = ceq_;
            end

            % Record the history of the nonlinear equality constraints only when `record_hist` is true.
            if nargin < 3 || record_hist
                obj.ceq_hist = [obj.ceq_hist, ceq_true];
            end
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

            % If the Feature is ``quantized'' and the option ``ground_truth'' is set to true, we should
            % use the modified constraint violation.
            if strcmp(obj.feature.name, FeatureName.QUANTIZED.value) && obj.feature.options.(FeatureOptionKey.GROUND_TRUTH.value)
                if strcmp(obj.ptype, 'u')
                    cv = 0;
                    return
                end
                
                if any(isfinite(obj.xl))
                    cv_bounds = max(max(obj.xl - x), 0);
                else
                    cv_bounds = 0;
                end
                if any(isfinite(obj.xu))
                    cv_bounds = max(max(x - obj.xu), cv_bounds);
                end
                if strcmp(obj.ptype, 'b')
                    cv = max(cv_bounds);
                    return
                end

                if ~isempty(obj.aub)
                    cv_linear = max(max(obj.aub * x - obj.bub), 0);
                else
                    cv_linear = 0;
                end
                if ~isempty(obj.aeq)
                    cv_linear = max(max(abs(obj.aeq * x - obj.beq)), cv_linear);
                end
                if strcmp(obj.ptype, 'l')
                    cv = max([cv_bounds; cv_linear]);
                    return
                end

                if ~isempty(obj.cub_)
                    cub_val = obj.cub(x, false);    % Do not record history
                    if ~isempty(cub_val)
                        cv_nonlinear = max(max(cub_val), 0);
                    else
                        cv_nonlinear = 0;
                    end
                else
                    cv_nonlinear = 0;
                end
                if ~isempty(obj.ceq_)
                    ceq_val = obj.ceq(x, false);    % Do not record history
                    if ~isempty(ceq_val)
                        cv_nonlinear = max(max(abs(ceq_val)), cv_nonlinear);
                    end
                end

                cv = max([cv_bounds; cv_linear; cv_nonlinear]);
            else
                % Generate the affine transformation.
                [A, b] = obj.feature.modifier_affine(obj.seed, obj.problem);
                cv = obj.problem.maxcv(A * x + b);
            end
        end


        % Note: We need to add methods `grad`, `hess`, `jcub`, and `jceq` to the FeaturedProblem class in the future.

    end
end