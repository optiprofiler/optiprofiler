classdef FeaturedProblem < Problem
%FEATUREDPROBLEM Problem interface and recorder for one featured trial.
%
%   FP = FeaturedProblem(P, F, MAX_EVAL, SEED) builds a fresh runtime for the
%   original Problem P and canonical Feature F. MAX_EVAL is a positive integer;
%   SEED is a nonnegative integer less than 2^32. There is no fifth
%   termination_eval argument. Reusing F does not reuse this trial's histories,
%   counters or random streams. Repetition policy belongs to benchmark, not FP.
%
%   FP is a Problem subclass with an initial point x0 and the usual objective,
%   bounds, linear constraints and nonlinear constraints. fun, cub and ceq
%   return the values observed by the solver. maxcv and the stored histories
%   use the scoring reference described below, not necessarily those observed
%   values.
%
%   .. rubric:: Execution strategy
%
%   Zero effective stages use execution_strategy='identity'; one uses
%   'legacy-single'. Both retain runtime_policy='matlab-legacy-single-v1' and
%   seed_policy='legacy-run-seed', including the established MATLAB numerical
%   kernels and payload-dependent random streams. Their constraint counters
%   count recorded queries, as described below.
%   Plain entries do not create extra stages or change these strategies.
%
%   Two or more effective stages use execution_strategy='composed-views',
%   runtime_policy='matlab-composed-views-v1' and
%   seed_policy='matlab-stage-horner32-v1'. One outer recorder owns the public
%   histories and budget; stage views own local execution state and receive
%   the immediate predecessor. A stage or custom callback may legitimately
%   query its predecessor more than once. This is not a promise of one original
%   callback per public query, statistical independence, or matching random
%   samples across MATLAB and Python. Composite derivative methods explicitly
%   raise UnsupportedCompositeDerivative; the legacy derivative behavior is
%   retained for identity/single execution.
%
%   .. rubric:: Observations, reference values and histories
%
%   Observation-only value changes such as noise do not change the reference
%   used for scoring. Spatial transformations transport both observation and
%   reference coordinates. A quantized stage with ground_truth=true also
%   snaps its local reference objective and nonlinear constraints; with false
%   it snaps observations only. Its bounds and linear constraints are checked
%   at the unsnapped point in that stage's coordinates. Reference values are
%   therefore not universally the original untransformed callback values.
%
%   Reference reads for initialization, history, maxcv and final scoring do
%   not consume the public fun/cub/ceq call budget or advance served-query
%   counters. They may still evaluate user callbacks. Observation-side custom
%   probes and constraint-gate probes are served predecessor queries and may
%   consume that predecessor's local samples without becoming extra outer
%   recorder entries.
%
%   Public read-only properties include:
%
%   - problem, feature: the original P and reusable F.
%   - max_eval, seed: this trial's budget and seed.
%   - execution_strategy, runtime_policy, seed_policy: the policies above.
%   - fun_hist: reference objective values recorded at admitted fresh objective
%     evaluations.
%   - cub_hist, ceq_hist: reference nonlinear-constraint histories, with one
%     column per recorded fresh evaluation (m-by-k; 0-by-k without nonlinear
%     constraints). Row and empty callback results are stored as one column.
%     The optional record_hist=false argument suppresses storage, not the
%     public call counter or evaluation.
%   - maxcv_hist: reference maximum violations at recorded objective points,
%     not a separate sequence of constraint calls; unavailable evaluations are
%     recorded as NaN.
%   - fun_init, maxcv_init: reference values at the featured initial point.
%
%   n_eval_fun is the length of fun_hist. On every execution strategy,
%   n_eval_cub and n_eval_ceq count recorded constraint queries, i.e. history
%   columns, independent of the number of constraint components and of the
%   orientation of the callback result. Earlier identity/legacy-single
%   versions returned length(cub_hist) and length(ceq_hist), which reported
%   the component count m while fewer than m queries were recorded, cached
%   stale values once m reached MAX_EVAL, counted m*k for row results and 0
%   without constraints. This is an intentional correctness change of the
%   inherited behavior. These accessors are not general counts of user
%   callbacks or solver calls.
%
%   .. rubric:: Cache and hard stop
%
%   Each channel checks its history-based n_eval accessor against MAX_EVAL.
%   Once that accessor is at least MAX_EVAL, a call returns the last observed
%   value cached for that channel, irrespective of the new point, and appends
%   no history. On the identity/legacy-single strategies, the same per-query
%   counter is the served index passed to stochastic constraint modifiers
%   (noise, random NaN, perturbed trailing digits, custom callbacks), so
%   repeated queries at one point receive distinct indices.
%
%   A separate public-call counter is maintained for each of fun, cub and ceq.
%   Calls served from cache still count. Before accepting a call, the channel
%   raises its ExceedTerminationEval error if its counter is already at least
%   2*MAX_EVAL. Thus, after 2*MAX_EVAL accepted public call attempts on one
%   channel, its next attempt raises before evaluating or returning the cache.
%   An evaluation that later errors can already have advanced this counter;
%   a history length is not evidence that the same number of calls succeeded.
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
        execution_strategy
        runtime_policy
        seed_policy
    end

    properties (GetAccess = private, SetAccess = private)

        real_n_eval_fun
        real_n_eval_cub
        real_n_eval_ceq
        last_fun
        last_cub
        last_ceq
        kernel
        final_view
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

            stages = feature.stages;
            kernel = [];
            view = [];
            if numel(stages) > 1
                view = optiprofiler_internal.FeatureProblemView(problem);
                for i_stage = 1:numel(stages)
                    view = optiprofiler_internal.FeatureProblemView(view, stages{i_stage}, seed);
                end
                pb_struct = struct('name',view.name,'x0',view.x0,'xl',view.xl,'xu',view.xu, ...
                    'aub',view.aub,'bub',view.bub,'aeq',view.aeq,'beq',view.beq,'fun',@(x) NaN);
            else
                if isempty(stages)
                    kernel = optiprofiler_internal.FeatureKernel('plain',struct());
                else
                    kernel = optiprofiler_internal.FeatureKernel(stages{1}.name,stages{1}.options);
                end
                pb_struct = struct();
                pb_struct.name = problem.name;
                % Modify the initial point.
                pb_struct.x0 = kernel.modifier_x0(seed, problem);
                % Modify the bounds.
                [pb_struct.xl, pb_struct.xu] = kernel.modifier_bounds(seed, problem);
                % Modify the linear inequality constraints.
                [pb_struct.aub, pb_struct.bub] = kernel.modifier_linear_ub(seed, problem);
                % Modify the linear equality constraints.
                [pb_struct.aeq, pb_struct.beq] = kernel.modifier_linear_eq(seed, problem);
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

            end

            % Initialize the FeaturedProblem object.
            obj@Problem(pb_struct);
            obj.kernel = kernel;
            obj.final_view = view;
            if isempty(view)
                if isempty(stages), obj.execution_strategy = 'identity';
                else, obj.execution_strategy = 'legacy-single'; end
                obj.runtime_policy = 'matlab-legacy-single-v1';
                obj.seed_policy = 'legacy-run-seed';
            else
                obj.execution_strategy = 'composed-views';
                obj.runtime_policy = 'matlab-composed-views-v1';
                obj.seed_policy = 'matlab-stage-horner32-v1';
            end
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

            if ~isempty(obj.final_view)
                [obj.fun_init,obj.maxcv_init] = obj.evaluateTruth(obj.x0);
                return
            end

            % Evaluate the objective function and the maximum constraint violation at the initial point.
            % Pay attention to the case when the feature is 'quantized' and the option ``ground_truth'' is set to true.
            % Note: For some problems (e.g., 'NOZZLEfp' from S2MPJ), the evaluation of the objective function or
            % constraints at the initial point may fail and return an empty value, especially when a feature
            % (e.g., 'linearly_transformed') is applied. This is likely due to the numerical sensitivity of the
            % underlying problem interfaces (e.g., MEX files) to the tiny numerical perturbations (around machine
            % precision) introduced by the affine transformation. A specific example of this was observed during
            % a random test with seed 2632 on problem 'NOZZLEfp'. To handle this, we check whether the evaluation
            % returns an empty value and, if so, attempt to use the evaluation at the original initial point.
            [A, b] = obj.kernel.modifier_affine(obj.seed, obj.problem);
            if strcmp(obj.kernel.name, FeatureName.QUANTIZED.value) && obj.kernel.options.(FeatureOptionKey.GROUND_TRUTH.value)
                val = obj.kernel.modifier_fun(A * obj.x0 + b, obj.seed, obj.problem, obj.n_eval_fun);
                if isempty(val)
                    val = obj.kernel.modifier_fun(obj.problem.x0, obj.seed, obj.problem, obj.n_eval_fun);
                    if isempty(val)
                        val = NaN;
                    end
                end
                obj.fun_init = val;
            else
                val = obj.problem.fun(A * obj.x0 + b);
                % We check whether `val` is empty.
                if isempty(val)
                    val = obj.problem.fun(obj.problem.x0); 
                    if isempty(val)
                        val = NaN;
                    end
                end
                obj.fun_init = val;
            end

            % Similar check for maxcv_init.
            % The initial violation must use the same truth as histories;
            % maxcv does not consume any solver constraint evaluations.
            if strcmp(obj.kernel.name, FeatureName.QUANTIZED.value) && obj.kernel.options.(FeatureOptionKey.GROUND_TRUTH.value)
                val = obj.maxcv(obj.x0);
            else
                val = obj.problem.maxcv(A * obj.x0 + b);
            end
            if isempty(val)
                val = obj.problem.maxcv(obj.problem.x0);
                if isempty(val)
                    val = NaN;
                end
            end
            obj.maxcv_init = val;

        end

        function value = get.n_eval_fun(obj)
            % Return number of objective function evaluations.

            value = length(obj.fun_hist);
        end

        function value = get.n_eval_cub(obj)
            % Number of recorded nonlinear inequality constraint queries: one
            % history column per query on every strategy (not length(), which
            % returned the component count for m-by-k histories).

            value = size(obj.cub_hist, 2);
        end

        function value = get.n_eval_ceq(obj)
            % Number of recorded nonlinear equality constraint queries: one
            % history column per query on every strategy.

            value = size(obj.ceq_hist, 2);
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

            if ~isempty(obj.final_view)
                f = obj.final_view.fun(x);
                obj.last_fun = f;
                obj.fun_hist = [obj.fun_hist,obj.final_view.reference('fun',x)];
                try
                    obj.maxcv_hist = [obj.maxcv_hist,obj.final_view.referenceMaxcv(x)];
                catch
                    obj.maxcv_hist = [obj.maxcv_hist,NaN];
                end
                return
            end

            % Generate the affine transformation.
            [A, b] = obj.kernel.modifier_affine(obj.seed, obj.problem);

            % Evaluate the modified the objective function value according to the feature and return the
            % modified value. We should not store the modified value because the performance 
            % of an optimization solver should be measured using the original objective function.
            f = obj.kernel.modifier_fun(A * x + b, obj.seed, obj.problem, obj.n_eval_fun);
            obj.last_fun = f;

            % Evaluate the objective function and store the results.
            f_true = obj.problem.fun(A * x + b);

            % If the feature is 'quantized' and the option ``ground_truth'' is set to true, we should
            % set f_true to f.
            if strcmp(obj.kernel.name, FeatureName.QUANTIZED.value) && obj.kernel.options.(FeatureOptionKey.GROUND_TRUTH.value)
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

            if ~isempty(obj.final_view)
                cub_ = obj.final_view.cub(x);
                obj.last_cub = cub_;
                reference_value = obj.final_view.reference('cub',x);
                if nargin < 3 || record_hist
                    obj.cub_hist = [obj.cub_hist,reference_value];
                end
                return
            end

            % Generate the affine transformation.
            [A, b] = obj.kernel.modifier_affine(obj.seed, obj.problem);

            % Evaluate the nonlinear inequality constraints and store the results.
            cub_ = obj.kernel.modifier_cub(A * x + b, obj.seed, obj.problem, obj.n_eval_cub);
            obj.last_cub = cub_;
            
            % Evaluate the nonlinear inequality constraints and store the results.
            cub_true = obj.problem.cub(A * x + b);

            % If the feature is 'quantized' and the option ``ground_truth'' is set to true, we should
            % use the modified constraint violation.
            if strcmp(obj.kernel.name, FeatureName.QUANTIZED.value) && obj.kernel.options.(FeatureOptionKey.GROUND_TRUTH.value)
                cub_true = cub_;
            end

            % Record the history of the nonlinear inequality constraints only when `record_hist` is true.
            % One column per recorded query: row and empty results become one column.
            if nargin < 3 || record_hist
                obj.cub_hist = [obj.cub_hist, cub_true(:)];
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

            if ~isempty(obj.final_view)
                ceq_ = obj.final_view.ceq(x);
                obj.last_ceq = ceq_;
                reference_value = obj.final_view.reference('ceq',x);
                if nargin < 3 || record_hist
                    obj.ceq_hist = [obj.ceq_hist,reference_value];
                end
                return
            end

            % Generate the affine transformation.
            [A, b] = obj.kernel.modifier_affine(obj.seed, obj.problem);

            % Evaluate the nonlinear equality constraints and store the results.
            ceq_ = obj.kernel.modifier_ceq(A * x + b, obj.seed, obj.problem, obj.n_eval_ceq);
            obj.last_ceq = ceq_;

            % Evaluate the nonlinear equality constraints and store the results.
            ceq_true = obj.problem.ceq(A * x + b);

            % If the Feature is ``quantized'' and the option ``ground_truth'' is set to true, we should
            % use the modified constraint violation.
            if strcmp(obj.kernel.name, FeatureName.QUANTIZED.value) && obj.kernel.options.(FeatureOptionKey.GROUND_TRUTH.value)
                ceq_true = ceq_;
            end

            % Record the history of the nonlinear equality constraints only when `record_hist` is true.
            % One column per recorded query: row and empty results become one column.
            if nargin < 3 || record_hist
                obj.ceq_hist = [obj.ceq_hist, ceq_true(:)];
            end
        end

        function varargout = maxcv(obj, x, detailed)
            if nargin < 3, detailed = false; end
            if ~isempty(obj.final_view)
                [varargout{1:nargout}] = obj.final_view.referenceMaxcv(x,detailed);
                return
            end
            % Preserve the scalar legacy path; detailed reference reads are new.
            if detailed
                [A,b] = obj.kernel.modifier_affine(obj.seed,obj.problem);
                if ~strcmp(obj.kernel.name,'quantized') || ~obj.kernel.options.ground_truth
                    [varargout{1:nargout}] = obj.problem.maxcv(A*x+b,true);
                else
                    error('MATLAB:FeaturedProblem:LegacyDetailedViolationUnsupported', ...
                        'Detailed quantized violation is not part of the legacy recorder interface.');
                end
                return
            end
            cv = obj.legacyMaxcv(x);
            varargout{1} = cv;
        end

        function cv = legacyMaxcv(obj, x)
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
            if strcmp(obj.kernel.name, FeatureName.QUANTIZED.value) && obj.kernel.options.(FeatureOptionKey.GROUND_TRUTH.value)
                if strcmp(obj.ptype, 'u')
                    cv = 0;
                    return
                end
                
                if any(isfinite(obj.xl))
                    cv_bounds = max([obj.xl - x; 0], [], 'includenan');
                else
                    cv_bounds = 0;
                end
                if any(isfinite(obj.xu))
                    cv_bounds = max([x - obj.xu; cv_bounds], [], 'includenan');
                end
                if strcmp(obj.ptype, 'b')
                    cv = cv_bounds;
                    return
                end

                if ~isempty(obj.aub)
                    cv_linear = max([obj.aub * x - obj.bub; 0], [], 'includenan');
                else
                    cv_linear = 0;
                end
                if ~isempty(obj.aeq)
                    cv_linear = max([abs(obj.aeq * x - obj.beq); cv_linear], [], 'includenan');
                end
                if strcmp(obj.ptype, 'l')
                    cv = max([cv_bounds; cv_linear], [], 'includenan');
                    return
                end

                if ~isempty(obj.cub_)
                    % Do not call the public oracle even with record_hist=false:
                    % it consumes the real budget and may return a cached value.
                    cub_val = obj.kernel.modifier_cub(x, obj.seed, obj.problem, obj.n_eval_cub);
                    if ~isempty(cub_val)
                        cv_nonlinear = max([cub_val(:); 0], [], 'includenan');
                    else
                        cv_nonlinear = 0;
                    end
                else
                    cv_nonlinear = 0;
                end
                if ~isempty(obj.ceq_)
                    ceq_val = obj.kernel.modifier_ceq(x, obj.seed, obj.problem, obj.n_eval_ceq);
                    if ~isempty(ceq_val)
                        cv_nonlinear = max([abs(ceq_val(:)); cv_nonlinear], [], 'includenan');
                    end
                end

                cv = max([cv_bounds; cv_linear; cv_nonlinear], [], 'includenan');
            else
                % Generate the affine transformation.
                [A, b] = obj.kernel.modifier_affine(obj.seed, obj.problem);
                cv = obj.problem.maxcv(A * x + b);
            end
        end


        function value=grad(obj,x)
            if isempty(obj.final_view), value=grad@Problem(obj,x); else, value=obj.final_view.grad(x); end
        end
        function value=hess(obj,x)
            if isempty(obj.final_view), value=hess@Problem(obj,x); else, value=obj.final_view.hess(x); end
        end
        function value=jcub(obj,x)
            if isempty(obj.final_view), value=jcub@Problem(obj,x); else, value=obj.final_view.jcub(x); end
        end
        function value=jceq(obj,x)
            if isempty(obj.final_view), value=jceq@Problem(obj,x); else, value=obj.final_view.jceq(x); end
        end
        function value=hcub(obj,x)
            if isempty(obj.final_view), value=hcub@Problem(obj,x); else, value=obj.final_view.hcub(x); end
        end
        function value=hceq(obj,x)
            if isempty(obj.final_view), value=hceq@Problem(obj,x); else, value=obj.final_view.hceq(x); end
        end

    end

    methods (Access = protected)
        function value = constraintDimension(obj,channel)
            if isempty(obj.final_view)
                value = constraintDimension@Problem(obj,channel);
            elseif strcmp(channel,'cub')
                value = obj.final_view.m_nonlinear_ub;
            else
                value = obj.final_view.m_nonlinear_eq;
            end
        end
        function x = constraintProbePoint(obj)
            % Inherited cub_/ceq_ are base callbacks, while obj.x0 is in solver
            % coordinates. Probe dimensions at the original problem's point,
            % or an affine map can falsely make a valid callback unavailable.
            x = obj.problem.x0;
        end
    end

    methods (Hidden)
        function x = toOriginalCoordinates(obj,x)
            if isempty(obj.final_view)
                [A,b] = obj.kernel.modifier_affine(obj.seed,obj.problem);
                x = A*x+b;
            else
                x = obj.final_view.toOriginalCoordinates(x);
            end
        end
        function receipt = runtimeReceipt(obj)
            stages = {};
            if ~isempty(obj.final_view), stages = obj.final_view.runtimeStages(); end
            receipt = struct('language','matlab','execution_strategy',obj.execution_strategy, ...
                'runtime_policy',obj.runtime_policy,'seed_policy',obj.seed_policy, ...
                'run_seed',obj.seed,'stages',{stages});
        end
        function [f, cv] = evaluateTruth(obj, x)
            if ~isempty(obj.final_view)
                f = obj.final_view.reference('fun',x);
                cv = obj.final_view.referenceMaxcv(x);
                return
            end
            % Internal scoring at solver coordinates, without an oracle call.
            % True uses the quantized problem itself; False and other features
            % retain base truth. Never snap the returned point, append history,
            % consume budget, or overwrite the last oracle value here.
            [A, b] = obj.kernel.modifier_affine(obj.seed, obj.problem);
            if strcmp(obj.kernel.name, FeatureName.QUANTIZED.value) && obj.kernel.options.(FeatureOptionKey.GROUND_TRUTH.value)
                f = obj.kernel.modifier_fun(A * x + b, obj.seed, obj.problem, obj.n_eval_fun);
                cv = obj.maxcv(x);
            else
                f = obj.problem.fun(A * x + b);
                cv = obj.problem.maxcv(A * x + b);
            end
        end
    end
end
