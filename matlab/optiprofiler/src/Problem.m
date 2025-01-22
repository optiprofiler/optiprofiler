classdef Problem < handle
%PROBLEM is a class that defines an optimization problem.
%
%   PROBLEM is a class that defines an optimization problem with the following
%   structure:
%
%       min fun(x)
%       s.t. xl <= x <= xu
%            aub @ x <= bub
%            aeq @ x = beq
%            cub(x) <= 0
%            ceq(x) = 0
%       with initial point x0.

    properties (GetAccess = public, SetAccess = public)

        name = 'Unnamed Problem'

    end

    properties (GetAccess = public, SetAccess = private)

        x_type = 'real'
        x0
        xl = []
        xu = []
        aub = []
        bub = []
        aeq = []
        beq = []

    end

    properties (Dependent)

        n
        m_linear_ub
        m_linear_eq
        p_type

    end

    properties (Access = protected)

        fun_
        cub_
        ceq_
        m_nonlinear_ub_
        m_nonlinear_eq_

    end

    methods

        function obj = Problem(varargin)
            %{
            Initialize a PROBLEM.

            Parameters
            ----------
            name : str, optional
                Name of the optimization problem.
            fun : callable
                Objective function to be minimized.

                    ``fun(x) -> float``

                where ``x`` is an array with shape (n,).
            x_type : str, optional
                Type of the variables.
            x0 : array_like, shape (n,)
                Initial guess.
            xl : array_like, shape (n,), optional
                Lower bounds on the variables ``xl <= x``.
            xu : array_like, shape (n,), optional
                Upper bounds on the variables ``x <= xu``.
            aub : array_like, shape (m_linear_ub, n), optional
                Left-hand side matrix of the linear constraints ``aub @ x <= bub``.
            bub : array_like, shape (m_linear_ub,), optional
                Right-hand side vector of the linear constraints ``aub @ x <= bub``.
            aeq : array_like, shape (m_linear_eq, n), optional
                Left-hand side matrix of the linear constraints ``aeq @ x == beq``.
            beq : array_like, shape (m_linear_eq,), optional
                Right-hand side vector of the linear constraints ``aeq @ x == beq``.
            cub : callable, optional
                Nonlinear inequality constraint ``cub(x) <= 0``.

                    ``cub(x) -> array_like, shape (m_nonlinear_ub,)``

                where ``x`` is an array with shape (n,).
            ceq : callable, optional
                Nonlinear equality constraint ``ceq(x) == 0``.

                    ``ceq(x) -> array_like, shape (m_nonlinear_eq,)``

                where ``x`` is an array with shape (n,).
            m_nonlinear_ub : int, optional
                Number of nonlinear inequality constraints.
            m_nonlinear_eq : int, optional
                Number of nonlinear equality constraints.
            p_type : str, optional
                Type of the optimization problem.

            Returns
            -------
            obj : Problem
                a PROBLEM object.
            %}

            if nargin < 1
                error("MATLAB:Problem:MissingArguments", "At least one argument is provided for `Problem`.")
            end
        
            if isstruct(varargin{1})
                % Handle the case when a struct is passed
                s = varargin{1};
        
                % Check if the struct contains the required fields
                if ~isfield(s, 'fun') || ~isfield(s, 'x0')
                    error("MATLAB:Problem:MissingFields", "The fields `fun` and `x0` for `Problem` are required.")
                end
                obj.fun_ = s.fun;

                % Check if the struct contains `cub` and `ceq` fields
                if isfield(s, 'cub')
                    obj.cub_ = s.cub;
                end
                if isfield(s, 'ceq')
                    obj.ceq_ = s.ceq;
                end

                % Check if the struct contains `m_nonlinear_ub` and `m_nonlinear_eq` fields
                if isfield(s, 'm_nonlinear_ub')
                    obj.m_nonlinear_ub_ = s.m_nonlinear_ub;
                end
                if isfield(s, 'm_nonlinear_eq')
                    obj.m_nonlinear_eq_ = s.m_nonlinear_eq;
                end
        
                % Iterate over the struct's fields and assign them to the object's properties
                expected_fields = {'name', 'fun', 'x_type', 'x0', 'xl', 'xu', 'aub', 'bub', 'aeq', 'beq', 'cub', 'ceq', 'm_nonlinear_ub', 'm_nonlinear_eq'};
                fields = fieldnames(s);
                for i = 1:numel(expected_fields)
                    if strcmp(expected_fields{i}, 'fun') || strcmp(expected_fields{i}, 'cub') || strcmp(expected_fields{i}, 'ceq') || strcmp(expected_fields{i}, 'm_nonlinear_ub') || strcmp(expected_fields{i}, 'm_nonlinear_eq')
                        continue
                    elseif ismember(expected_fields{i}, fields)
                        obj.(expected_fields{i}) = s.(expected_fields{i});
                    end
                end
            else
                error("MATLAB:Problem:NotStruct", "Invalid input for `Problem`. A struct argument is expected.")
            end

            % Check that the arguments are legal and consistent.

            % Check that `x_type` belongs to {'real', 'integer', 'binary'}.
            if ~ismember(obj.x_type, {'real', 'integer', 'binary'})
                error("MATLAB:Problem:x_type_NotValid", "The argument `x_type` for `Problem` must be one of {'real', 'integer', 'binary'}.")
            end

            % Check that `x0` is a real vector.
            if ~isrealvector(obj.x0)
                error("MATLAB:Problem:x0_NotRealVector", "The argument `x0` for `Problem` must be a real vector.")
            end

            % Check that `xl` is a real vector having the same size as `x0`.
            if ~isrealvector(obj.xl) || numel(obj.xl) ~= obj.n
                error("MATLAB:Problem:xl_x0_NotConsistent", "The argument `xl` for `Problem` must be a real vector having size %d.", obj.n);
            end
            % Check that `xu` a real vector having the same size as `x0`.
            if ~isrealvector(obj.xu) || numel(obj.xu) ~= obj.n
                error("MATLAB:Problem:xu_x0_NotConsistent", "The argument `xu` for `Problem` must have size %d.", obj.n);
            end
            % Check that `bub` is a real vector.
            if ~isrealvector(obj.bub)
                error("MATLAB:Problem:bub_NotRealVector", "The argument `bub` for `Problem` must be a real vector.")
            end
            % Check that `aub` is a real matrix with shape (numel(obj.bub), n).
            if ~isrealmatrix(obj.aub) || ~isequal(size(obj.aub), [numel(obj.bub), obj.n])
                error("MATLAB:Problem:aub_m_linear_ub_n_NotConsistent", "The argument `aub` for `Problem` must be a real matrix having shape (%d, %d).", obj.m_linear_ub, obj.n);
            end
            % Check that `beq` is a real vector.
            if ~isrealvector(obj.beq)
                error("MATLAB:Problem:beq_NotRealVector", "The argument `beq` for `Problem` must be a real vector.")
            end
            % Check that `aeq` is a real matrix with shape (numel(obj.beq), n).
            if ~isrealmatrix(obj.aeq) || ~isequal(size(obj.aeq), [obj.m_linear_eq, obj.n])
                error("MATLAB:Problem:aeq_m_linear_eq_n_NotConsistent", "The argument `aeq` for `Problem` must have shape (%d, %d).", obj.m_linear_eq, obj.n);
            end
            % Check whether `m_nonlinear_ub` is empty or zero if `cub` is empty.
            if isempty(obj.cub_) && ~isempty(obj.m_nonlinear_ub_) && obj.m_nonlinear_ub_ > 0
                error("MATLAB:Problem:m_nonlinear_ub_cub_NotConsistent", "The argument `m_nonlinear_ub` for `Problem` must be empty or zero if the argument `cub` is empty.");
            end
            % Check whether `m_nonlinear_eq` is empty or zero if `ceq` is empty.
            if isempty(obj.ceq_) && ~isempty(obj.m_nonlinear_eq_) && obj.m_nonlinear_eq_ > 0
                error("MATLAB:Problem:m_nonlinear_eq_ceq_NotConsistent", "The argument `m_nonlinear_eq` for `Problem` must be empty or zero if the argument `ceq` is empty.");
            end
        end


        % Setter functions.

        % Preprocess the type of the variables.
        function set.x_type(obj, value)
            if ~ismember(value, {'real', 'integer', 'binary'})
                error("MATLAB:Problem:x_type_NotValid", "The argument `x_type` for `Problem` must be one of {'real', 'integer', 'binary'}.")
            end
            obj.x_type = value;
        end

        % Preprocess the initial guess.
        function set.x0(obj, value)
            if ~isvector(value)
                error("MATLAB:Problem:x0_NotVector", "The argument `x0` for `Problem` must be a vector.")
            end
            obj.x0 = reshape(value, [], 1);
        end

        % Preprocess the objective function.
        function set.fun_(obj, fun)
            % This function only accepts a new function handle.
            if ~isa(fun, 'function_handle')
                error("MATLAB:Problem:fun_NotFunctionHandle", "The argument `fun` for `Problem` must be a function handle.")
            end
            % Check if `fun` can accept only one argument.
            if nargin(fun) ~= 1
                error("MATLAB:Problem:fun_NotOneArguementFun", "The argument `fun` for `Problem` must accept exactly one argument.")
            end
            
            obj.fun_ = fun;
        end

        % Preprocess the bound constraints.
        function set.xl(obj, xl)
            if ~isempty(xl)
                if ~isvector(xl)
                    error("MATLAB:Problem:xl_NotVector", "The argument `xl` for `Problem` must be a vector.")
                end
                obj.xl = reshape(xl, [], 1);
            else
                obj.xl = [];
            end
        end

        function set.xu(obj, xu)
            if ~isempty(xu)
                if ~isvector(xu)
                    error("MATLAB:Problem:xu_NotVector", "The argument `xu` for `Problem` must be a vector.")
                end
                obj.xu = reshape(xu, [], 1);
            else
                obj.xu = [];
            end
        end

        % Preprocess the linear constraints.
        function set.aub(obj, aub)
            if ~isempty(aub)
                if ~ismatrix(aub)
                    error("MATLAB:Problem:aub_NotMatrix", "The argument `aub` for `Problem` must be a matrix.")
                end
                obj.aub = aub;
            else
                obj.aub = [];
            end
        end

        function set.bub(obj, bub)
            if ~isempty(bub)
                if ~isvector(bub)
                    error("MATLAB:Problem:bub_NotVector", "The argument `bub` for `Problem` must be a vector.")
                end
                obj.bub = reshape(bub, [], 1);
            else
                obj.bub = [];
            end
        end

        function set.aeq(obj, aeq)
            if ~isempty(aeq)
                if ~ismatrix(aeq)
                    error("MATLAB:Problem:aeq_NotMatrix", "The argument `aeq` for `Problem` must be a matrix.")
                end
                obj.aeq = aeq;
            else
                obj.aeq = [];
            end
        end

        function set.beq(obj, beq)
            if ~isempty(beq)
                if ~isvector(beq)
                    error("MATLAB:Problem:beq_NotVector", "The argument `beq` for `Problem` must be a vector.")
                end
                obj.beq = reshape(beq, [], 1);
            else
                obj.beq = [];
            end
        end

        % Preprocess the nonlinear constraints.
        function set.cub_(obj, cub_)
            if ~isa(cub_, 'function_handle') && ~isempty(cub_)
                error("MATLAB:Problem:cub_NotFunctionHandle", "The argument `cub` for `Problem` must be a function handle.")
            end
            obj.cub_ = cub_;
        end

        function set.ceq_(obj, ceq_)
            if ~isa(ceq_, 'function_handle') && ~isempty(ceq_)
                error("MATLAB:Problem:ceq_NotFunctionHandle", "The argument `ceq` for `Problem` must be a function handle.")
            end
            
            obj.ceq_ = ceq_;
        end

        % Preprocess the number of nonlinear constraints.
        function set.m_nonlinear_ub_(obj, m_nonlinear_ub)
            if ~isempty(m_nonlinear_ub)
                if ~isintegerscalar(m_nonlinear_ub) || m_nonlinear_ub < 0
                    error("MATLAB:Problem:m_nonlinear_ub_NotPositiveScalar", "The argument `m_nonlinear_ub` for `Problem` must be a nonnegative integer.")
                end
                obj.m_nonlinear_ub_ = m_nonlinear_ub;
            else
                obj.m_nonlinear_ub_ = [];
            end
        end

        function set.m_nonlinear_eq_(obj, m_nonlinear_eq)
            if ~isempty(m_nonlinear_eq)
                if ~isintegerscalar(m_nonlinear_eq) || m_nonlinear_eq < 0
                    error("MATLAB:Problem:m_nonlinear_eq_NotPositiveScalar", "The argument `m_nonlinear_eq` for `Problem` must be a nonnegative integer.")
                end
                obj.m_nonlinear_eq_ = m_nonlinear_eq;
            else
                obj.m_nonlinear_eq_ = [];
            end
        end

        % Getter functions for dependent properties.

        function value = get.n(obj)
            value = numel(obj.x0);
        end

        function value = get.m_linear_ub(obj)
            value = sum(~isinf(obj.bub));
        end

        function value = get.m_linear_eq(obj)
            value = numel(obj.beq);
        end

        function value = get.p_type(obj)
            try
                if obj.m_nonlinear_ub + obj.m_nonlinear_eq > 0
                    value = 'nonlinearly constrained';
                elseif obj.m_linear_ub + obj.m_linear_eq > 0
                    value = 'linearly constrained';
                elseif any(obj.xl > -inf) || any(obj.xu < inf)
                    value = 'bound-constrained';
                else
                    value = 'unconstrained';
                end
            catch ME
                value = 'nonlinearly constrained';
            end
        end

        % Other getter functions.

        function value = get.xl(obj)
            if isempty(obj.xl)
                value = -inf(obj.n, 1);
            else
                value = obj.xl;
            end
        end
    
        function value = get.xu(obj)
            if isempty(obj.xu)
                value = inf(obj.n, 1);
            else
                value = obj.xu;
            end
        end
    
        function value = get.aub(obj)
            if isempty(obj.aub)
                value = NaN(0, obj.n);
            else
                value = obj.aub;
            end
        end
    
        function value = get.bub(obj)
            if isempty(obj.bub)
                value = NaN(0, 1);
            else
                value = obj.bub;
            end
        end
    
        function value = get.aeq(obj)
            if isempty(obj.aeq)
                value = NaN(0, obj.n);
            else
                value = obj.aeq;
            end
        end
    
        function value = get.beq(obj)
            if isempty(obj.beq)
                value = NaN(0, 1);
            else
                value = obj.beq;
            end
        end

        % Other functions.

        function varargout = maxcv(obj, x, detailed)
            % Evaluate the maximum constraint violation.

            if nargin < 3
                detailed = false;
            end

            % Check if x is a vector.
            if ~isrealvector(x)
                error("MATLAB:Problem:InvalidInputForMaxCV", "The input `x` for method `maxcv` in `Problem` must be a real vector.")
            end

            % Check if x has the same size as the dimension of the problem.
            if numel(x) ~= obj.n
                error("MATLAB:Problem:WrongSizeInputForMaxCV", "The input `x` for method `maxcv` in `Problem` must have size %d.", obj.n)
            end

            if strcmp(obj.p_type, 'unconstrained')
                cv = 0;
                if detailed
                    varargout{1} = cv;
                    varargout{2} = 0;
                    varargout{3} = 0;
                    varargout{4} = 0;
                else
                    varargout{1} = cv;
                end
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
            if strcmp(obj.p_type, 'bound-constrained')
                cv = max(cv_bounds);
                if detailed
                    varargout{1} = cv;
                    varargout{2} = cv_bounds;
                    varargout{3} = 0;
                    varargout{4} = 0;
                else
                    varargout{1} = cv;
                end
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
            if strcmp(obj.p_type, 'linearly constrained')
                cv = max([cv_bounds; cv_linear]);
                if detailed
                    varargout{1} = cv;
                    varargout{2} = cv_bounds;
                    varargout{3} = cv_linear;
                    varargout{4} = 0;
                else
                    varargout{1} = cv;
                end
                return
            end

            if ~isempty(obj.cub_)
                cub_val = obj.cub(x);
                if ~isempty(cub_val)
                    cv_nonlinear = max(max(cub_val), 0);
                else
                    cv_nonlinear = 0;
                end
            else
                cv_nonlinear = 0;
            end
            if ~isempty(obj.ceq_)
                ceq_val = obj.ceq(x);
                if ~isempty(ceq_val)
                    cv_nonlinear = max(max(abs(ceq_val)), cv_nonlinear);
                end
            end

            cv = max([cv_bounds; cv_linear; cv_nonlinear]);

            if detailed
                varargout{1} = cv;
                varargout{2} = cv_bounds;
                varargout{3} = cv_linear;
                varargout{4} = cv_nonlinear;
            else
                varargout{1} = cv;
            end

        end

        function f = fun(obj, x)
            % Check if x is a vector.
            if ~isrealvector(x)
                error("MATLAB:Problem:InvalidInputForFUN", "The input `x` for method `fun` in `Problem` must be a vector.")
            end

            if numel(x) ~= obj.n
                error("MATLAB:Problem:WrongSizeInputForFUN", "The input `x` for method `fun` in `Problem` must have size %d.", obj.n)
            end

            try
                % Try to evaluate the function at x.
                f = obj.fun_(x);
            catch ME
                % If an error occurred, issue a warning and set f to NaN.
                warning(ME.identifier, '%s', ME.message);
                f = NaN;
            end

            try
                f = double(f);
            catch ME
                % If an error occurred, issue a warning and set f to NaN.
                warning(ME.identifier, '%s', ME.message);
                f = NaN;
            end
        end

        function f = cub(obj, x)
            if ~isrealvector(x)
                error("MATLAB:Problem:InvalidInputForCUB", "The input `x` for method `cub` in `Problem` must be a vector.")
            end

            if numel(x) ~= obj.n
                error("MATLAB:Problem:WrongSizeInputForCUB", "The input `x` for method `cub` in `Problem` must have size %d.", obj.n)
            end

            if isempty(obj.cub_)
                f = NaN(0, 1);
            else
                try
                    f = obj.cub_(x);
                catch ME
                    warning(ME.identifier, '%s', ME.message);
                    f = NaN(obj.m_nonlinear_ub, 1);
                end
                if ~isrealvector(f) && ~isempty(f)
                    error("MATLAB:Problem:InvalidOutputForCUB", "The output of method `cub` in `Problem` must be a real vector.")
                end
                try
                    f = double(f);
                catch ME
                    warning(ME.identifier, '%s', ME.message);
                    f = NaN;
                end
            end

            if isempty(obj.m_nonlinear_ub_)
                obj.m_nonlinear_ub_ = numel(f);
            end
            if numel(f) ~= obj.m_nonlinear_ub
                error("MATLAB:Problem:cubx_m_nonlinear_ub_NotConsistent", "The size of `cub(x)` is not consistent with `m_nonlinear_ub`=%d.", obj.m_nonlinear_ub)
            end
        end

        function f = ceq(obj, x)
            if ~isrealvector(x)
                error("MATLAB:Problem:InvalidInputForCEQ", "The input `x` for method `ceq` in `Problem` must be a real vector.")
            end

            if numel(x) ~= obj.n
                error("MATLAB:Problem:WrongSizeInputForCEQ", "The input `x` for method `ceq` in `Problem` must have size %d.", obj.n)
            end

            if isempty(obj.ceq_)
                f = NaN(0, 1);
            else
                try
                    f = obj.ceq_(x);
                catch ME
                    warning(ME.identifier, '%s', ME.message);
                    f = NaN(obj.m_nonlinear_eq, 1);
                end
                if ~isrealvector(f) && ~isempty(f)
                    error("MATLAB:Problem:InvalidOutputForCEQ", "The output of method `ceq` in `Problem` must be a vector.")
                end
                try
                    f = double(f);
                catch ME
                    warning(ME.identifier, '%s', ME.message);
                    f = NaN;
                end
            end

            if isempty(obj.m_nonlinear_eq_)
                obj.m_nonlinear_eq_ = numel(f);
            end
            if numel(f) ~= obj.m_nonlinear_eq
                error("MATLAB:Problem:ceqx_m_nonlinear_eq_NotConsistent", "The size of `ceq(x)` is not consistent with `m_nonlinear_eq`=%d.", obj.m_nonlinear_eq)
            end
        end

        function m = m_nonlinear_ub(obj)
            if isempty(obj.m_nonlinear_ub_)
                if isempty(obj.cub_)
                    m = 0;
                else
                    error("MATLAB:Problem:m_nonlinear_ub_Unknown", "The number of nonlinear inequality constraints for `Problem` is unknown.");
                end
            else
                m = obj.m_nonlinear_ub_;
            end
        end

        function m = m_nonlinear_eq(obj)
            if isempty(obj.m_nonlinear_eq_)
                if isempty(obj.ceq_)
                    m = 0;
                else
                    error("MATLAB:Problem:m_nonlinear_eq_Unknown", "The number of nonlinear equality constraints for `Problem` is unknown.");
                end
            else
                m = obj.m_nonlinear_eq_;
            end
        end

        function project_x0(obj)
            % Project the initial guess onto the feasible region.
            function val = dist_x0_sq(x)
                val = norm(x - obj.x0)^2;
            end
            function [c, ceq] = nonlcon(x)
                c = obj.cub(x);
                ceq = obj.ceq(x);
            end
            if strcmp(obj.p_type, 'unconstrained')
                return
            elseif strcmp(obj.p_type, 'bound-constrained')
                obj.x0 = min(max(obj.x0, obj.xl), obj.xu);
            elseif strcmp(obj.p_type, 'linearly constrained') && obj.m_linear_ub == 0 && all(obj.xl == -Inf) && all(obj.xu == Inf)
                try
                    [~, res] = evalc('lsqr(obj.aeq, obj.beq - obj.aeq * obj.x0)');
                    obj.x0 = obj.x0 + res;
                catch
                end
            elseif strcmp(obj.p_type, 'linearly constrained') && exist('lsqlin') == 2
                try
                    [~, res] = evalc('lsqlin(eye(obj.n), obj.x0, obj.aub, obj.bub, obj.aeq, obj.beq, obj.xl, obj.xu, obj.x0)');
                    obj.x0 = res;
                catch
                end
            elseif ~strcmp(obj.p_type, 'unconstrained') && exist('fmincon') == 2
                try
                    [~, res] = evalc('fmincon(@(x) dist_x0_sq(x), obj.x0, obj.aub, obj.bub, obj.aeq, obj.beq, obj.xl, obj.xu, @(x) nonlcon(x))');
                    obj.x0 = res;
                catch
                end
            end
        end
    end

end