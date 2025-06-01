classdef Problem < handle
%PROBLEM is a class that defines an optimization problem.
%
%   PROBLEM describes an optimization problem with the following structure:
%
%       min fun(x)
%       s.t. xl <= x <= xu,
%            aub * x <= bub,
%            aeq * x = beq,
%            cub(x) <= 0,
%            ceq(x) = 0,
%       with initial point x0,
%
%   where `fun` is the objective function, `x` is the variable to optimize,
%   `xl` and `xu` are the lower and upper bounds, `aub` and `bub` are the
%   coefficient matrix and right-hand side vector of the linear inequality
%   constraints, `aeq` and `beq` are the coefficient matrix and right-hand side
%   vector of the linear equality constraints, `cub` is the function of
%   nonlinear inequality constraints, `ceq` is the function of nonlinear
%   equality constraints.
%
%   PROBLEM should be initialized by the following signature:
%
%       P = PROBLEM(P_STRUCT);
%
%   where the return P is an instance of the class PROBLEM and the input
%   P_STRUCT is a struct.
%
%   The input struct P_STRUCT should contain the following fields:
%
%       1. compulsory fields:
%       
%       - fun: the objective function ``fun(x) -> float``.
%       - x0: the initial guess.
%
%       2. optional fields:
%
%       - name: the name of the problem. It should be a string or a char.
%       Default is 'Unnamed Problem'.
%       - xl: the lower bounds on the variable `x` in the form of ``xl <= x``.
%         It should be a vector of the same size as `x`. Default is -Inf.
%       - xu: the upper bounds on the variable `x` in the form of ``xl >= x``.
%         It should be a vector of the same size as `x`. Default is Inf.
%       - aub, bub: the coefficient matrix and right-hand side vector of the
%         linear inequality constraints ``aub * x <= bub``. The default setting
%         of `aub` and `bub` are empty matrix and vector.
%       - aeq, beq: the coefficient matrix and right-hand side vector of the
%         linear equality constraints ``aeq * x <= beq``. The default setting
%         of `aeq` and `beq` are empty matrix and vector.
%       - cub: the function of nonlinearly inequality constraints
%         ``cub(x) <= 0``, where ``cub(x) -> float vector``. By default, 
%         ``cub(x)`` will return an empty vector.
%       - ceq: the function of nonlinearly equality constraints
%         ``ceq(x) <= 0``, where ``ceq(x) -> float vector``. By default, 
%         ``ceq(x)`` will return an empty vector.
%       - grad: the gradient of the objective function
%         ``grad(x) -> float vector``. By default, `grad(x)` will return an
%         empty vector.
%       - hess: the Hessian of the objective function
%         ``hess(x) -> float matrix``. By default, `hess(x)` will return an
%         empty matrix.
%       - jcub: the Jacobian of the nonlinearly inequality constraints
%         ``jcub(x) -> float matrix``. Note that the column size of `jcub(x)`
%         should be the same as the length of `x` while the row size should be
%         the same as the length of `cub(x)`. By default, `jcub(x)` will return
%         an empty matrix.
%       - jceq: the Jacobian of the nonlinearly equality constraints
%         ``jceq(x) -> float matrix``. Note that the column size of `jceq(x)`
%         should be the same as the length of `x` while the row size should be
%         the same as the length of `ceq(x)`. By default, `jceq(x)` will return
%         an empty matrix.
%       - hcub: the Hessian of the nonlinearly inequality constraints
%         ``hcub(x) -> cell array of float matrices``. The i-th element of
%         `hcub(x)` should be the Hessian of the i-th function in `cub`. By
%         default, `hcub(x)` will return an empty cell.
%       - hceq: the Hessian of the nonlinearly equality constraints
%         ``hceq(x) -> cell array of float matrices``. The i-th element of
%         `hceq(x)` should be the Hessian of the i-th function in `ceq`. By
%         default, `hceq(x)` will return an empty cell.
%
%   The output P contains following properties:
%
%       1. properties inherited from the input struct:
%
%       name, x0, xl, xu, aub, bub, aeq, beq
%
%       2. properties dependent on the input struct:
%
%       - ptype: type of the problem. It should be 'u' (unconstrained), 'b'
%         (bound-constrained), 'l' (linearly constrained), or 'n' (nonlinearly
%         constrained).
%       - n: dimension of the problem, which is the length of the variable `x`.
%       - mb: number of the bound constraints, which is the length of finite
%         elements in `xl` and `xu`.
%       - m_linear_ub: number of the linear inequality constraints, which is
%         the length of `bub`.
%       - m_linear_eq: number of the linear equality constraints, which is the
%         length of `beq`.
%       - m_nonlinear_ub: number of the nonlinear inequality constraints, which
%         is the length of `cub(x)`.
%       - m_nonlinear_eq: number of the nonlinear equality constraints, which
%         is the length of `ceq(x)`.
%       - mlcon: number of the linear constraints, which is the sum of
%         `m_linear_ub` and `m_linear_eq`.
%       - mnlcon: number of the nonlinear constraints, which is the sum of
%         `m_nonlinear_ub` and `m_nonlinear_eq`.
%       - mcon: number of the constraints, which is the sum of `mlcon` and
%         `mnlcon`.
%
%   The output P contains following methods:
%
%       1. methods inherited from the input struct:
%
%       fun, grad, hess, cub, ceq, jcub, jceq, hcub, hceq
%
%       2. other methods:
%
%       - maxcv: the maximum constraint violation ``maxcv(x) -> float``, which
%         is defined as the maximum of the infinity norms of
%         `max(xl - x, 0)`, `max(x - xu, 0)`, `max(aub * x - bub, 0)`,
%         `aeq * x - beq`, `max(cub(x), 0)`, `ceq(x)`.
%       - project_x0: trying to project the initial guess `x0` onto the
%         feasible region if it is not feasible (but it may fail).
%


    properties (GetAccess = public, SetAccess = public)

        name = 'Unnamed Problem'
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
        mb
        m_linear_ub
        m_linear_eq
        m_nonlinear_ub
        m_nonlinear_eq
        mlcon
        mnlcon
        mcon
        ptype

    end

    properties (Access = protected)

        fun_
        grad_
        hess_
        cub_
        ceq_
        jcub_
        jceq_
        hcub_
        hceq_

    end

    methods

        % Initialize a PROBLEM.
        function obj = Problem(varargin)

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

                % Check if the struct contains `grad` and `hess` fields
                if isfield(s, 'grad')
                    obj.grad_ = s.grad;
                end
                if isfield(s, 'hess')
                    obj.hess_ = s.hess;
                end

                % Check if the struct contains `cub` and `ceq` fields
                if isfield(s, 'cub')
                    obj.cub_ = s.cub;
                end
                if isfield(s, 'ceq')
                    obj.ceq_ = s.ceq;
                end

                % Check if the struct contains `jcub`, `jceq`, `hcub`, and `hceq` fields
                if isfield(s, 'jcub')
                    obj.jcub_ = s.jcub;
                end
                if isfield(s, 'jceq')
                    obj.jceq_ = s.jceq;
                end
                if isfield(s, 'hcub')
                    obj.hcub_ = s.hcub;
                end
                if isfield(s, 'hceq')
                    obj.hceq_ = s.hceq;
                end
        
                % Iterate over the struct's fields and assign them to the object's properties
                expected_fields = {'name', 'fun', 'grad', 'hess', 'x0', 'xl', 'xu', 'aub', 'bub', 'aeq', 'beq', 'cub', 'ceq', 'jcub', 'jceq', 'hcub', 'hceq'};
                fields = fieldnames(s);
                for i = 1:numel(expected_fields)
                    if strcmp(expected_fields{i}, 'fun') || strcmp(expected_fields{i}, 'grad') || strcmp(expected_fields{i}, 'hess') || strcmp(expected_fields{i}, 'cub') || strcmp(expected_fields{i}, 'ceq') || strcmp(expected_fields{i}, 'jcub') || strcmp(expected_fields{i}, 'jceq') || strcmp(expected_fields{i}, 'hcub') || strcmp(expected_fields{i}, 'hceq')
                        continue
                    elseif ismember(expected_fields{i}, fields)
                        obj.(expected_fields{i}) = s.(expected_fields{i});
                    end
                end
            else
                error("MATLAB:Problem:NotStruct", "Invalid input for `Problem`. A struct argument is expected.")
            end

            % Check that the arguments are legal and consistent.

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
        end

        % Setter functions.

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

        % Preprocess the gradient and the Hessian of the objective function.
        function set.grad_(obj, grad_)
            if ~isa(grad_, 'function_handle') && ~isempty(grad_)
                error("MATLAB:Problem:grad_NotFunctionHandle", "The argument `grad` for `Problem` must be a function handle.")
            end
            obj.grad_ = grad_;
        end

        function set.hess_(obj, hess_)
            if ~isa(hess_, 'function_handle') && ~isempty(hess_)
                error("MATLAB:Problem:hess_NotFunctionHandle", "The argument `hess` for `Problem` must be a function handle.")
            end
            obj.hess_ = hess_;
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

        % Preprocess the Jacobian and the Hessian of the nonlinear constraints.
        function set.jcub_(obj, jcub_)
            if ~isa(jcub_, 'function_handle') && ~isempty(jcub_)
                error("MATLAB:Problem:jcub_NotFunctionHandle", "The argument `jcub` for `Problem` must be a function handle.")
            end
            obj.jcub_ = jcub_;
        end

        function set.jceq_(obj, jceq_)
            if ~isa(jceq_, 'function_handle') && ~isempty(jceq_)
                error("MATLAB:Problem:jceq_NotFunctionHandle", "The argument `jceq` for `Problem` must be a function handle.")
            end
            obj.jceq_ = jceq_;
        end

        function set.hcub_(obj, hcub_)
            if ~isa(hcub_, 'function_handle') && ~isempty(hcub_)
                error("MATLAB:Problem:hcub_NotFunctionHandle", "The argument `hcub` for `Problem` must be a function handle.")
            end
            obj.hcub_ = hcub_;
        end

        function set.hceq_(obj, hceq_)
            if ~isa(hceq_, 'function_handle') && ~isempty(hceq_)
                error("MATLAB:Problem:hceq_NotFunctionHandle", "The argument `hceq` for `Problem` must be a function handle.")
            end
            obj.hceq_ = hceq_;
        end

        % Getter functions for dependent properties.

        function value = get.n(obj)
            value = numel(obj.x0);
        end

        function value = get.mb(obj)
            value = sum(~isinf(-obj.xl)) + sum(~isinf(obj.xu));
        end

        function value = get.m_linear_ub(obj)
            value = sum(~isinf(obj.bub));
        end

        function value = get.m_linear_eq(obj)
            value = numel(obj.beq);
        end

        function value = get.m_nonlinear_ub(obj)
            try
                cub = obj.cub_(obj.x0);
                value = numel(cub);
            catch ME
                value = 0;
            end
        end

        function value = get.m_nonlinear_eq(obj)
            try
                ceq = obj.ceq_(obj.x0);
                value = numel(ceq);
            catch ME
                value = 0;
            end
        end

        function value = get.mlcon(obj)
            value = obj.m_linear_ub + obj.m_linear_eq;
        end

        function value = get.mnlcon(obj)
            value = obj.m_nonlinear_ub + obj.m_nonlinear_eq;
        end

        function value = get.mcon(obj)
            value = obj.mlcon + obj.mnlcon;
        end

        function value = get.ptype(obj)
            try
                if obj.m_nonlinear_ub + obj.m_nonlinear_eq > 0
                    value = 'n';
                elseif obj.m_linear_ub + obj.m_linear_eq > 0
                    value = 'l';
                elseif any(obj.xl > -inf) || any(obj.xu < inf)
                    value = 'b';
                else
                    value = 'u';
                end
            catch ME
                value = 'n';
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

            if strcmp(obj.ptype, 'u')
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
            if strcmp(obj.ptype, 'b')
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
            if strcmp(obj.ptype, 'l')
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

        function g = grad(obj, x)
            % Check if x is a vector.
            if ~isrealvector(x)
                error("MATLAB:Problem:InvalidInputForGRAD", "The input `x` for method `grad` in `Problem` must be a vector.")
            end

            if numel(x) ~= obj.n
                error("MATLAB:Problem:WrongSizeInputForGRAD", "The input `x` for method `grad` in `Problem` must have size %d.", obj.n)
            end

            if isempty(obj.grad_)
                g = NaN(0, 1);
            else
                try
                    % Try to evaluate the gradient at x.
                    g = obj.grad_(x);
                catch ME
                    % If an error occurred, issue a warning and set g to NaN.
                    warning(ME.identifier, '%s', ME.message);
                    g = NaN(obj.n, 1);
                end
                if ~(isrealvector(g) && numel(g) == obj.n) && ~isempty(g)
                    error("MATLAB:Problem:InvalidOutputForGRAD", "The output of method `grad` in `Problem` must be a real vector having size %d.", obj.n)
                end
                try
                    g = double(g);
                catch ME
                    % If an error occurred, issue a warning and set g to NaN.
                    warning(ME.identifier, '%s', ME.message);
                    g = NaN(obj.n, 1);
                end
            end
        end

        function h = hess(obj, x)
            % Check if x is a vector.
            if ~isrealvector(x)
                error("MATLAB:Problem:InvalidInputForHESS", "The input `x` for method `hess` in `Problem` must be a vector.")
            end

            if numel(x) ~= obj.n
                error("MATLAB:Problem:WrongSizeInputForHESS", "The input `x` for method `hess` in `Problem` must have size %d.", obj.n)
            end

            if isempty(obj.hess_)
                h = NaN(0, 1);
            else
                try
                    % Try to evaluate the Hessian at x.
                    h = obj.hess_(x);
                catch ME
                    % If an error occurred, issue a warning and set h to NaN.
                    warning(ME.identifier, '%s', ME.message);
                    h = NaN(obj.n, obj.n);
                end
                if ~(isrealmatrix(h) && isequal(size(h), [obj.n, obj.n])) && ~isempty(h)
                    error("MATLAB:Problem:InvalidOutputForHESS", "The output of method `hess` in `Problem` must be a real matrix having shape (%d, %d).", obj.n, obj.n)
                end
                try
                    h = double(h);
                catch ME
                    % If an error occurred, issue a warning and set h to NaN.
                    warning(ME.identifier, '%s', ME.message);
                    h = NaN(obj.n, obj.n);
                end
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

            if numel(f) ~= obj.m_nonlinear_ub
                error("MATLAB:Problem:cubx_m_nonlinear_ub_NotConsistent", "The size of cub(x) is not consistent with m_nonlinear_ub=%d.\n\nCurrent x:\n%s\n\nCurrent cub(x):\n%s", obj.m_nonlinear_ub, mat2str(x(:), 16), mat2str(f(:), 16));
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

            if numel(f) ~= obj.m_nonlinear_eq
                error("MATLAB:Problem:ceqx_m_nonlinear_eq_NotConsistent", "The size of `ceq(x)` is not consistent with `m_nonlinear_eq`=%d.\n\nCurrent x:\n%s\n\nCurrent ceq(x):\n%s", obj.m_nonlinear_eq, mat2str(x(:), 16), mat2str(f(:), 16));
            end
        end

        function J = jcub(obj, x)
            if ~isrealvector(x)
                error("MATLAB:Problem:InvalidInputForJCUB", "The input `x` for method `jcub` in `Problem` must be a real vector.")
            end

            if numel(x) ~= obj.n
                error("MATLAB:Problem:WrongSizeInputForJCUB", "The input `x` for method `jcub` in `Problem` must have size %d.", obj.n)
            end

            if isempty(obj.jcub_)
                J = NaN(0, 1);
            else
                try
                    J = obj.jcub_(x);
                catch ME
                    warning(ME.identifier, '%s', ME.message);
                    J = NaN(obj.m_nonlinear_ub, obj.n);
                end
                if ~isrealmatrix(J) && ~isempty(J)
                    error("MATLAB:Problem:InvalidOutputForJCUB", "The output of method `jcub` in `Problem` must be a real matrix.")
                end
                if ~isequal(size(J), [obj.m_nonlinear_ub, obj.n]) && ~isempty(J)
                    error("MATLAB:Problem:jcubx_m_nonlinear_ub_n_NotConsistent", "The shape of `jcub(x)` is not consistent with `m_nonlinear_ub`=%d and `n`=%d.", obj.m_nonlinear_ub, obj.n)
                end
                try
                    J = double(J);
                catch ME
                    warning(ME.identifier, '%s', ME.message);
                    J = NaN(obj.m_nonlinear_ub, obj.n);
                end
            end
        end

        function J = jceq(obj, x)
            if ~isrealvector(x)
                error("MATLAB:Problem:InvalidInputForJCEQ", "The input `x` for method `jceq` in `Problem` must be a real vector.")
            end

            if numel(x) ~= obj.n
                error("MATLAB:Problem:WrongSizeInputForJCEQ", "The input `x` for method `jceq` in `Problem` must have size %d.", obj.n)
            end

            if isempty(obj.jceq_)
                J = NaN(0, 1);
            else
                try
                    J = obj.jceq_(x);
                catch ME
                    warning(ME.identifier, '%s', ME.message);
                    J = NaN(obj.m_nonlinear_eq, obj.n);
                end
                if ~isrealmatrix(J) && ~isempty(J)
                    error("MATLAB:Problem:InvalidOutputForJCEQ", "The output of method `jceq` in `Problem` must be a real matrix.")
                end
                if ~isequal(size(J), [obj.m_nonlinear_eq, obj.n]) && ~isempty(J)
                    error("MATLAB:Problem:jceqx_m_nonlinear_eq_n_NotConsistent", "The shape of `jceq(x)` is not consistent with `m_nonlinear_eq`=%d and `n`=%d.", obj.m_nonlinear_eq, obj.n)
                end
                try
                    J = double(J);
                catch ME
                    warning(ME.identifier, '%s', ME.message);
                    J = NaN(obj.m_nonlinear_eq, obj.n);
                end
            end
        end

        function H = hcub(obj, x)
            % H is a cell whose i-th element is the Hessian of the i-th nonlinear inequality constraint.

            if ~isrealvector(x)
                error("MATLAB:Problem:InvalidInputForHCUB", "The input `x` for method `hcub` in `Problem` must be a real vector.")
            end

            if numel(x) ~= obj.n
                error("MATLAB:Problem:WrongSizeInputForHCUB", "The input `x` for method `hcub` in `Problem` must have size %d.", obj.n)
            end

            if isempty(obj.hcub_)
                H = cell(0, 1);
            else
                try
                    H = obj.hcub_(x);
                catch ME
                    warning(ME.identifier, '%s', ME.message);
                    H = cell(obj.m_nonlinear_ub, 1);
                    for i = 1:obj.m_nonlinear_ub
                        H{i} = NaN(obj.n, obj.n);
                    end
                end
                if ~(iscell(H) && numel(H) == obj.m_nonlinear_ub) && ~isempty(H)
                    error("MATLAB:Problem:InvalidOutputForHCUB", "The output of method `hcub` in `Problem` must be a cell with %d elements.", obj.m_nonlinear_ub)
                end
                for i = 1:obj.m_nonlinear_ub
                    if ~isrealmatrix(H{i}) && ~isempty(H{i})
                        error("MATLAB:Problem:InvalidOutputForHCUB", "The output of method `hcub` in `Problem` must be a cell whose i-th element is a real matrix.")
                    end
                    if ~isequal(size(H{i}), [obj.n, obj.n]) && ~isempty(H{i})
                        error("MATLAB:Problem:hcubx_m_nonlinear_ub_n_NotConsistent", "The shape of the i-th Hessian in `hcub(x)` is not consistent with `n`=%d.", obj.n)
                    end
                    try
                        H{i} = double(H{i});
                    catch ME
                        warning(ME.identifier, '%s', ME.message);
                        H{i} = NaN(obj.n, obj.n);
                    end
                end
            end
        end

        function H = hceq(obj, x)
            % H is a cell whose i-th element is the Hessian of the i-th nonlinear equality constraint.

            if ~isrealvector(x)
                error("MATLAB:Problem:InvalidInputForHCEQ", "The input `x` for method `hceq` in `Problem` must be a real vector.")
            end

            if numel(x) ~= obj.n
                error("MATLAB:Problem:WrongSizeInputForHCEQ", "The input `x` for method `hceq` in `Problem` must have size %d.", obj.n)
            end

            if isempty(obj.hceq_)
                H = cell(0, 1);
            else
                try
                    H = obj.hceq_(x);
                catch ME
                    warning(ME.identifier, '%s', ME.message);
                    H = cell(obj.m_nonlinear_eq, 1);
                    for i = 1:obj.m_nonlinear_eq
                        H{i} = NaN(obj.n, obj.n);
                    end
                end
                if ~(iscell(H) && numel(H) == obj.m_nonlinear_eq) && ~isempty(H)
                    error("MATLAB:Problem:InvalidOutputForHCEQ", "The output of method `hceq` in `Problem` must be a cell with %d elements.", obj.m_nonlinear_eq)
                end
                for i = 1:obj.m_nonlinear_eq
                    if ~isrealmatrix(H{i}) && ~isempty(H{i})
                        error("MATLAB:Problem:InvalidOutputForHCEQ", "The output of method `hceq` in `Problem` must be a cell whose i-th element is a real matrix.")
                    end
                    if ~isequal(size(H{i}), [obj.n, obj.n]) && ~isempty(H{i})
                        error("MATLAB:Problem:hceqx_m_nonlinear_eq_n_NotConsistent", "The shape of the i-th Hessian in `hceq(x)` is not consistent with `n`=%d.", obj.n)
                    end
                    try
                        H{i} = double(H{i});
                    catch ME
                        warning(ME.identifier, '%s', ME.message);
                        H{i} = NaN(obj.n, obj.n);
                    end
                end
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
            if strcmp(obj.ptype, 'u')
                return
            elseif strcmp(obj.ptype, 'b')
                obj.x0 = min(max(obj.x0, obj.xl), obj.xu);
            elseif strcmp(obj.ptype, 'l') && obj.m_linear_ub == 0 && all(obj.xl == -Inf) && all(obj.xu == Inf)
                try
                    [~, res] = evalc('lsqr(obj.aeq, obj.beq - obj.aeq * obj.x0)');
                    obj.x0 = obj.x0 + res;
                catch
                end
            elseif strcmp(obj.ptype, 'l') && exist('lsqlin') == 2
                try
                    [~, res] = evalc('lsqlin(eye(obj.n), obj.x0, obj.aub, obj.bub, obj.aeq, obj.beq, obj.xl, obj.xu, obj.x0)');
                    obj.x0 = res;
                catch
                end
            elseif ~strcmp(obj.ptype, 'u') && exist('fmincon') == 2
                try
                    [~, res] = evalc('fmincon(@(x) dist_x0_sq(x), obj.x0, obj.aub, obj.bub, obj.aeq, obj.beq, obj.xl, obj.xu, @(x) nonlcon(x))');
                    obj.x0 = res;
                catch
                end
            end
        end
    end

end