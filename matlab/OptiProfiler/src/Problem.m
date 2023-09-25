classdef Problem < handle
    %PROBLEM used to define test problems for optimization solvers.
    %

    properties (GetAccess = public, SetAccess = private)

        x0
        xl
        xu
        aub
        bub
        aeq
        beq

    end

    properties (GetAccess = public, SetAccess = private)

        m_nonlinear_ub
        m_nonlinear_eq

    end

    properties (Dependent)

        n
        m_linear_ub
        m_linear_eq

    end

    properties (Access = private)

        fun_
        cub_
        ceq_

    end

    methods

        %{
            Initialize an optimization problem.

            Parameters
            ----------
            fun : callable
                Objective function to be minimized.

                    ``fun(x) -> float``

                where ``x`` is an array with shape (n,).
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

            Raises
            ------
            ValueError
                If an argument received an invalid value.
        %}
        function obj = Problem(varargin)
            if nargin < 1
                error("MATLAB:Problem:MissingArguments", "At least one argument is required.")
            end
        
            if isstruct(varargin{1})
                % Handle the case when a struct is passed
                s = varargin{1};
        
                % Check if the struct contains the required fields
                if ~isfield(s, 'fun') || ~isfield(s, 'x0')
                    error("MATLAB:Problem:MissingFields", "The `fun` and `x0` fields are required.")
                end
                obj.fun_ = s.fun;

                % Check if the struct contains `cub` and `ceq` fields
                if isfield(s, 'cub')
                    obj.cub_ = s.cub;
                end
                if isfield(s, 'ceq')
                    obj.ceq_ = s.ceq;
                end
        
                % Iterate over the struct's fields and assign them to the object's properties
                expected_fields = {'fun', 'x0', 'xl', 'xu', 'aub', 'bub', 'aeq', 'beq', 'cub', 'ceq', 'm_nonlinear_ub', 'm_nonlinear_eq'};
                fields = fieldnames(s);
                for i = 1:numel(expected_fields)
                    if strcmp(expected_fields{i}, 'fun') || strcmp(expected_fields{i}, 'cub') || strcmp(expected_fields{i}, 'ceq')
                        continue
                    elseif ismember(expected_fields{i}, fields)
                        obj.(expected_fields{i}) = s.(expected_fields{i});
                    else
                        obj.(expected_fields{i}) = [];
                    end
                end
            else
                error("MATLAB:Problem:NotStruct", "Invalid input. A struct argument is expected.")
            end

            % Check that the arguments are consistent.

            % Check that `xl` has the same size as `x0`.
            if numel(obj.xl) ~= obj.n
                error("MATLAB:Problem:xl_x0_NotConsistent", "The argument `xl` must have size %d.", obj.n);
            end
            % Check that `xu` has the same size as `x0`.
            if numel(obj.xu) ~= obj.n
                error("MATLAB:Problem:xu_x0_NotConsistent", "The argument `xu` must have size %d.", obj.n);
            end
            % Check that `aub` is a matrix with shape (m_linear_ub, n).
            if ~isequal(size(obj.aub), [obj.m_linear_ub, obj.n])
                error("MATLAB:Problem:aub_m_linear_ub_n_NotConsistent", "The argument `aub` must have shape (%d, %d).", obj.m_linear_ub, obj.n);
            end
            % Check that `aeq` is a matrix with shape (m_linear_eq, n).
            if ~isequal(size(obj.aeq), [obj.m_linear_eq, obj.n])
                error("MATLAB:Problem:aeq_m_linear_eq_n_NotConsistent", "The argument `aeq` must have shape (%d, %d).", obj.m_linear_eq, obj.n);
            end
            % Check that whether `m_nonlinear_ub` is empty or zero if `cub` is empty.
            if isempty(obj.cub_) && ~isempty(obj.m_nonlinear_ub) && obj.m_nonlinear_ub > 0
                error("MATLAB:Problem:m_nonlinear_ub_cub_NotConsistent", "The argument `m_nonlinear_ub` must be empty or zero if the argument `cub` is empty.");
            end
            % Check that whether `m_nonlinear_eq` is empty or zero if `ceq` is empty.
            if isempty(obj.ceq_) && ~isempty(obj.m_nonlinear_eq) && obj.m_nonlinear_eq > 0
                error("MATLAB:Problem:m_nonlinear_eq_ceq_NotConsistent", "The argument `m_nonlinear_eq` must be empty or zero if the argument `ceq` is empty.");
            end
        end



        % Setter functions.

        % Preprocess the initial guess.
        function obj = set.x0(obj, value)
            if ~isvector(value)
                error("MATLAB:Problem:x0_NotVector", "The argument `x0` must be a vector.")
            end
            obj.x0 = reshape(value, [], 1);
        end

        % Preprocess the objective function.
        function obj = set.fun_(obj, fun)
            % This function only accepts a new function handle.
            if ~isa(fun, "function_handle")
                error("MATLAB:Problem:fun_NotFunctionHandle", "The argument `fun` must be a function handle.")
            end
            % Check if `fun` can accept only one argument.
            if nargin(fun) ~= 1
                error("MATLAB:Problem:fun_NotOneArguementFun", "The function must accept exactly one argument.")
            end
            
            obj.fun_ = fun;
        end

        % Preprocess the bound constraints.
        function obj = set.xl(obj, xl)
            if ~isempty(xl)
                if ~isvector(xl)
                    error("MATLAB:Problem:xl_NotVector", "The argument `xl` must be a vector.")
                end
                obj.xl = reshape(xl, [], 1);
            else
                obj.xl = [];
            end
        end

        function obj = set.xu(obj, xu)
            if ~isempty(xu)
                if ~isvector(xu)
                    error("MATLAB:Problem:xu_NotVector", "The argument `xu` must be a vector.")
                end
                obj.xu = reshape(xu, [], 1);
            else
                obj.xu = [];
            end
        end

        % Preprocess the linear constraints.
        function obj = set.aub(obj, aub)
            if ~isempty(aub)
                if ~ismatrix(aub)
                    error("MATLAB:Problem:aub_NotMatrix", "The argument `aub` must be a matrix.")
                end
                obj.aub = aub;
            else
                obj.aub = [];
            end
        end

        function obj = set.bub(obj, bub)
            if ~isempty(bub)
                if ~isvector(bub)
                    error("MATLAB:Problem:bub_NotVector", "The argument `bub` must be a vector.")
                end
                obj.bub = reshape(bub, [], 1);
            else
                obj.bub = [];
            end
        end

        function obj = set.aeq(obj, aeq)
            if ~isempty(aeq)
                if ~ismatrix(aeq)
                    error("MATLAB:Problem:aeq_NotMatrix", "The argument `aeq` must be a matrix.")
                end
                obj.aeq = aeq;
            else
                obj.aeq = [];
            end
        end

        function obj = set.beq(obj, beq)
            if ~isempty(beq)
                if ~isvector(beq)
                    error("MATLAB:Problem:beq_NotVector", "The argument `beq` must be a vector.")
                end
                obj.beq = reshape(beq, [], 1);
            else
                obj.beq = [];
            end
        end

        % Preprocess the nonlinear constraints.
        function obj = set.cub_(obj, cub)
            if ~isa(cub, "function_handle") && ~isempty(cub)
                error("MATLAB:Problem:cub_NotFunctionHandle", "The argument `cub` must be a function handle.")
            end
            
            obj.cub_ = cub;
        end

        function obj = set.ceq_(obj, ceq)
            if ~isa(ceq, "function_handle") && ~isempty(ceq)
                error("MATLAB:Problem:ceq_NotFunctionHandle", "The argument `ceq` must be a function handle.")
            end
            
            obj.ceq_ = ceq;
        end

        % Preprocess the number of nonlinear constraints.
        function obj = set.m_nonlinear_ub(obj, m_nonlinear_ub)
            if ~isempty(m_nonlinear_ub)
                if ~(isnumeric(m_nonlinear_ub) && m_nonlinear_ub >= 0 && mod(m_nonlinear_ub, 1) == 0)
                    error("MATLAB:Problem:m_nonlinear_ub_NotPositiveScalar", "The argument `m_nonlinear_ub` must be a nonnegative integer.")
                end
                obj.m_nonlinear_ub = m_nonlinear_ub;
            else
                obj.m_nonlinear_ub = [];
            end
        end

        function obj = set.m_nonlinear_eq(obj, m_nonlinear_eq)
            if ~isempty(m_nonlinear_eq)
                if ~(isnumeric(m_nonlinear_eq) && m_nonlinear_eq >= 0 && mod(m_nonlinear_eq, 1) == 0)
                    error("MATLAB:Problem:m_nonlinear_eq_NotPositiveScalar", "The argument `m_nonlinear_eq` must be a nonnegative integer.")
                end
                obj.m_nonlinear_eq = m_nonlinear_eq;
            else
                obj.m_nonlinear_eq = [];
            end
        end

        % Getter functions for dependent properties.

        function value = get.n(obj)
            value = numel(obj.x0);
        end

        function value = get.m_linear_ub(obj)
            value = numel(obj.bub);
        end

        function value = get.m_linear_eq(obj)
            value = numel(obj.beq);
        end

        % Other getter functions.

        function value = get.m_nonlinear_ub(obj)
            try
                if isempty(obj.m_nonlinear_ub)
                    if isempty(obj.cub_)
                        value = 0;
                    else
                        error("MATLAB:Problem:m_nonlinear_ub_Unknown", "The number of nonlinear inequality constraints is unknown.");
                    end
                else
                    value = obj.m_nonlinear_ub;
                end
            catch ME
                warning("Error occurred while getting m_nonlinear_ub: %s", ME.message);
                value = [];
            end
        end

        function value = get.m_nonlinear_eq(obj)
            try
                if isempty(obj.m_nonlinear_eq)
                    if isempty(obj.ceq_)
                        value = 0;
                    else
                        error("MATLAB:Problem:m_nonlinear_eq_Unknown", "The number of nonlinear equality constraints is unknown.");
                    end
                else
                    value = obj.m_nonlinear_eq;
                end
            catch ME
                warning("Error occurred while getting m_nonlinear_eq: %s", ME.message);
                value = [];
            end
        end        

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

        function cv = maxcv(obj, x)
            % Evaluate the maximum constraint violation.
            
            % Initialize cv to 0
            cv = 0;
            
            % Check lower bound constraint violation
            if ~isempty(obj.xl)
                cv = max(max(obj.xl - x), cv);
            end
            
            % Check upper bound constraint violation
            if ~isempty(obj.xu)
                cv = max(max(x - obj.xu), cv);
            end
            
            % Check linear inequality constraint violation
            if ~isempty(obj.aub)
                cv = max(max(obj.aub * x - obj.bub), cv);
            end
            
            % Check linear equality constraint violation
            if ~isempty(obj.aeq)
                cv = max(max(abs(obj.aeq * x - obj.beq)), cv);
            end
            
            % Check nonlinear inequality constraint violation
            if ~isempty(obj.cub_)
                cv = max(max(obj.cub(x)), cv);
            end
            
            % Check nonlinear equality constraint violation
            if ~isempty(obj.ceq_)
                cv = max(max(abs(obj.ceq(x))), cv);
            end
        end

        function f = fun(obj, x)
            % Check if x is a vector.
            if ~isvector(x)
                error("MATLAB:Problem:InvalidInputForFUN", "The input `x` must be a vector.")
            end

            if numel(x) ~= obj.n
                error("MATLAB:Problem:WrongSizeInputForFUN", "The input `x` must have size %d.", obj.n)
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
            if ~isvector(x)
                error("MATLAB:Problem:InvalidInputForCUB", "The input `x` must be a vector.")
            end

            if numel(x) ~= obj.n
                error("MATLAB:Problem:WrongSizeInputForCUB", "The input `x` must have size %d.", obj.n)
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
                if ~isvector(f)
                    error("MATLAB:Problem:InvalidOutputForCUB", "The output of the nonlinear inequality constraint must be a vector.")
                end
                try
                    f = double(f);
                catch ME
                    warning(ME.identifier, '%s', ME.message);
                    f = NaN;
                end
            end

            if isempty(obj.m_nonlinear_ub)
                obj.m_nonlinear_ub = numel(f);
            end
            if numel(f) ~= obj.m_nonlinear_ub
                error("MATLAB:Problem:cubx_m_nonlinear_ub_NotConsistent", "The size of `cub(x)` is not consistent with `m_nonlinear_ub`=%d.", obj.m_nonlinear_ub)
            end
        end

        function f = ceq(obj, x)
            if ~isvector(x)
                error("MATLAB:Problem:InvalidInputForCEQ", "The input `x` must be a vector.")
            end

            if numel(x) ~= obj.n
                error("MATLAB:Problem:WrongSizeInputForCEQ", "The input `x` must have size %d.", obj.n)
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
                if ~isvector(f)
                    error("MATLAB:Problem:InvalidOutputForCEQ", "The output of the nonlinear equality constraint must be a vector.")
                end
                try
                    f = double(f);
                catch ME
                    warning(ME.identifier, '%s', ME.message);
                    f = NaN;
                end
            end

            if isempty(obj.m_nonlinear_eq)
                obj.m_nonlinear_eq = numel(f);
            end
            if numel(f) ~= obj.m_nonlinear_eq
                error("MATLAB:Problem:ceqx_m_nonlinear_eq_NotConsistent", "The size of `ceq(x)` is not consistent with `m_nonlinear_eq`=%d.", obj.m_nonlinear_eq)
            end
        end
    end

end