classdef Problem
    %{
    Optimization problem.
    %}

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
    properties (GetAccess = public, SetAccess = private)

        fun
        x0
        xl = []
        xu = []
        aub = []
        bub = []
        aeq = []
        beq = []
        cub = []
        ceq = []
        m_nonlinear_ub = []
        m_nonlinear_eq = []

    end

    properties (Dependent)

        n
        m_linear_ub
        m_linear_eq

    end

    methods

        function obj = Problem(varargin)
            if nargin < 1
                error("At least one argument is required.")
            end
        
            if isstruct(varargin{1})
                % Handle the case when a struct is passed
                s = varargin{1};
        
                % Check if the struct contains the required fields
                if ~isfield(s, 'fun') || ~isfield(s, 'x0')
                    error("The `fun` and `x0` fields are required.")
                end
        
                % Iterate over the struct's fields and assign them to the object's properties
                fields = fieldnames(s);
                for i = 1:numel(fields)
                    obj.(fields{i}) = s.(fields{i});
                end
            else
                error("Invalid input. A struct argument is expected.")
            end

            % Check that the arguments are consistent.

            % Check that `xl` has the same size as `x0`.
            if numel(obj.xl) ~= obj.n
                error('The argument `xl` must have size %d.', obj.n);
            end
            % Check that `xu` has the same size as `x0`.
            if numel(obj.xu) ~= obj.n
                error('The argument `xu` must have size %d.', obj.n);
            end
            % Check that `aub` is a matrix with shape (m_linear_ub, n).
            if ~isequal(size(obj.aub), [obj.m_linear_ub, obj.n])
                error('The argument `aub` must have shape (%d, %d).', obj.m_linear_ub, obj.n);
            end
            % Check that `aeq` is a matrix with shape (m_linear_eq, n).
            if ~isequal(size(obj.aeq), [obj.m_linear_eq, obj.n])
                error('The argument `aeq` must have shape (%d, %d).', obj.m_linear_eq, obj.n);
            end
            % Check that whether `m_nonlinear_ub` is empty or zero if `cub` is empty.
            if isempty(obj.cub) && ~isempty(obj.m_nonlinear_ub) && obj.m_nonlinear_ub > 0
                error('The argument `m_nonlinear_ub` must be empty or zero if the argument `cub` is empty.');
            end
            % Check that whether `m_nonlinear_eq` is empty or zero if `ceq` is empty.
            if isempty(obj.ceq) && ~isempty(obj.m_nonlinear_eq) && obj.m_nonlinear_eq > 0
                error('The argument `m_nonlinear_eq` must be empty or zero if the argument `ceq` is empty.');
            end
        end



        % Setter functions.

        % Preprocess the initial guess.
        function obj = set.x0(obj, value)
            if ~isvector(value)
                error("The argument `x0` must be a one-dimensional array.")
            end
            obj.x0 = reshape(value, [], 1);
        end

        % Preprocess the objective function.
        function obj = set.fun(obj, fun)
            % This function only accepts a new function handle.
            if ~isa(fun, "function_handle")
                error("The argument `fun` must be a function handle.")
            end
            % Check if `fun` can accept only one argument.
            if nargin(fun) ~= 1
                error('The function must accept exactly one argument.')
            end
            % Try to assign `fun` to the object's `fun` field.
            try
                obj.fun = fun;
            catch
                error("Error occurred while assigning 'fun' field. Please check the input.")
            end
        end

        % Preprocess the bound constraints.
        function obj = set.xl(obj, xl)
            if ~isempty(xl)
                if ~isvector(xl)
                    error("The argument `xl` must be a one-dimensional array.")
                end
                obj.xl = reshape(xl, [], 1);
            else
                obj.xl = [];
            end
        end

        function obj = set.xu(obj, xu)
            if ~isempty(xu)
                if ~isvector(xu)
                    error("The argument `xu` must be a one-dimensional array.")
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
                    error("The argument `aub` must be a two-dimensional array.")
                end
                obj.aub = aub;
            else
                obj.aub = [];
            end
        end

        function obj = set.bub(obj, bub)
            if ~isempty(bub)
                if ~isvector(bub)
                    error("The argument `bub` must be a one-dimensional array.")
                end
                obj.bub = reshape(bub, [], 1);
            else
                obj.bub = [];
            end
        end

        function obj = set.aeq(obj, aeq)
            if ~isempty(aeq)
                if ~ismatrix(aeq)
                    error("The argument `aeq` must be a two-dimensional array.")
                end
                obj.aeq = aeq;
            else
                obj.aeq = [];
            end
        end

        function obj = set.beq(obj, beq)
            if ~isempty(beq)
                if ~isvector(beq)
                    error("The argument `beq` must be a one-dimensional array.")
                end
                obj.beq = reshape(beq, [], 1);
            else
                obj.beq = [];
            end
        end

        % Preprocess the nonlinear constraints.
        function obj = set.cub(obj, cub)
            if ~isempty(cub)
                if ~isa(cub, "function_handle")
                    error("The argument `cub` must be a function handle.")
                end
                try
                    cub(obj.x0);
                    obj.cub = cub;
                catch
                    error("The argument `cub` must accept an array with shape (%d, 1) as input.", obj.n);
                end
            else
                obj.cub = [];
            end
        end

        function obj = set.ceq(obj, ceq)
            if ~isempty(ceq)
                if ~isa(ceq, "function_handle")
                    error("The argument `ceq` must be a function handle.")
                end
                try
                    ceq(obj.x0);
                    obj.ceq = ceq;
                catch
                    error("The argument `ceq` must accept an array with shape (%d, 1) as input.", obj.n);
                end
            else
                obj.ceq = [];
            end
        end

        % Preprocess the number of nonlinear constraints.
        function obj = set.m_nonlinear_ub(obj, m_nonlinear_ub)
            if ~isempty(m_nonlinear_ub)
                if ~(isnumeric(m_nonlinear_ub) && m_nonlinear_ub >= 0 && mod(m_nonlinear_ub, 1) == 0)
                    error('The argument `m_nonlinear_ub` must be a nonnegative integer.')
                end
                obj.m_nonlinear_ub = m_nonlinear_ub;
            else
                obj.m_nonlinear_ub = [];
            end
        end

        function obj = set.m_nonlinear_eq(obj, m_nonlinear_eq)
            if ~isempty(m_nonlinear_eq)
                if ~(isnumeric(m_nonlinear_eq) && m_nonlinear_eq >= 0 && mod(m_nonlinear_eq, 1) == 0)
                    error('The argument `m_nonlinear_eq` must be a nonnegative integer.')
                end
                obj.m_nonlinear_eq = m_nonlinear_eq;
            else
                obj.m_nonlinear_eq = [];
            end
        end

        % Getter for dependent properties.

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
            if isempty(obj.m_nonlinear_ub)
                if isempty(obj.cub)
                    value = 0;
                else
                    error('The number of nonlinear inequality constraints is unknown.');
                end
            else
                value = obj.m_nonlinear_ub;
            end
        end

        function value = get.m_nonlinear_eq(obj)
            if isempty(obj.m_nonlinear_eq)
                if isempty(obj.ceq)
                    value = 0;
                else
                    error('The number of nonlinear equality constraints is unknown.');
                end
            else
                value = obj.m_nonlinear_eq;
            end
        end


        % Getter functions.

        function value = get.fun(obj)
            value = @(x) obj.checkFUN(x);
        end

        function f = checkFUN(obj, x)
            try
                % Try to evaluate the function at x.
                f = obj.fun(x);
                f = double(f);
            catch ME
                % If an error occurred, issue a warning and set f to NaN.
                warning(ME.identifier, '%s', ME.message);
                f = NaN;
            end
            if ~(isnumeric(f) || islogical(f))
                % Check if the output is a float/int/boolean. If not, set f to NaN.
                f = NaN;
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

        function value = calc_fun(obj, x)
            value = double(obj.fun(x));
        end
    
        function value = calc_cub(obj, x)
            if isempty(obj.cub)
                value = NaN(0, 1);
            else
                value = double(obj.cub(x));
            end
            if isempty(obj.m_nonlinear_ub)
                obj.m_nonlinear_ub = numel(value);
            end
        end
    
        function value = calc_ceq(obj, x)
            if isempty(obj.ceq)
                value = NaN(0, 1);
            else
                value = double(obj.ceq(x));
            end
            if isempty(obj.m_nonlinear_eq)
                obj.m_nonlinear_eq = numel(value);
            end
        end

    end

end