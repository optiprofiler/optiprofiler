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
    properties

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

        function obj = Problem(fun, x0, varargin)
            
            if nargin < 2
                error("The arguments `fun` and `x0` are required.")
            else
                obj.fun = fun;
                obj.x0 = reshape(x0, [], 1);
            end
            if nargin > 2
                obj.xl = varargin{1};
            end
            if nargin > 3
                obj.xu = varargin{2};
            end
            if nargin > 4
                obj.aub = varargin{3};
            end
            if nargin > 5
                obj.bub = varargin{4};
            end
            if nargin > 6
                obj.aeq = varargin{5};
            end
            if nargin > 7
                obj.beq = varargin{6};
            end
            if nargin > 8
                obj.cub = varargin{7};
            end
            if nargin > 9
                obj.ceq = varargin{8};
            end
            if nargin > 10
                obj.m_nonlinear_ub = varargin{9};
            end
            if nargin > 11
                obj.m_nonlinear_eq = varargin{10};
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
            if ~isa(fun, "function_handle")
                error("The argument `fun` must be a function handle.")
            end
            try
                fun(obj.x0);
                obj.fun = fun;
            catch
                error("The argument `fun` must accept an array with shape (%d, 1) as input.", obj.n);
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
        function obj = set._m_nonlinear_ub(obj, m_nonlinear_ub)
            if ~isempty(m_nonlinear_ub)
                if ~(isnumeric(m_nonlinear_ub) && m_nonlinear_ub >= 0 && mod(m_nonlinear_ub, 1) == 0)
                    error('The argument `m_nonlinear_ub` must be a nonnegative integer.')
                end
                obj._m_nonlinear_ub = m_nonlinear_ub;
            else
                obj._m_nonlinear_ub = [];
            end
        end

        function obj = set._m_nonlinear_eq(obj, m_nonlinear_eq)
            if ~isempty(m_nonlinear_eq)
                if ~(isnumeric(m_nonlinear_eq) && m_nonlinear_eq >= 0 && mod(m_nonlinear_eq, 1) == 0)
                    error('The argument `m_nonlinear_eq` must be a nonnegative integer.')
                end
                obj._m_nonlinear_eq = m_nonlinear_eq;
            else
                obj._m_nonlinear_eq = [];
            end
        end

        % Setter for dependent properties.

        function value = set.n(obj)
            value = numel(obj.x0);
        end

        function value = set.m_linear_ub(obj)
            value = numel(obj.bub);
        end

        function value = set.m_linear_eq(obj)
            value = numel(obj.beq);
        end

        function value = set.m_nonlinear_ub(obj)
            if isempty(obj.cub)
                value = 0;
            else
                error('The number of nonlinear inequality constraints is unknown.');
            end
        end

        function value = set.m_nonlinear_eq(obj)
            if isempty(obj.ceq)
                value = 0;
            else
                error('The number of nonlinear equality constraints is unknown.');
            end
        end


        % Getter functions.

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