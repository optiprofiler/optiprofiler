.. _matproblem:

Problem
=======

**Problem(**\ *p_struct*\ **)**
    **Problem** is a class that defines an optimization problem.

-------------------------------------------------------------------------

**Problem** describes an optimization problem with the following structure:

.. math::

    \min \quad & \mathrm{fun}(x) \\
    \text{s.t.} \quad & x_l \le x \le x_u, \\
    & A_{\mathrm{ub}} x \le b_{\mathrm{ub}}, \\
    & A_{\mathrm{eq}} x = b_{\mathrm{eq}}, \\
    & c_{\mathrm{ub}}(x) \le 0, \\
    & c_{\mathrm{eq}}(x) = 0, \\
    & \text{with initial point } x_0,

where ``fun`` is the objective function, ``x`` is the variable to optimize, ``xl`` and ``xu`` are the lower and upper bounds, ``aub`` and ``bub`` are the coefficient matrix and right-hand side vector of the linear inequality constraints, ``aeq`` and ``beq`` are the coefficient matrix and right-hand side vector of the linear equality constraints, ``cub`` is the function of nonlinear inequality constraints, ``ceq`` is the function of nonlinear equality constraints.

-------------------------------------------------------------------------

**Problem** should be initialized by the following signature:

**P** = **Problem**\(**p_struct**) creates an instance of the class **Problem** and the input **p_struct** is a struct.

The input struct **p_struct** should contain the following fields:

Compulsory fields:

    - **fun**: The objective function to be minimized. It should accept a vector and return a real number: ``fun(x) -> float``, where ``x`` is a vector.

    - **x0**: the initial guess.

Optional fields:

    - **name**: the name of the problem. It should be a string or a char. Default is ``'Unnamed Problem'``.
    
    - **xl**: the lower bounds on the variable ``x`` in the form of ``xl <= x``. It should be a vector of the same size as ``x``. Default is ``-Inf``.

    - **xu**: the upper bounds on the variable ``x`` in the form of ``xl >= x``. It should be a vector of the same size as ``x``. Default is ``Inf``.

    - **aub**, **bub**: the coefficient matrix and right-hand side vector of the linear inequality constraints ``aub * x <= bub``. The default setting of ``aub`` and ``bub`` are empty matrix and vector.

    - **aeq**, **beq**: the coefficient matrix and right-hand side vector of the linear equality constraints ``aeq * x = beq``. The default setting of `aeq` and `beq` are empty matrix and vector.

    - **cub**: the function of nonlinearly inequality constraints ``cub(x) <= 0``, where ``cub(x) -> float vector``. By default, ``cub(x)`` will return an empty vector.

    - **ceq**: the function of nonlinear equality constraints ``ceq(x) = 0``, where ``ceq(x) -> float vector``. By default, ``ceq(x)`` will return an empty vector.

    - **grad**: the gradient of the objective function ``grad(x) -> float vector``. By default, ``grad(x)`` will return an empty vector.

    - **hess**: the Hessian of the objective function ``hess(x) -> float matrix``. By default, ``hess(x)`` will return an empty matrix.

    - **jcub**: the Jacobian of the nonlinearly inequality constraints ``jcub(x) -> float matrix``. Note that the column size of ``jcub(x)`` should be the same as the length of ``x`` while the row size should be the same as the length of ``cub(x)``. By default, ``jcub(x)`` will return an empty matrix.

    - **jceq**: the Jacobian of the nonlinearly equality constraints ``jceq(x) -> float matrix``. Note that the column size of ``jceq(x)`` should be the same as the length of ``x`` while the row size should be the same as the length of ``ceq(x)``. By default, ``jceq(x)`` will return an empty matrix.

    - **hcub**: the Hessian of the nonlinearly inequality constraints ``hcub(x) -> cell array of float matrices``. The i-th element of ``hcub(x)`` should be the Hessian of the i-th function in ``cub``. By default, ``hcub(x)`` will return an empty cell.

    - **hceq**: the Hessian of the nonlinearly equality constraints ``hceq(x) -> cell array of float matrices``. The i-th element of ``hceq(x)`` should be the Hessian of the i-th function in ``ceq``. By default, ``hceq(x)`` will return an empty cell.

    - **reference**: an optional author/provider fact with exactly five fields: ``fun`` (finite objective), ``maxcv`` (finite nonnegative total violation), ``kind``, ``source`` (non-empty provenance text), and ``mapping``. Kinds are ``'lower_bound'``, ``'optimum'``, ``'best_known'`` and ``'target'``. The first two require ``maxcv=0``; the latter two may describe an infeasible point or a chosen target. The canonical mapping token is ``'objective_violation/1'``. The earlier ``'feasible_objective/1'`` draft is accepted only as a zero-violation objective, including its four-field ``merit`` layout; newly saved records always use the five raw fields. Points and mapping callbacks are not stored. Validation evaluates no objective or constraint.

      Standard profiles convert a retained featured reference using the current ``merit_fun(fun, maxcv, maxcv_init)`` separately for each run, including on load. This replaces the empirical scoring target, not the raw histories. Omitted or safely dropped reference keeps the V1 minimum over the selected solver histories, initial values and optional plain run. This default is not a ``best_known`` claim or a certified optimum. Custom merit functions need not preserve objective lower-bound ordering: ``kind`` describes the raw fact, not a theorem about arbitrary callbacks. Actual run merits are never clamped to the target. See :ref:`FeaturedProblem <matfeaturedproblem>` for feature propagation.

The output **P** contains following properties:

Properties inherited from the input struct:

    **name**, **x0**, **xl**, **xu**, **aub**, **bub**, **aeq**, **beq**, **reference**

    (``reference`` is ``[]`` when unknown, otherwise a struct with ``fun``, ``maxcv``, ``kind``, ``source``, ``mapping``. It is set at construction only and cannot be assigned. An invalid stored object record reads as ``[]`` with a warning; supported legacy zero-violation records are migrated explicitly.)

Properties dependent on the input struct:

    - **ptype**: type of the problem. It should be ``'u'`` (unconstrained), ``'b'`` (bound-constrained), ``'l'`` (linearly constrained), or ``'n'`` (nonlinearly constrained).

    - **n**: dimension of the problem, which is the length of the variable ``x``.

    - **mb**: number of the bound constraints, which is the length of finite elements in ``xl`` and ``xu``.

    - **m_linear_ub**: number of the linear inequality constraints, which is the length of ``bub``.

    - **m_linear_eq**: number of the linear equality constraints, which is the length of ``beq``.

    - **m_nonlinear_ub**: number of the nonlinear inequality constraints, which is the length of ``cub(x)``.

    - **m_nonlinear_eq**: number of the nonlinear equality constraints, which is the length of ``ceq(x)``.

    - **mlcon**: number of the linear constraints, which is the sum of ``m_linear_ub`` and ``m_linear_eq``.

    - **mnlcon**: number of the nonlinear constraints, which is the sum of ``m_nonlinear_ub`` and ``m_nonlinear_eq``.

    - **mcon**: number of the constraints, which is the sum of ``mlcon`` and ``mnlcon``.

The output **P** contains following methods:

Methods inherited from the input struct:

    **fun**, **grad**, **hess**, **cub**, **ceq**, **jcub**, **jceq**, **hcub**, **hceq**

Other methods:

    - **maxcv**: the maximum constraint violation ``maxcv(x) -> float``, which is defined as the maximum of the infinity norms of ``max(xl - x, 0)``, ``max(x - xu, 0)``, ``max(aub * x - bub, 0)``, ``aeq * x - beq``, ``max(cub(x), 0)``, ``ceq(x)``.

    - **project_x0**: trying to project the initial guess ``x0`` onto the feasible region if it is not feasible (but it may fail).
