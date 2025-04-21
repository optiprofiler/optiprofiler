.. _matproblem:

Problem
=======

**Problem(**\ *p_struct*\ **)**
    **Problem** is a class that defines an optimization problem.

-------------------------------------------------------------------------

**Problem** describes an optimization problem with the following structure:

.. parsed-literal::

    min  **fun**\(**x**)
    s.t. **xl** <= **x** <= **xu**,
         **aub** * **x** <= **bub**,
         **aeq** * **x** = **beq**,
         **cub**\(**x**) <= 0,
         **ceq**\(**x**) = 0,
    with initial point **x0**,

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

    - **aeq**, **beq**: the coefficient matrix and right-hand side vector of the linear equality constraints ``aeq * x <= beq``. The default setting of `aeq` and `beq` are empty matrix and vector.

    - **cub**: the function of nonlinearly inequality constraints ``cub(x) <= 0``, where ``cub(x) -> float vector``. By default, ``cub(x)`` will return an empty vector.

    - **ceq**: the function of nonlinearly equality constraints ``ceq(x) <= 0``, where ``ceq(x) -> float vector``. By default, ``ceq(x)`` will return an empty vector.

    - **grad**: the gradient of the objective function ``grad(x) -> float vector``. By default, ``grad(x)`` will return an empty vector.

    - **hess**: the Hessian of the objective function ``hess(x) -> float matrix``. By default, ``hess(x)`` will return an empty matrix.

    - **jcub**: the Jacobian of the nonlinearly inequality constraints ``jcub(x) -> float matrix``. Note that the column size of ``jcub(x)`` should be the same as the length of ``x`` while the row size should be the same as the length of ``cub(x)``. By default, ``jcub(x)`` will return an empty matrix.

    - **jceq**: the Jacobian of the nonlinearly equality constraints ``jceq(x) -> float matrix``. Note that the column size of ``jceq(x)`` should be the same as the length of ``x`` while the row size should be the same as the length of ``ceq(x)``. By default, ``jceq(x)`` will return an empty matrix.

    - **hcub**: the Hessian of the nonlinearly inequality constraints ``hcub(x) -> cell array of float matrices``. The i-th element of ``hcub(x)`` should be the Hessian of the i-th function in ``cub``. By default, ``hcub(x)`` will return an empty cell.

    - **hceq**: the Hessian of the nonlinearly equality constraints ``hceq(x) -> cell array of float matrices``. The i-th element of ``hceq(x)`` should be the Hessian of the i-th function in ``ceq``. By default, ``hceq(x)`` will return an empty cell.

The output **P** contains following properties:

Properties inherited from the input struct:

    **name**, **x0**, **xl**, **xu**, **aub**, **bub**, **aeq**, **beq**

Properties dependent on the input struct:

    - **n**: dimension of the problem, which is the length of the variable ``x``.

    - **m_linear_ub**: number of the linear inequality constraints, which is the length of ``bub``.

    - **m_linear_eq**: number of the linear equality constraints, which is the length of ``beq``.

    - **m_nonlinear_ub**: number of the nonlinear inequality constraints, which is the length of ``cub(x)``.

    - **m_nonlinear_eq**: number of the nonlinear equality constraints, which is the length of ``ceq(x)``.

    - **ptype**: type of the problem. It should be ``'u'`` (unconstrained), ``'b'`` (bound-constrained), ``'l'`` (linearly constrained), or ``'n'`` (nonlinearly constrained).

The output **P** contains following methods:

Methods inherited from the input struct:

    **fun**, **grad**, **hess**, **cub**, **ceq**, **jcub**, **jceq**, **hcub**, **hceq**

Other methods:

    - **maxcv**: the maximum constraint violation ``maxcv(x) -> float``, which is defined as the maximum of the infinity norms of ``max(xl - x, 0)``, ``max(x - xu, 0)``, ``max(aub * x - bub, 0)``, ``aeq * x - beq``, ``max(cub(x), 0)``, ``ceq(x)``.

    - **project_x0**: trying to project the initial guess ``x0`` onto the feasible region if it is not feasible (but it may fail).