.. _matfeature:

Feature
=======

**Feature(**\ *name*, *options*\ **)**
    **Feature** is a class that defines a mapping from an optimization problem to a new one with specified features.

-------------------------------------------------------------------------

We are interested to test solvers on problems with different features. For example, we want to test the performance of solvers under the case where the objective function is noisy. For this purpose, we define **Feature** class.

Suppose we have an optimization problem

.. parsed-literal::

    min  **fun**\(**x**)
    s.t. **xl** <= **x** <= **xu**,
         **aub** * **x** <= **bub**,
         **aeq** * **x** = **beq**,
         **cub**\(**x**) <= 0,
         **ceq**\(**x**) = 0,
    with initial point **x0**。

Then, **Feature** maps the above optimization problem to the following one.

.. parsed-literal::

    min  **fun_mod**\(**A** * **x** + **b**)
    s.t. **xl_mod** <= **A** * **x** + **b** <= **xu_mod**,
         **aub_mod** * (**A** * **x** + **b**) <= **bub_mod**,
         **aeq_mod** * (**A** * **x** + **b**) = **beq_mod**,
         **cub_mod**\(**A** * **x** + **b**) <= 0,
         **ceq_mod**\(**A** * **x** + **b**) = 0,
    with initial guess **x0_mod**.

-------------------------------------------------------------------------------

**Feature** should be initialized by the following signature:

**F** = **Feature**\(**name**) creates an instance of **Feature** with the feature name **name**;

**F** = **Feature**\(**name**, **options**) creates an instance of **Feature** with the feature name **name** and a struct **options**.

The input **name** should be one of the following char or string:

    1. **\'plain\'**:
        do nothing to the optimization problem.

    2. **\'perturbed_x0\'**:
        perturb the initial guess ``x0``.

    3. **\'noisy\'**:
        add noise to the objective function and nonlinear constraints.

    4. **\'truncated\'**:
        truncate values of the objective function and nonlinear constraints to a given number of significant digits.

    5. **\'permuted\'**:
        randomly permute the variables. Note that ``x0``, ``xl``, ``xu``, ``aub``, ``bub``, ``aeq``, ``beq`` will be modified since we want to keep the new problem mathematically equivalent to the original one.

    6. **\'linearly_transformed\'**:
        generate a invertible linear transformation in the form D * Q' with D being a diagonal matrix and Q being a random orthogonal matrix, and apply the transformation to the variables. In this way, the Hessian of the objective function becomes Q * D * H * D * Q' with H being the original Hessian. Note that ``x0``, ``xl``, ``xu``, ``aub``, ``bub``, ``aeq``, ``beq`` will be modified since we want to keep the new problem mathematically equivalent to the original one.

    7. **\'random_nan\'**:
        randomly replace values of the objective function and nonlinear constraints with ``NaN``.

    8. **\'unrelaxable_constraints\'**:
        set the objective function to ``Inf`` outside the feasible region.

    9. **\'nonquantifiable_constraints\'**:
        replace values of nonlinear constraints with either ``0`` (if the constraint is satisfied) or ``1`` (if the constraint is violated).

    10. **\'quantized\'**:
            quantize the objective function and nonlinear constraints.

    11. **\'custom\'**:
            user-defined **Feature**.

The optional input **options** should be a struct which can contain the following fields:

  - **n_runs**: the number of runs of the experiments under the given feature. Default is ``5`` for stochastic features and ``1`` for deterministic features.

  - **distribution**: the distribution of perturbation in ``'perturbed_x0'`` feature or noise in ``'noisy'`` feature. It should be either a string (or char), or a function handle

      ``(random_stream, dimension) -> random vector``

    that accepts a ``random_stream`` and the ``dimension`` of a problem and returning a ``random vector`` with the given ``dimension``. In ``'perturbed_x0'`` case, the char should be either ``'spherical'`` or ``'gaussian'`` (default is ``'spherical'``). In ``'noisy'`` case, the char should be either ``'gaussian'`` or ``'uniform'`` (default is ``'gaussian'``).

  - **perturbation_level**: the magnitude of the perturbation to the initial guess in the ``'perturbed_x0'`` feature. Default is ``1e-3``.

  - **noise_level**: the magnitude of the noise in the ``'noisy'`` feature. Default is ``1e-3``.

  - **noise_type**: the type of the noise in the ``'noisy'`` features. It should be either ``'absolute'``, ``'relative'``, or ``'mixed'``. Default is ``'mixed'``.

  - **significant_digits**: the number of significant digits in the ``'truncated'`` feature. Default is ``6``.

  - **perturbed_trailing_digits**: whether we will randomize the trailing zeros of the objective function value in the ``'perturbed_x0'`` feature. Default is ``false``.

  - **rotated**: whether to use a random or given rotation matrix to rotate the coordinates of a problem in the ``'linearly_transformed'`` feature. Default is ``true``.

  - **condition_factor**: the scaling factor of the condition number of the linear transformation in the ``'linearly_transformed'`` feature. More specifically, the condition number of the linear transformation will be ``2 ^ (condition_factor * n / 2)``, where ``n`` is the dimension of the problem. Default is ``0``.

  - **nan_rate**: the probability that the evaluation of the objective function will return NaN in the ``'random_nan'`` feature. Default is ``0.05``.

  - **unrelaxable_bounds**: whether the bound constraints are unrelaxable or not in the ``'unrelaxable_constraints'`` feature. Default is ``true``.

  - **unrelaxable_linear_constraints**: whether the linear constraints are unrelaxable or not in the ``'unrelaxable_constraints'`` feature. Default is ``false``.

  - **unrelaxable_nonlinear_constraints**: whether the nonlinear constraints are unrelaxable or not in the ``'unrelaxable_constraints'`` feature. Default is ``false``.

  - **mesh_size**: the size of the mesh in the ``'quantized'`` feature. Default is ``1e-3``.

  - **mesh_type**: the type of the mesh in the ``'quantized'`` feature. It should be either ``'absolute'`` or ``'relative'``. Default is ``'absolute'``.

  - **ground_truth**: whether the featured problem is the ground truth or not in the ``'quantized'`` feature. Default is ``true``.

  - **mod_x0**: the modifier function to modify the inital guess in the ``'custom'`` feature. It should be a function handle

      ``(random_stream, problem) -> modified_x0``,

    where ``problem`` is an instance of the class Problem, and ``modified_x0`` is the modified initial guess. No default.

  - **mod_affine**: the modifier function to generate the affine transformation applied to the variables in the ``'custom'`` feature. It should be a function handle

      ``(random_stream, problem) -> (A, b, inv)``,

    where ``problem`` is an instance of the class Problem, ``A`` is the matrix of the affine transformation, ``b`` is the vector of the affine transformation, and ``inv`` is the inverse of matrix ``A``. No default.

  - **mod_bounds**: the modifier function to modify the bound constraints in the ``'custom'`` feature. It should be a function handle

      ``(random_stream, problem) -> (modified_xl, modified_xu)``,

    where ``problem`` is an instance of the class Problem, ``modified_xl`` is the modified lower bound, and ``modified_xu`` is the modified upper bound. No default.

  - **mod_linear_ub**: the modifier function to modify the linear inequality constraints in the ``'custom'`` feature. It should be a function handle

      ``(random_stream, problem) -> (modified_aub, modified_bub)``,

    where ``problem`` is an instance of the class Problem, ``modified_aub`` is the modified matrix of the linear inequality constraints, and ``modified_bub`` is the modified vector of the linear inequality constraints. No default.

  - **mod_linear_eq**: the modifier function to modify the linear equality constraints in the ``'custom'`` feature. It should be a function handle

      ``(random_stream, problem) -> (modified_aeq, modified_beq)``,

    where ``problem`` is an instance of the class Problem, ``modified_aeq`` is the modified matrix of the linear equality constraints, and ``modified_beq`` is the modified vector of the linear equality constraints. No default.

  - **mod_fun**: the modifier function to modify the objective function in the ``'custom'`` feature. It should be a function handle

      ``(x, random_stream, problem) -> modified_fun``,

    where ``x`` is the evaluation point, ``problem`` is an instance of the class Problem, and ``modified_fun`` is the modified objective function value. No default.

  - **mod_cub**: the modifier function to modify the nonlinear inequality constraints in the ``'custom'`` feature. It should be a function handle

      ``(x, random_stream, problem) -> modified_cub``,

    where ``x`` is the evaluation point, ``problem`` is an instance of the class Problem, and ``modified_cub`` is the modified vector of the nonlinear inequality constraints. No default.

  - **mod_ceq**: the modifier function to modify the nonlinear equality constraints in the ``'custom'`` feature. It should be a function handle

      ``(x, random_stream, problem) -> modified_ceq``,

    where ``x`` is the evaluation point, ``problem`` is an instance of the class Problem, and ``modified_ceq`` is the modified vector of the nonlinear equality constraints. No default.

-------------------------------------------------------------------------------

Different input **name** will have different valid fields of **options**. We list the valid fields for each input **name** as

    1. **\'plain\'**:
        **n_runs**

    2. **\'perturbed_x0\'**:
        **n_runs**, **distribution**, **perturbation_level**

    3. **\'noisy\'**:
        **n_runs**, **distribution**, **noise_level**, **noise_type**

    4. **\'truncated\'**:
        **n_runs**, **significant_digits**, **perturbed_trailing_digits**

    5. **\'permuted\'**:
        **n_runs**

    6. **\'linearly_transformed\'**:
        **n_runs**, **rotated**, **condition_factor**

    7. **\'random_nan\'**:
        **n_runs**, **nan_rate**

    8. **\'unrelaxable_constraints\'**:
        **n_runs**, **unrelaxable_bounds**, **unrelaxable_linear_constraints**, **unrelaxable_nonlinear_constraints**

    9. **\'nonquantifiable_constraints\'**:
        **n_runs**

    10. **\'quantized\'**:
            **n_runs**, **mesh_size**, **mesh_type**, **ground_truth**

    11. **\'custom\'**:
            **n_runs**, **mod_x0**, **mod_affine**, **mod_bounds**, **mod_linear_ub**, **mod_linear_eq**, **mod_fun**, **mod_cub**, **mod_ceq**

-------------------------------------------------------------------------------

The output **F** contains the following methods:

    1. **modifier_x0**:
        a function handle to modify the initial guess.

    2. **modifier_affine**:
        a function handle to generate an invertible matrix ``A`` (and its inverse) and a vector ``b`` for the affine transformation applied to the variables.

    3. **modifier_bounds**:
        a function handle to modify the lower and upper bounds.
    
    4. **modifier_linear_ub**:
        a function handle to modify the linear inequality constraints.

    5. **modifier_linear_eq**:
        a function handle to modify the linear equality constraints.

    6. **modifier_fun**:
        a function handle to modify the objective function value.

    7. **modifier_cub**:
        a function handle to modify the values of the nonlinear inequality constraints.

    8. **modifier_ceq**:
        a function handle to modify the values of the nonlinear equality constraints.

All the methods of **F** will be used later to modify the optimization
problem.