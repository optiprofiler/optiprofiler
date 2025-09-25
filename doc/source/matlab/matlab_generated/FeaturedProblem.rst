.. _matfeaturedproblem:

FeaturedProblem
===============

**FeaturedProblem(**\ *P*, *F*, *seed*, *max_eval*, *termination_eval*\ **)**
    **FeaturedProblem** is a subclass of **Problem** class and defines an optimization problem with a specific feature.

-------------------------------------------------------------------------

**Problem** and its subclass **FeaturedProblem** describe the following optimization problem:

.. parsed-literal::

    min  **fun**\(**x**)
    s.t. **xl** <= **x** <= **xu**,
         **aub** * **x** <= **bub**,
         **aeq** * **x** = **beq**,
         **cub**\(**x**) <= 0,
         **ceq**\(**x**) = 0,
    with initial point **x0**,

**FeaturedProblem** should be initialized by the following signature:

**FP** = **FeaturedProblem**\(**P**, **F**, **seed**, **max_eval**, **termination_eval**) creates an instance of **FeaturedProblem**, where the input **P** is an instance of **Problem**, the input **F** is an instance of **Feature**, the input **seed** is a nonnegative integer seed less than 2^32, and the inputs **max_eval** and **termination_eval** are positive integers.

The output **FP** contains the following properties:

    - **problem**: the original optimization problem.

    - **feature**: the feature applied to the optimization problem.

    - **max_eval**: the maximum number of function evaluations.

    - **seed**: the seed for the random number generator.

    - **fun_hist**: the history of the evaluated objective function values.

    - **cub_hist**: the history of the evaluated nonlinear inequality constraints.

    - **ceq_hist**: the history of the evaluated nonlinear equality constraints.

    - **maxcv_hist**: the history of the maximum constraint violation.

    - **n_eval_fun**: the minimum between the number of objective function evaluations and **max_eval**.

    - **n_eval_cub**: the minimum between the number of nonlinear inequality constraint evaluations and **max_eval**.

    - **n_eval_ceq**: the minimum between the number of nonlinear equality constraint evaluations and **max_eval**.

    - **fun_init**: the objective function value at the initial point.

    - **maxcv_init**: the maximum constraint violation at the initial point.

The output **FP** contains all the methods of **Problem**, but the methods **fun**, **cub**, **ceq**, and **maxcv** are modified by the input Feature.

Note the following two points.
1. When the number of function evaluations reaches the input **max_eval**, the methods **fun**, **cub**, and **ceq** will return the values of the objective function and constraints at the point where the maximum number of function evaluations is reached, respectively.
2. When the number of function evaluations reaches the input **termination_eval**, the methods **fun**, **cub**, and **ceq** will raise an error to terminate the optimization process.