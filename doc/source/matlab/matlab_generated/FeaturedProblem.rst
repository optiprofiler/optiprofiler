.. _matfeaturedproblem:

FeaturedProblem
===============

**FeaturedProblem(**\ *P*, *F*, *max_eval*, *seed*\ **)**
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

**FP** = **FeaturedProblem**\(**P**, **F**, **max_eval**, **seed**) creates an instance of **FeaturedProblem**, where the input **P** is an instance of **Problem**, the input **F** is an instance of **Feature**, the input **max_eval** is a positive integer that specifies the maximum number of function evaluations allowed, and the input **seed** is a nonnegative integer seed less than 2^32.

The output **FP** contains the following properties:

    - **problem**: the original optimization problem.

    - **feature**: the feature applied to the optimization problem.

    - **max_eval**: the maximum number of function evaluations.

    - **seed**: the seed for the random number generator.

    - **fun_hist**: the history of the evaluated objective function values.

    - **cub_hist**: the history of the evaluated nonlinear inequality constraints.

    - **ceq_hist**: the history of the evaluated nonlinear equality constraints.

    - **maxcv_hist**: the history of the maximum constraint violation.

    - **n_eval**: the number of objective function evaluations.

The output **FP** contains all the methods of **Problem**, but the methods **fun**, **cub**, **ceq**, and **maxcv** are modified by the input Feature.

Pay attention that when the number of function evaluations reaches the input **max_eval**, the methods **fun**, **cub**, and **ceq** will return the values of the objective function and constraints at the point where the maximum number of function evaluations is reached, respectively.