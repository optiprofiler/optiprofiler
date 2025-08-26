.. _matbenchmark:

benchmark
=========

**benchmark(**\ *varargin*\ **)**
    **benchmark** creates multiple profiles for benchmarking optimization solvers on a set of problems with different features.

-------------------------------------------------------------------------

*Signatures*:

**solver_scores** = **benchmark**\(**solvers**) creates performance profiles, data profiles, and log-ratio profiles for the given **solvers** on the default unconstrained problem set, returning **solver_scores** based on the profiles. **solvers** is a cell array of function handles. We require **solvers** to accept specified inputs and return specified outputs. Details can be found in the following *Cautions* part.

**solver_scores** = **benchmark**\(**solvers**, **feature_name**) creates profiles for the given **solvers** on the default unconstrained problem set with the specified feature **feature_name**.

**solver_scores** = **benchmark**\(**solvers**, **options**) creates profiles for the given **solvers** with options specified in the struct **options**. See *Options* part for more details.

**solver_scores** = **benchmark**\(**options**) creates profiles with options specified in the struct **options**. Note that the struct **options** should at least contain the field **load** with the value **\'latest\'** or a time stamp of an experiment in the format of 'yyyyMMdd_HHmmss'. In this case, we will load the data from the specified experiment and draw the profiles.

[**solver_scores**, **profile_scores**] = **benchmark**\(...) returns a 4D tensor **profile_scores** containing scores for all profiles. See **score_fun** in *Options* part for more details.

[**solver_scores**, **profile_scores**, **curves**] = **benchmark**\(...) returns a cell array **curves** containing the curves of all the profiles.

-----------------------------------------------------------------------------

*Cautions*:

1. Each **solver** in **solvers** should accept corresponding signature(s) depending on the test suite you choose:

- For an unconstrained problem:

  .. code-block:: matlab
  
    x = solver(fun, x0)
        
  where ``fun`` is a function handle of the objective function accepting a column vector and returning a real number, and ``x0`` is the initial guess which is a column vector;

- For a bound-constrained problem:

  .. code-block:: matlab
  
    x = solver(fun, x0, xl, xu)
        
  where ``xl`` and ``xu`` are the lower and upper bounds of the variables which are column vectors (they can contain Inf or -Inf);

- For a linearly constrained problem:

  .. code-block:: matlab
  
    x = solver(fun, x0, xl, xu, aub, bub, aeq, beq)
        
  where ``aub`` and ``aeq`` are the matrices of the linear inequality and equality constraints, and ``bub`` and ``beq`` are the vectors of the linear inequality and equality constraints;

- For a nonlinearly constrained problem:

  .. code-block:: matlab
  
    x = solver(fun, x0, xl, xu, aub, bub, aeq, beq, cub, ceq)
        
  where ``cub`` and ``ceq`` are the functions of the nonlinear inequality and equality constraints accepting a column vector and returning a column vector.

2. The log-ratio profiles are available only when there are exactly two solvers. For more information of performance and data profiles, see [1]_, [2]_, [5]_. For that of log-ratio profiles, see [4]_, [6]_.

-------------------------------------------------------------------------

*Options*:

Options should be specified in a struct. The following are the available fields of the struct:

1. *Options for profiles and plots*:

  - **bar_colors**: two different colors for the bars of two solvers in the log-ratio profiles. It can be a cell array of short names of colors ``('r', 'g', 'b', 'c', 'm', 'y', 'k')`` or a 2-by-3 matrix with each row being a RGB triplet. Default is set to the first two colors in the ``line_colors`` option.

  - **benchmark_id**: the identifier of the test. It is used to create the specific directory to store the results. Default is ``'out'`` if the option ``load`` is not provided, otherwise default is ``'.'``.

  - **errorbar_type**: the type of the uncertainty interval that can be either ``'minmax'`` or ``'meanstd'``. When ``n_runs`` is greater than 1, we run several times of the experiments and get average curves and get average curves and uncertainty intervals. Default is ``'minmax'``, meaning that we takes the pointwise minimum and maximum of the curves.

  - **feature_stamp**: the stamp of the feature with the given options. It is used to create the specific directory to store the results. Default depends on features.

  - **line_colors**: the colors of the lines in the plots. It can be a cell array of short names of colors ``('r', 'g', 'b', 'c', 'm', 'y', 'k')`` or a matrix with each row being a RGB triplet. Default line colors are those in the palettename named "gem" (see MATLAB documentation for 'colororder'). Note that if the number of solvers is greater than the number of colors, we will cycle through the colors.

  - **line_styles**: the styles of the lines in the plots. It can be a cell array of chars that are the combinations of line styles ``('-', '-.', '--', ':')`` and markers ``('o', '+', '*', '.', 'x', 's', 'd', '^', 'v', '>', '<', 'p', 'h')``. Default line style order is ``{'-', '-.', '--', ':'}``. Note that if the number of solvers is greater than the number of line styles, we will cycle through the styles.

  - **line_widths**: the widths of the lines in the plots. It should be a positive scalar or a vector. Default is ``1.5``. Note that if the number of solvers is greater than the number of line widths, we will cycle through the widths.

  - **load**: loading the stored data from a completed experiment and draw profiles. It can be either ``'latest'`` or a time stamp of an experiment in the format of 'yyyyMMdd_HHmmss'. No default.

  - **max_eval_factor**: the factor multiplied to each problem's dimension to get the maximum number of evaluations for each problem. Default is ``500``.

  - **max_tol_order**: the maximum order of the tolerance. In any profile (performance profiles, data profiles, and log-ratio profiles), we need to set a group of 'tolerances' to define the convergence test of the solvers. (Details can be found in the references.) We will set the tolerances as ``10^(-1:-1:-max_tol_order)``. Default is ``10``.

  - **merit_fun**: the merit function to measure the quality of a point using the objective function value and the maximum constraint violation. It should be a function handle

      ``(fun_value, maxcv_value, maxcv_init) -> merit_value``,

    where ``fun_value`` is the objective function value, ``maxcv_value`` is the maximum constraint violation, and ``maxcv_init`` is the maximum constraint violation at the initial guess. The size of ``fun_values`` and ``maxcv_values`` is the same, and the size of ``maxcv_init`` is the same as the second to last dimensions of ``fun_values``. The default merit function ``varphi(x)`` is defined by the objective function ``f(x)`` and the maximum constraint violation ``v(x)`` as

    .. parsed-literal::

        **varphi**\(**x**) = **f**\(**x**),                      if **v**\(**x**) <= v1,
        **varphi**\(**x**) = **f**\(**x**) + 1e5 * (**v**\(**x**) - v1),  if v1 < **v**\(**x**) <= v2,
        **varphi**\(**x**) = Inf,                       if **v**\(**x**) > v2,

    where ``v1 = max(1e-5, v0)`` and ``v2 = min(0.01, 1e-10 * max(1, v0))``, and ``v0`` is the initial maximum constraint violation.

  - **n_jobs**: the number of parallel jobs to run the test. Default is the default number of workers in the default local cluster.
  
  - **normalized_scores**: whether to normalize the scores of the solvers by the maximum score of the solvers. Default is ``false``.

  - **project_x0**: whether to project the initial point to the feasible set. Default is ``false``.

  - **run_plain**: whether to run an extra experiment with the ``'plain'`` feature. Default is ``false``.

  - **savepath**: the path to store the results. Default is ``'pwd'``, the current working directory.

  - **score_fun**: the scoring function to calculate the scores of the solvers. It should be a function handle

      ``profile_scores -> solver_scores``,

    where ``profile_scores`` is a 4D tensor containing scores for all profiles. The first dimension of ``profile_scores`` corresponds to the index of the solver, the second corresponds to the index of tolerance starting from 1, the third represents history-based or output-based profiles, and the fourth represents performance profiles, data profiles, or log-ratio profiles. The default scoring function takes the average of the history-based performance profiles under all the tolerances.

  - **score_only**: whether to only calculate the scores of the solvers without drawing the profiles and saving the data. Default is ``false``.

  - **score_weight_fun**: the weight function to calculate the scores of the solvers in the performance and data profiles. It should be a function handle representing a nonnegative function in R^+. Default is ``1``.

  - **seed**: the seed of the random number generator. Default is ``0``.

  - **semilogx**: whether to use the semilogx scale during plotting profiles (performance profiles and data profiles). Default is ``true``.

  - **silent**: whether to show the information of the progress. Default is ``false``.

  - **solver_isrand**: whether the solvers are randomized or not. It is a logical array of the same length as the number of solvers, where the value is true if the solver is randomized, and false otherwise. Default is all false. Note that if ``n_runs`` is not specified, we will set it 5 for the randomized solvers.

  - **solver_names**: the names of the solvers. Default is the names of the function handles in **solvers**.

  - **solver_verbose**: the level of the verbosity of the solvers. ``0`` means no verbosity, ``1`` means some verbosity, and ``2`` means full verbosity. Default is ``1``.

  - **solvers_to_load**: the indices of the solvers to load when the ``load`` option is provided. It can be a vector of different integers selected from 1 to the total number of solvers of the loading experiment. At least two indices should be provided. Default is all the solvers.

  - **summarize_data_profiles**: whether to add all the data profiles to the summary PDF. Default is ``true``.

  - **summarize_log_ratio_profiles**: whether to add all the log-ratio profiles to the summary PDF. Default is ``false``.

  - **summarize_output_based_profiles**: whether to add all the output-based profiles of the selected profiles to the summary PDF. Default is ``true``.

  - **summarize_performance_profiles**: whether to add all the performance profiles to the summary PDF. Default is ``true``.

  - **xlabel_data_profile**: the label of the x-axis of the data profiles. Default is ``'Number of simplex gradients'``. Note: the ``'Interpreter'`` property is set to ``'latex'``, so LaTeX formatting is supported. The same applies to the options ``xlabel_log_ratio_profile``, ``xlabel_performance_profile``, ``ylabel_data_profile``, ``ylabel_log_ratio_profile``, and ``ylabel_performance_profile``.

  - **xlabel_log_ratio_profile**: the label of the x-axis of the log-ratio profiles. Default is ``'Problem'``.

  - **xlabel_performance_profile**: the label of the x-axis of the performance profiles. Default is ``'Performance ratio'``.

  - **ylabel_data_profile**: the label of the y-axis of the data profiles. Default is ``'Data profiles ($\\mathrm{tol} = %s$)'``, where ``%s`` will be replaced by the current tolerance in LaTeX format. You can also use ``%s`` in your custom label, and it will be replaced accordingly. The same applies to the options ``ylabel_log_ratio_profile`` and ``ylabel_performance_profile``.

  - **ylabel_log_ratio_profile**: the label of the y-axis of the log-ratio profiles. Default is ``'Log-ratio profiles ($\\mathrm{tol} = %s$)'``, where ``%s`` will be replaced by the current tolerance in LaTeX format.

  - **ylabel_performance_profile**: ylabel_performance_profile: the label of the y-axis of the performance profiles. Default is ``'Performance profiles ($\\mathrm{tol} = %s$)'``, where ``%s`` will be replaced by the current tolerance in LaTeX format.

2. *Options for features*:

  - **feature_name**: the name of the feature. The available features are ``'plain'``, ``'perturbed_x0'``, ``'noisy'``, ``'truncated'``, ``'permuted'``, ``'linearly_transformed'``, ``'random_nan'``, ``'unrelaxable_constraints'``, ``'nonquantifiable_constraints'``, ``'quantized'``, and ``'custom'``. Default is ``'plain'``.

  - **n_runs**: the number of runs of the experiments with the given feature. Default is ``5`` for stochastic features and ``1`` for deterministic features.

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

3. *Options for problems*:

Options in this part are used to select problems for benchmarking. First select which problem libraries to use based on the ``plibs`` option. Then select problems from these libraries according to the given options (``problem_names``, ``ptype``, ``mindim``, ``maxdim``, ``minb``, ``maxb``, ``minlcon``, ``maxlcon``, ``minnlcon``, ``maxnlcon``, ``mincon``, ``maxcon``, and ``excludelist``). Following is the list of available options:

  - **plibs**: the problem libraries to be used. It should be a cell array of strings or chars. The available choices are subfolder names in the ``problems`` directory. There are three subfolders after installing the package: ``s2mpj``, ``matcutest``, and ``custom``. Default setting is ``'s2mpj'``.

  - **ptype**: the type of the problems to be selected. It should be a string or char consisting of any combination of ``'u'`` (unconstrained), ``'b'`` (bound constrained), ``'l'`` (linearly constrained), and ``'n'`` (nonlinearly constrained), such as ``'b'``, ``'ul'``, ``'ubn'``. Default is ``'u'``.

  - **mindim**: the minimum dimension of the problems to be selected. Default is ``1``.

  - **maxdim**: the maximum dimension of the problems to be selected. Default is ``mindim + 10``.

  - **minb**: the minimum number of bound constraints of the problems to be selected. Default is ``0``.

  - **maxb**: the maximum number of bound constraints of the problems to be selected. Default is ``minb + 10``.

  - **minlcon**: the minimum number of linear constraints of the problems to be selected. Default is ``0``.

  - **maxlcon**: the maximum number of linear constraints of the problems to be selected. Default is ``minlcon + 10``.

  - **minnlcon**: the minimum number of nonlinear constraints of the problems to be selected. Default is ``0``.

  - **maxnlcon**: the maximum number of nonlinear constraints of the problems to be selected. Default is ``minnlcon + 10``.

  - **mincon**: the minimum number of linear and nonlinear constraints of the problems to be selected. Default is ``min(minlcon, minnlcon)``.

  - **maxcon**: the maximum number of linear and nonlinear constraints of the problems to be selected. Default is ``max(maxlcon, maxnlcon)``.

  - **excludelist**: the list of problems to be excluded. Default is not to exclude any problem.

  - **problem_names**: the names of the problems to be selected. It should be a cell array of strings or chars. Default is not to select any problem by name but by the options above.

You may also pass an instance of the class Problem by the option

  - **problem**: a problem to be benchmarked. It should be an instance of the class Problem. If it is provided, we will only run the test on this problem with the given feature and draw the history plots. Default is not to set any problem.

*Several points to note*:

  1. The information about two problem libraries is available in the following links:
    S2MPJ (see [3]_) <https://github.com/GrattonToint/S2MPJ> and MatCUTEst <https://github.com/matcutest>.

  2. If you want to use your own problem library, please check the README.txt in the directory ``problems`` or the :ref:`guidance <use>` in our website for more details.

  3. The problem library MatCUTEst is only available when the OS is Linux.

  4. If the ``load`` option is provided, we will use the provided options to select data from the specified experiment for plotting the profiles. Available options are:

    - *Options for profiles and plots*: ``benchmark_id``, ``solver_names``, ``feature_stamp``, ``errorbar_type``, ``savepath``, ``max_tol_order``, ``merit_fun``, ``run_plain``, ``score_only``, ``summarize_performance_profiles``, ``summarize_data_profiles``, ``summarize_log_ratio_profiles``, ``summarize_output_based_profiles``, ``silent``, ``semilogx``, ``normalized_scores``, ``score_weight_fun``, ``score_fun``, ``solvers_to_load``, ``line_colors``, ``line_styles``, ``line_widths``, ``bar_colors``.

    - *Options for features*: none.

    - *Options for problems*: ``plibs``, ``ptype``, ``mindim``, ``maxdim``, ``minb``, ``maxb``, ``minlcon``, ``maxlcon``, ``minnlcon``, ``maxnlcon``, ``mincon``, ``maxcon``, ``excludelist``.

-----------------------------------------------------------------------------

*References*:

.. [1] E. D. Dolan and J. J. Moré. Benchmarking optimization software with
    performance profiles. *Math. Program.*, 91(2):201--213, 2002.
    doi:10.1007/s101070100263
    <https://doi.org/10.1007/s101070100263>.

.. [2] N. Gould and J. Scott. A note on performance profiles for
    benchmarking software. *ACM Trans. Math. Software*, 43(2):15:1--5,
    2016. doi:10.1145/2950048 <https://doi.org/10.1145/2950048>.

.. [3] S. Gratton and Ph. L. Toint. S2MPJ and CUTEst optimization problems
    for Matlab, Python and Julia. arXiv:2407.07812, 2024.

.. [4] J. L. Morales. A numerical study of limited memory BFGS methods.
    *Appl. Math. Lett.*, 15(4):481--487, 2002.
    doi:10.1016/S0893-9659(01)00162-8
    <https://doi.org/10.1016/S0893-9659(01)00162-8>.

.. [5] J. J. Moré and S. M. Wild. Benchmarking derivative-free optimization
    algorithms. *SIAM J. Optim.*, 20(1):172--191, 2009.
    doi:10.1137/080724083 <https://doi.org/10.1137/080724083>.

.. [6] H.-J. M. Shi, M. Q. Xuan, F. Oztoprak, and J. Nocedal. On the
    numerical performance of finite-difference-based methods for
    derivative-free optimization. *Optim. Methods Softw.*,
    38(2):289--311, 2023. doi:10.1080/10556788.2022.2121832
    <https://doi.org/10.1080/10556788.2022.2121832>.