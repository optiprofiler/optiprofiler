.. _use_python:

Usage for Python
================

OptiProfiler provides a :func:`~optiprofiler.benchmark` function. This is the main entry point to the package. It benchmarks given solvers on the selected test suite.

We provide below simple examples on how to use OptiProfiler in Python. For more details on the signature of the :func:`~optiprofiler.benchmark` function, please refer to the :ref:`Python API documentation <pythonapi>`.

Examples
--------

.. _py_example1:

Example 1: first example to try out
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

Let us first try to benchmark two callable optimization solvers **solver1** and **solver2** on the default test suite.
(Note that each **solver** must accept signatures mentioned in the `Cautions` part of the :func:`~optiprofiler.benchmark` function according to the type of problems you want to solve.)

To do this, run:

.. code-block:: python

    from optiprofiler import benchmark

    scores = benchmark([solver1, solver2])

This will benchmark the two solvers under the default test setting, which means ``'plain'`` feature (see :class:`~optiprofiler.Feature`) and unconstrained problems from the default problem library whose dimension is smaller or equal to 2. It will also return the scores of the two solvers based on the profiles.

There will be a new folder named ``out`` in the current working directory, which contains a subfolder named ``plain_<timestamp>`` with all the detailed results. Additionally, a PDF file named ``summary.pdf`` is generated, summarizing all the performance profiles and data profiles.

.. _py_example2:

Example 2: one step further by adding options
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

You can also add options to the benchmark function. For example, if you want to benchmark three solvers **solver1**, **solver2**, and **solver3** on the test suite with the ``'noisy'`` feature and all the unconstrained and bound-constrained problems with dimension between 6 and 10 from the default problem set, you can run:

.. code-block:: python

    from optiprofiler import benchmark

    scores = benchmark(
        [solver1, solver2, solver3],
        ptype='ub',
        mindim=6,
        maxdim=10,
        feature_name='noisy',
    )

This will create the corresponding folders ``out/noisy_<timestamp>`` and files as in :ref:`Example 1 <py_example1>`. More details on the options can be found in the :func:`~optiprofiler.benchmark` function documentation.

.. _py_example3:

Example 3: useful option **load**
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

OptiProfiler provides a practically useful option named ``load``. This option allows you to load the results from a previous benchmarking run (without solving all the problems again) and use them to draw new profiles with different options. For example, if you have just run :ref:`Example 2 <py_example2>` and OptiProfiler has finished the job and successfully created the folder ``out`` in the current working directory, you can run:

.. code-block:: python

    from optiprofiler import benchmark

    scores = benchmark(
        load='latest',
        solvers_to_load=[0, 2],
        ptype='u',
        mindim=7,
        maxdim=9,
    )

This will directly draw the profiles for the **solver1** and **solver3** with the ``'noisy'`` feature and all the unconstrained problems with dimension between 7 and 9 selected from the previous run. The results will also be saved under the current directory with a new subfolder named ``noisy_<timestamp>`` with the new timestamp.

.. _py_example4:

Example 4: testing parametrized solvers
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

If you want to benchmark a solver with one variable parameter, you can define callables by looping over the parameter values. For example, if **solver** accepts the signature ``solver(fun, x0, para)``, and you want to benchmark it with the parameter ``para`` taking values from 1 to 3, you can run:

.. code-block:: python

    from optiprofiler import benchmark

    def make_solver(para):
        def solver_wrapper(fun, x0):
            return solver(fun, x0, para)
        return solver_wrapper

    solvers = [make_solver(i) for i in range(1, 4)]
    solver_names = [f'solver{i}' for i in range(1, 4)]
    scores = benchmark(solvers, solver_names=solver_names)

.. note::

    We use named functions (``def``) instead of lambda expressions here.
    Lambda functions cannot be serialized by Python's ``multiprocessing``
    module, which would force OptiProfiler to fall back to sequential
    execution. Using named functions ensures that the benchmark can run in
    parallel when ``n_jobs > 1``.

.. _py_example5:

Example 5: customizing the test suite
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

OptiProfiler allows you to customize the test suite by creating your own feature and loading your own problem library.
For example, if you want to create a new feature that adds noise to the objective function and perturbs the initial guess at the same time, you can try the following:

.. code-block:: python

    from optiprofiler import benchmark

    def mod_fun(x, rand_stream, problem):
        return problem.fun(x) + 1e-3 * rand_stream.standard_normal()

    def mod_x0(rand_stream, problem):
        return problem.x0 + 1e-3 * rand_stream.standard_normal(problem.n)

    scores = benchmark(
        [solver1, solver2],
        feature_name='custom',
        mod_fun=mod_fun,
        mod_x0=mod_x0,
    )

.. note::

    In this example, we use named functions (``def``) instead of lambda
    expressions. This is recommended because lambda functions cannot be
    serialized by Python's ``multiprocessing`` module, which would prevent
    parallel execution. Using named functions ensures that the benchmark can
    run in parallel when ``n_jobs > 1``.

If you want to benchmark solvers based on your own problem library, you should do the following three steps:

1. Create a directory anywhere on your system (e.g., ``'/path/to/my_libs'``), and create a subfolder inside it for your problem library (e.g., ``'myproblems'``), so the structure looks like::

       /path/to/my_libs/
       └── myproblems/
           └── myproblems_tools.py

2. In ``myproblems_tools.py``, implement two functions:

   - **myproblems_load**: A function that accepts a string representing the optimization problem name and returns a :class:`~optiprofiler.Problem` instance.

   - **myproblems_select**: A function that accepts a dictionary to specify desired problem characteristics and returns a list of problem names that satisfy the requirements.

   In general, the module should be named ``<library_name>_tools.py`` and the two functions should be named ``<library_name>_load`` and ``<library_name>_select``.

3. Use the benchmark function with the ``custom_problem_libs_path`` option pointing to your directory. For example, to use both the default S2MPJ library and your custom library ``'myproblems'``, you can run:

.. code-block:: python

    scores = benchmark(
        [solver1, solver2],
        plibs=['s2mpj', 'myproblems'],
        custom_problem_libs_path='/path/to/my_libs',
    )

For a detailed guide on the required structure of a custom problem library, please refer to the `guide on our website <https://www.optprof.com>`_ or the `README on GitHub <https://github.com/optiprofiler/optiprofiler/blob/main/python/optiprofiler/problem_libs/README.txt>`_.
