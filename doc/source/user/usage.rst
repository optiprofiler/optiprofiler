.. _use:

Usage for MATLAB
================

OptiProfiler provides a :ref:`benchmark <matbenchmark>` function. This is the main entry point to the package. It benchmarks given solvers on the selected test suite.

We provide below simple examples on how to use OptiProfiler in MATLAB. For more details on the signature of the :ref:`benchmark <matbenchmark>` function, please refer to the :ref:`MATLAB API documentation <matlabapi>`.

Examples
--------

Ex1: first example to try out
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

Let us first try to benchmark two callable optimization solvers **solver1** and **solver2** (e.g., **fminsearch** and **fminunc** in MATLAB Optimization Toolbox) on the default test suite.
(Note that each **solver** must accept signatures mentioned in the `Cautions` part of the :ref:`benchmark <matbenchmark>` function according to the type of problems you want to solve.)

To do this, run:

.. code-block:: matlab

    benchmark({@solver1, @solver2})

This will benchmark the two solvers on the default test suite, which means ``'plain'`` feature (see :ref:`Feature <matfeature>`) and 49 unconstrained problems from the default problem set (see :ref:`s_load <matsload>` or :ref:`s_select <matsselect>`).

There will be a new folder named ``'out'`` in the current working directory, which contains a subfolder named ``plain_<timestamp>`` with all the detailed results and a pdf file named ``summary.pdf`` summarizing all the performance profiles and data profiles.


Ex2: one step further by adding options
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

You can also add options to the benchmark function. For example, if you want to benchmark three solvers **solver1**, **solver2**, and **solver3** on the test suite with the ``'noisy'`` feature and all the unconstrained and bound-constrained problems with dimension between 6 and 10 from the default problem set, you can run:

.. code-block:: matlab

    options.ptype = 'ub';
    options.mindim = 6;
    options.maxdim = 10;
    options.feature_name = 'noisy';
    benchmark({@solver1, @solver2, @solver3}, options)

This will create the corresponding folders ``'out/noisy_<timestamp>'`` and files as in the previous example Ex1, but with the new options. More details on the options can be found in the :ref:`benchmark <matbenchmark>` function documentation.


Ex3: useful option **load**
^^^^^^^^^^^^^^^^^^^^^^^^^^^

OptiProfiler provides a practically useful option named **load**. This option allows you to load the results from a previous benchmarking run (without solving all the problems again) and use them to draw new profiles with different options. For example, if you have just run the second example Ex2 and OptiProfiler has finished the job and successfully created the folder ``'out'`` in the current working directory, you can run:

.. code-block:: matlab

    options.load = 'latest';
    options.solvers_to_load = [1 3];
    options.ptype = 'u';
    options.mindim = 7;
    options.maxdim = 9;
    benchmark(options)

This will directly draw the profiles for the **solver1** and **solver3** with the ``'noisy'`` feature and all the unconstrained problems with dimension between 7 and 9 selected from the previous run. The results will also be saved under the current directory with a new subfolder named ``noisy_<timestamp>`` with the new timestamp.


Ex4: scoring solvers
^^^^^^^^^^^^^^^^^^^^

Besides generating profiles, OptiProfiler also returns scores for the solvers based in the profiles. For example, if you want to score the solvers in Ex1, you can run:

.. code-block:: matlab

    scores = benchmark({@solver1, @solver2})

This will return a vector of scores for the two solvers. The scores by default are calculated based on all the history-based performance profiles and can be modified by the option **score_fun** (see :ref:`benchmark <matbenchmark>` function documentation for more details).