.. _use:

Usage for MATLAB
================

OptiProfiler provides a :ref:`benchmark <matbenchmark>` function. This is the main entry point to the package. It benchmarks given solvers on the selected test suite.

We provide below simple examples on how to use OptiProfiler in MATLAB. For more details on the signature of the :ref:`benchmark <matbenchmark>` function, please refer to the :ref:`MATLAB API documentation <matlabapi>`.

Examples
--------

First example to try out
^^^^^^^^^^^^^^^^^^^^^^^^

Let us first try to benchmark two callable optimization solvers **solver1** and **solver2** (e.g., **fminsearch** and **fminunc** in MATLAB Optimization Toolbox) on the default test suite.
(Note that each **solver** must accept signatures mentioned in the `Cautions` part of the :ref:`benchmark <matbenchmark>` function according to the type of problem you want to solve.)

To do this, run:

.. code-block:: matlab

    benchmark({@solver1, @solver2})

This will benchmark the two solvers on the default test suite, which means ``'plain'`` feature (see :ref:`Feature <matfeature>`) and 49 unconstrained problems from the default problem set (see :ref:`s_load <matsload>` or :ref:`s_select <matsselect>`).

There will be a new folder named ``'out'`` in the current working directory, which contains a subfolder named ``plain_<timestamp>`` with all the detailed results and a pdf file named ``summary.pdf`` summarizing all the performance profiles and data profiles.


A step further: adding options
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

Useful option: **load**
^^^^^^^^^^^^^^^^^^^^^^^

Scoring solvers
^^^^^^^^^^^^^^^