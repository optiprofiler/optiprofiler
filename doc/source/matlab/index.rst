.. _matlabapi:

MATLAB API documentation
========================

This section references the MATLAB API of OptiProfiler. It contains the detailed documentation of the main function and classes in OptiProfiler. For installation or simple usage examples, please refer to the :ref:`User guide <guide>`.

Main function for most users
----------------------------

:doc:`matlab_generated/benchmark`
    Function to benchmark solvers under problem libraries with specific features.

.. toctree::
    :hidden:

    matlab_generated/benchmark

Main classes
------------

:doc:`matlab_generated/Problem`
    Class representing a general optimization problem.

:doc:`matlab_generated/Feature`
    Class representing a feature used for performance analysis.

:doc:`matlab_generated/FeaturedProblem`
    Class representing a problem equipped with a specific feature.

.. toctree::
    :hidden:

    matlab_generated/Problem
    matlab_generated/Feature
    matlab_generated/FeaturedProblem

Other tools
-----------

:doc:`matlab_generated/s2mpj_load`
    Function to load a problem from the S2MPJ collection.

:doc:`matlab_generated/s2mpj_select`
    Function to select problems from the S2MPJ collection.

:doc:`matlab_generated/matcutest_load`
    Function to load a problem from the MatCUTEst collection.

:doc:`matlab_generated/matcutest_select`
    Function to select problems from the MatCUTEst collection.

.. toctree::
    :hidden:

    matlab_generated/s2mpj_load
    matlab_generated/s2mpj_select
    matlab_generated/matcutest_load
    matlab_generated/matcutest_select