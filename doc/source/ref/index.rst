.. _pythonapi:

.. module:: optiprofiler

Python API documentation
========================

:Release: |release|
:Date: |today|

This section references the Python API of OptiProfiler. It contains the detailed documentation of the main function and classes in OptiProfiler. For installation or simple usage examples, please refer to the :ref:`User guide <guide>`.

.. currentmodule:: optiprofiler

Main function for most users
-----------------------------

.. autosummary::
    :toctree: generated/

    benchmark

Main classes
------------

.. autosummary::
    :toctree: generated/

    Problem
    Feature
    FeaturedProblem

Problem library configuration
------------------------------

.. autosummary::
    :toctree: generated/

    get_plib_config
    set_plib_config

Other tools
-----------

.. autosummary::
    :toctree: generated/

    show_versions

.. currentmodule:: optiprofiler.problem_libs.s2mpj

.. autosummary::
    :toctree: generated/

    s2mpj_load
    s2mpj_select

.. currentmodule:: optiprofiler.problem_libs.pycutest

.. autosummary::
    :toctree: generated/

    pycutest_load
    pycutest_select
