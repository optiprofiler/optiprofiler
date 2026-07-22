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

Problem-library protocol
------------------------

.. autosummary::
    :toctree: generated/

    ProblemLibraryPlugin
    list_problem_libraries

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

External provider APIs are documented by their independently installed
packages: `PyCUTEst <https://github.com/optiprofiler/pycutest>`_,
`SOLAR <https://github.com/optiprofiler/solar_python>`_, and experimental
`RS13 <https://github.com/optiprofiler/rs13>`_.  See
:ref:`python_problem_libraries` for their development installation and
uninstall boundaries.
