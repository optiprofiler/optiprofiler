OptiProfiler: a platform for benchmarking optimization solvers
==============================================================

|docs_badge| |codecov_badge|

|py_unit| |py_multi_os| |py_random| |py_stress|

|ml_unit| |ml_multi_os| |ml_random| |ml_stress|

.. |docs_badge| image:: https://img.shields.io/readthedocs/optiprofiler/latest?logo=readthedocs&style=for-the-badge
    :target: http://www.optprof.com
.. |codecov_badge| image:: https://img.shields.io/codecov/c/github/optiprofiler/optiprofiler?style=for-the-badge&logo=codecov
   :target: https://app.codecov.io/github/optiprofiler/optiprofiler/tree/main
.. |py_unit| image:: https://github.com/optiprofiler/optiprofiler/actions/workflows/python-unit_test.yml/badge.svg
   :target: https://github.com/optiprofiler/optiprofiler/actions/workflows/python-unit_test.yml
.. |py_multi_os| image:: https://github.com/optiprofiler/optiprofiler/actions/workflows/python-test_multi-os.yml/badge.svg
   :target: https://github.com/optiprofiler/optiprofiler/actions/workflows/python-test_multi-os.yml
.. |py_random| image:: https://github.com/optiprofiler/optiprofiler/actions/workflows/python-random_test.yml/badge.svg
   :target: https://github.com/optiprofiler/optiprofiler/actions/workflows/python-random_test.yml
.. |py_stress| image:: https://github.com/optiprofiler/optiprofiler/actions/workflows/python-stress_test.yml/badge.svg
   :target: https://github.com/optiprofiler/optiprofiler/actions/workflows/python-stress_test.yml
.. |ml_unit| image:: https://github.com/optiprofiler/optiprofiler/actions/workflows/matlab-unit_test.yml/badge.svg
   :target: https://github.com/optiprofiler/optiprofiler/actions/workflows/matlab-unit_test.yml
.. |ml_multi_os| image:: https://github.com/optiprofiler/optiprofiler/actions/workflows/matlab-test_multi-os.yml/badge.svg
   :target: https://github.com/optiprofiler/optiprofiler/actions/workflows/matlab-test_multi-os.yml
.. |ml_random| image:: https://github.com/optiprofiler/optiprofiler/actions/workflows/matlab-random_test.yml/badge.svg
   :target: https://github.com/optiprofiler/optiprofiler/actions/workflows/matlab-random_test.yml
.. |ml_stress| image:: https://github.com/optiprofiler/optiprofiler/actions/workflows/matlab-stress_test.yml/badge.svg
   :target: https://github.com/optiprofiler/optiprofiler/actions/workflows/matlab-stress_test.yml

🌐 **Website:** `www.optprof.com <https://www.optprof.com>`_

Python
------

Install OptiProfiler from PyPI:

.. code-block:: bash

    pip install optiprofiler

You can also install OptiProfiler from conda-forge:

.. code-block:: bash

    conda install conda-forge::optiprofiler

.. note::

    OptiProfiler includes the `S2MPJ <https://github.com/GrattonToint/S2MPJ>`_
    problem library by default. The Python `PyCUTEst
    <https://github.com/optiprofiler/pycutest>`_, `SOLAR
    <https://github.com/optiprofiler/solar_python>`_, and experimental `RS13
    <https://github.com/optiprofiler/rs13>`_ providers are independently
    installed plugin packages. PyCUTEst and CUTEst remain a separately managed
    system runtime; SOLAR and RS13 ship the problem definitions needed by their
    adapters.

Python package removal follows distribution ownership. Running
``python -m pip uninstall optiprofiler`` removes the core package and its
bundled S2MPJ contents, but it does not remove independently installed
problem-library distributions, their external runtimes, or their caches.
Those plugins require OptiProfiler and become usable again after the core is
reinstalled. Remove one adapter separately with, for example,
``python -m pip uninstall optiprofiler-pycutest``. This does not remove the
independently managed PyCUTEst or CUTEst runtime and problem definitions.

MATLAB
------

Download the MATLAB-only ZIP from the OptiProfiler release page and extract it
to a writable directory. In MATLAB, navigate to the extracted directory, where
``setup.m`` is located, and run:

.. code-block:: matlab

    setup

If the optional SOLAR MATLAB adapter is not installed locally, ``setup``
asks whether to download it. For non-interactive scripts, pass the choice
explicitly:

.. code-block:: matlab

    setup(struct('install_solar', true))
    setup(struct('install_solar', false))

When accepted or requested, ``setup`` clones `optiprofiler/solar_matlab
<https://github.com/optiprofiler/solar_matlab>`_ under
``prefdir/optiprofiler/problem_libraries`` by default, adds the adapter root to
the MATLAB path, and registers its canonical functions. Pass
``problem_library_root`` to choose another writable parent directory. The
adapter is not a submodule of the OptiProfiler repository. Passing
``install_solar=false`` skips the prompt and optional download when the adapter
is not installed locally, which is useful in CI and batch scripts.
The public problem-library name remains ``solar``; use ``plibs={'solar'}``
in MATLAB benchmark options. The library-specific ``solar_select`` and
``solar_load`` functions remain callable after setup. Generic code may instead
call ``library = resolveProblemLibrary('solar')`` and then use
``library.select`` and ``library.load``.

To detach an optional or custom library without deleting its source, runtime,
cache, or compiled artifacts, call
``unregisterProblemLibrary('solar')``. Running ``setup uninstall`` removes the
OptiProfiler core and setup-managed search paths, including registered adapter
roots, but preserves external-library registrations and data.
The complete workflow, including storage locations, updates, and the
difference between detaching and deleting files, is documented under
``doc/source/user/problem_libraries/``.
