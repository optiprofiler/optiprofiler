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

    OptiProfiler includes the `S2MPJ <https://github.com/GrattonToint/S2MPJ>`_ problem library by default. The Python `PyCUTEst <https://github.com/optiprofiler/pycutest>`_, `SOLAR <https://github.com/optiprofiler/solar_python>`_, and experimental `RS13 <https://github.com/optiprofiler/rs13>`_ providers are independently installed plugin packages; their source, runtime, and license files are not included in the OptiProfiler core distribution.

MATLAB
------

1. Clone the repository using the following command:

.. code-block:: bash

    git clone https://github.com/optiprofiler/optiprofiler.git

2. In MATLAB, navigate to the folder where the source code is located, and you will see a file named ``setup.m``. Run the following command in the MATLAB command window:

.. code-block:: matlab

    setup

If the optional SOLAR MATLAB adapter is not installed locally, ``setup``
asks whether to download it. For non-interactive scripts, pass the choice
explicitly:

.. code-block:: matlab

    setup(struct('install_solar', true))
    setup(struct('install_solar', false))

When accepted or requested, ``setup`` clones `optiprofiler/solar_matlab
<https://github.com/optiprofiler/solar_matlab>`_ into the local directory
``matlab/optiprofiler/problem_libs/solar`` and adds it to the MATLAB path.
This directory is a local optional installation, not a submodule of the
OptiProfiler repository. Passing ``install_solar=false`` skips the prompt and
optional download when the adapter is not installed locally, which is useful
in CI and batch scripts.
The public problem-library name remains ``solar``; use ``plibs={'solar'}``
in MATLAB benchmark options.
