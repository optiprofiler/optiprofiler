.. _install:

Installation
============

Python
------

Install OptiProfiler from PyPI:

.. code-block:: bash

    pip install optiprofiler

You can also install OptiProfiler from conda-forge:

.. code-block:: bash

    conda install conda-forge::optiprofiler

.. note::

    OptiProfiler includes `S2MPJ <https://github.com/GrattonToint/S2MPJ>`_
    by default. Optional Python problem libraries are installed as independent
    distributions. See :ref:`python_problem_libraries` for the provider list,
    installation commands, runtime requirements, and uninstall boundaries.

MATLAB
------

1. Download the MATLAB-only ZIP from the OptiProfiler release page and extract
   it to a writable directory.  The archive contains ``setup.m``, the MATLAB
   engine, documentation, licenses, and bundled S2MPJ; it does not clone the
   Python package or repository-development files.

2. In MATLAB, navigate to the extracted directory, where ``setup.m`` is
   located, and run:

.. code-block:: matlab

    setup

The ``setup`` function performs the following tasks:

- Adds the necessary directories to the MATLAB search path.
- Uses the bundled `S2MPJ <https://github.com/GrattonToint/S2MPJ>`_ library and
  can optionally set up `MatCUTEst <https://github.com/matcutest>`_ and SOLAR.

Note that the installation of MatCUTEst and the SOLAR MATLAB adapter is optional. During the setup process, you will be asked whether you want to install them. By default, optional repositories are cloned under ``prefdir/optiprofiler/problem_libraries`` rather than into the OptiProfiler source tree. Passing ``install_solar=false`` or ``install_matcutest=false`` skips the corresponding adapter even when it is already present. Please be aware that MatCUTEst is only supported on Linux systems and is not available on macOS or Windows. The SOLAR adapter builds a local C++ executable from its slim SOLAR runtime when it is first used.

For automated environments (e.g., CI/CD scripts) where interactive input is not possible, you can bypass the prompt by providing an additional option to the ``setup`` function:

.. code-block:: matlab

    setup(struct('install_matcutest', false, 'install_solar', true))

To place optional libraries in a specific writable directory, set
``problem_library_root`` explicitly:

.. code-block:: matlab

    setup(struct( ...
        'install_matcutest', false, ...
        'install_solar', true, ...
        'problem_library_root', '/path/to/optiprofiler-problem-libraries'))

3. If you want to uninstall the package, you can run:

.. code-block:: matlab

    setup uninstall

This removes OptiProfiler and setup-managed MATLAB path entries, including
registered adapter roots, but preserves the registry, optional checkouts,
runtimes, caches, executables, and user data.  Use
``unregisterProblemLibrary`` to detach one provider without deleting its files.
See :ref:`matlab_problem_libraries` for provider installation, storage,
updates, and the full removal matrix.
