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

    OptiProfiler includes the `S2MPJ <https://github.com/GrattonToint/S2MPJ>`_ problem library by default. The Python `PyCUTEst <https://github.com/optiprofiler/pycutest>`_, `SOLAR <https://github.com/optiprofiler/solar_python>`_, and experimental `RS13 <https://github.com/optiprofiler/rs13>`_ providers are installed independently and discovered through the versioned problem-library entry-point protocol. They are not part of the OptiProfiler core wheel or source distribution; follow each provider repository for its runtime and license requirements.

MATLAB
------

1. Clone the repository using the following command:

.. code-block:: bash

    git clone https://github.com/optiprofiler/optiprofiler.git

2. In MATLAB, navigate to the root directory of this repository, where you can see a file named ``setup.m``. Run the following command in the MATLAB command window:

.. code-block:: matlab

    setup

The ``setup`` function performs the following tasks:

- Adds the necessary directories to the MATLAB search path.
- Uses the bundled `S2MPJ <https://github.com/GrattonToint/S2MPJ>`_ library and can optionally set up `MatCUTEst <https://github.com/matcutest>`_ and SOLAR at the commits recorded in ``matlab_problem_libraries.lock``.

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
