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

    OptiProfiler includes the `S2MPJ <https://github.com/GrattonToint/S2MPJ>`_ problem library by default. Optional problem libraries include `PyCUTEst <https://jfowkes.github.io/pycutest/>`_ and the SOLAR adapter. PyCUTEst requires a separate installation. SOLAR uses a slim runtime derived from `bbopt/solar <https://github.com/bbopt/solar>`_ under LGPL-2.1; see the adapter README and runtime manifest for provenance.

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
- Clones the default problem libraries, including `S2MPJ <https://github.com/GrattonToint/S2MPJ>`_, and can optionally set up `MatCUTEst <https://github.com/matcutest>`_ and SOLAR.

Note that the installation of MatCUTEst and the SOLAR MATLAB adapter is optional. During the setup process, you will be asked whether you want to install them. If the ``solar_matlab`` submodule is already populated, ``setup`` adds it to the MATLAB path automatically; passing ``install_solar=false`` skips it even when present. Please be aware that MatCUTEst is only supported on Linux systems and is not available on macOS or Windows. The SOLAR adapter builds a local C++ executable from the vendored slim SOLAR runtime when it is first used.

For automated environments (e.g., CI/CD scripts) where interactive input is not possible, you can bypass the prompt by providing an additional option to the ``setup`` function:

.. code-block:: matlab

    setup(struct('install_matcutest', false, 'install_solar', true))

3. If you want to uninstall the package, you can run:

.. code-block:: matlab

    setup uninstall
