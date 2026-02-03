.. _install:

Installation
============

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
- Clones the default problem libraries, including `S2MPJ <https://github.com/GrattonToint/S2MPJ>`_ and `MatCUTEst <https://github.com/matcutest>`_.

Note that the installation of MatCUTEst is optional. During the setup process, you will be asked whether you want to install it. Please be aware that MatCUTEst is only supported on Linux systems and is not available on macOS or Windows.

For automated environments (e.g., CI/CD scripts) where interactive input is not possible, you can bypass the prompt by providing an additional option to the ``setup`` function:

.. code-block:: matlab

    setup(struct('install_matcutest', true))  % Or false if you do not need MatCUTEst

3. If you want to uninstall the package, you can run:

.. code-block:: matlab

    setup uninstall

Python
------

Work in progress.
