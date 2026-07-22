.. _matlab_problem_library_providers:

MATLAB providers
================

.. list-table:: Provider ownership and installation
    :header-rows: 1
    :widths: 12 26 28 34

    * - Public name
      - Adapter source
      - Runtime or problem definitions
      - Removal boundary
    * - ``s2mpj``
      - Bundled with OptiProfiler
      - Bundled S2MPJ subset
      - Installed and removed with the core; it cannot be unregistered.
    * - ``matcutest``
      - Setup-managed ``optiprofiler/matcutest`` checkout
      - MatCUTEst runtime under the optional-library root
      - Linux only; registry removal and file removal are separate actions.
    * - ``solar``
      - Setup-managed ``optiprofiler/solar_matlab`` checkout
      - Slim SOLAR source and locally compiled executable
      - Detaching preserves the checkout, runtime, executable, and output.

S2MPJ
-----

S2MPJ is the bundled default and needs no separate download.  Use
``s2mpj_select`` and ``s2mpj_load`` directly, or select the public name
``'s2mpj'`` through ``options.plibs``.

MatCUTEst
---------

MatCUTEst is optional and supported only on Linux.  Request it during setup
with ``install_matcutest=true``.  Setup installs the locked adapter and runtime
under the optional-library root and registers ``matcutest_select`` and
``matcutest_load``.

SOLAR
-----

SOLAR is optional and supported on all MATLAB platforms targeted by
OptiProfiler.  Request it with ``install_solar=true``.  The adapter includes a
slim runtime and builds its executable locally when first needed.

Sources and provenance
----------------------

The OptiProfiler repositories below are integration layers or maintained
subsets, not replacements for the original projects.  Please follow the
upstream projects' citation and license guidance when using their problem
collections.

* **S2MPJ:** `OptiProfiler adapter and synchronized MATLAB subset
  <https://github.com/optiprofiler/s2mpj_matlab>`_; `original S2MPJ repository
  <https://github.com/GrattonToint/S2MPJ>`_.
* **MatCUTEst:** `OptiProfiler adapter and maintained precompiled mirror
  <https://github.com/optiprofiler/matcutest>`_; `official precompiled
  MatCUTEst distribution <https://github.com/matcutest/matcutest_compiled>`_;
  `MatCUTEst source project <https://github.com/matcutest/matcutest>`_;
  `underlying CUTEst engine <https://github.com/ralna/CUTEst>`_.
* **SOLAR:** `OptiProfiler MATLAB adapter and pinned slim runtime snapshot
  <https://github.com/optiprofiler/solar_matlab>`_; `original SOLAR repository
  <https://github.com/bbopt/solar>`_.
