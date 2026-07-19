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

The provider repositories contain detailed provenance, licenses, maintenance
checks, and provider-specific implementation notes:

* `S2MPJ MATLAB adapter <https://github.com/optiprofiler/s2mpj_matlab>`_;
* `MatCUTEst adapter <https://github.com/optiprofiler/matcutest>`_; and
* `SOLAR MATLAB adapter <https://github.com/optiprofiler/solar_matlab>`_.
