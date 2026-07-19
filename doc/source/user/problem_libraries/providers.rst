.. _python_problem_library_providers:

Python providers
================

.. list-table:: Provider ownership and installation
    :header-rows: 1
    :widths: 12 24 30 34

    * - Public name
      - Distribution
      - Runtime or problem definitions
      - Removal boundary
    * - ``s2mpj``
      - Bundled in ``optiprofiler``
      - Bundled S2MPJ subset
      - Installed and removed with the core.
    * - ``pycutest``
      - ``optiprofiler-pycutest``
      - User-managed PyCUTEst, CUTEst, MASTSIF, and compiled problems
      - Removing the adapter preserves the complete upstream runtime.
    * - ``solar``
      - ``optiprofiler-solar``
      - Slim SOLAR source in the adapter; executable in a user cache
      - Removing the adapter preserves its external cache and user output.
    * - ``rs13``
      - ``optiprofiler-rs13`` (experimental)
      - Official RS13 Python problem definitions shipped with the adapter
      - Removing the adapter removes its packaged problems but preserves
        benchmark output and unrelated user data.

S2MPJ
-----

S2MPJ is the default provider and needs no separate installation.  Its public
selection and loading functions are :func:`~optiprofiler.problem_libs.s2mpj.s2mpj_select`
and :func:`~optiprofiler.problem_libs.s2mpj.s2mpj_load`.

PyCUTEst
--------

The PyCUTEst adapter connects OptiProfiler to an existing upstream PyCUTEst
installation.  Install and validate PyCUTEst first, then install
``optiprofiler-pycutest``.  Platform-specific compiler, CUTEst, and MASTSIF
setup remains entirely upstream-owned.

SOLAR
-----

The SOLAR adapter includes the slim runtime source used by its wrapper.  It
builds the executable lazily and caches the result outside the Python
distribution.  Set ``SOLAR_CACHE_DIR`` only when a non-default cache root is
required.

RS13
----

RS13 is an experimental provider for the Rios--Sahinidis derivative-free
optimization test set.  Its distribution includes the official Python problem
definitions required by ``rs13_load`` together with adapter-maintained
selection metadata.  Complete solution vectors and independent upstream audit
records are maintainer test inputs; they are not required for a benchmark and
are not installed for users.

The external provider repositories contain their detailed provenance,
licenses, maintenance checks, and provider-specific APIs:

* `PyCUTEst adapter <https://github.com/optiprofiler/pycutest>`_;
* `SOLAR adapter <https://github.com/optiprofiler/solar_python>`_; and
* `RS13 adapter <https://github.com/optiprofiler/rs13>`_.
