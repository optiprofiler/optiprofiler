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

Sources and provenance
----------------------

The OptiProfiler repositories below are integration layers or maintained
subsets, not replacements for the original projects.  Please follow the
upstream projects' citation and license guidance when using their problem
collections.

* **S2MPJ:** `OptiProfiler adapter and synchronized Python subset
  <https://github.com/optiprofiler/s2mpj_python>`_; `original S2MPJ repository
  <https://github.com/GrattonToint/S2MPJ>`_.
* **PyCUTEst:** `OptiProfiler PyCUTEst adapter
  <https://github.com/optiprofiler/pycutest>`_; `upstream PyCUTEst
  <https://github.com/jfowkes/pycutest>`_; `CUTEst
  <https://github.com/ralna/CUTEst>`_; `SIFDecode
  <https://github.com/ralna/SIFDecode>`_; `MASTSIF problem definitions on
  Bitbucket <https://bitbucket.org/optrove/sif>`_.
* **SOLAR:** `OptiProfiler adapter and pinned slim runtime snapshot
  <https://github.com/optiprofiler/solar_python>`_; `original SOLAR repository
  <https://github.com/bbopt/solar>`_.
* **RS13:** `OptiProfiler adapter and vendored official 502-problem Python runtime
  <https://github.com/optiprofiler/rs13>`_; `official problem-collection source
  <https://minlp.com/black-box-optimization-test-problems>`_; `original paper
  <https://doi.org/10.1007/s10898-012-9951-y>`_.  RS13's official source is a
  project page and archive rather than a public upstream GitHub repository.
