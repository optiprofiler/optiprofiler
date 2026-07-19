.. _python_problem_libraries:

Python problem-library lifecycle
================================

OptiProfiler keeps the benchmarking engine separate from optional problem
libraries.  The core ``optiprofiler`` distribution includes S2MPJ as its
bundled default.  PyCUTEst, SOLAR, and the experimental RS13 provider use
independent adapter distributions, so installing or removing one of them does
not rewrite the core package.

.. warning::

    The independent Python adapter distributions described on this page are
    development artifacts and are not part of the current stable PyPI release.
    Until they are published together with a compatible OptiProfiler release,
    install the core and adapters from matching source checkouts.  The tested
    repository commits are recorded in ``problem_libraries.lock`` at the
    OptiProfiler repository root.

What owns what
--------------

Four kinds of files may be involved in one problem library:

1. the OptiProfiler core engine;
2. the OptiProfiler adapter distribution;
3. an upstream runtime or problem-data installation; and
4. generated executables, caches, and experiment output.

These have separate owners and therefore separate update and removal steps.

.. list-table:: Python provider ownership
    :header-rows: 1
    :widths: 12 22 30 36

    * - Public name
      - Adapter distribution
      - Runtime or problem data
      - Lifecycle boundary
    * - ``s2mpj``
      - Bundled in ``optiprofiler``
      - Bundled S2MPJ subset
      - Installed and removed with the core distribution.
    * - ``pycutest``
      - ``optiprofiler-pycutest``
      - User-managed PyCUTEst, CUTEst, MASTSIF, and compiled problems
      - The adapter never installs or removes the upstream runtime.
    * - ``solar``
      - ``optiprofiler-solar``
      - Slim SOLAR source shipped by the adapter; executable built in a user
        cache by default
      - Removing the adapter removes its packaged source but preserves the
        external cache and any user-supplied executable.
    * - ``rs13``
      - ``optiprofiler-rs13`` (experimental)
      - User-supplied official ``rs13pm`` files selected by ``RS13PM_DIR``
      - The adapter does not download or delete the upstream archive, extracted
        source, solution data, or audit data.

Install a development adapter
-----------------------------

Use a fresh environment and install a matching OptiProfiler checkout first.
For example, with sibling source checkouts:

.. code-block:: bash

    python -m pip install /path/to/optiprofiler
    python -m pip install /path/to/pycutest

Replace the second path with ``solar_python`` or ``rs13`` for another adapter.
The source repository name and installed distribution name are intentionally
different.  For example, the ``pycutest`` repository installs
``optiprofiler-pycutest`` and the import package ``optiprofiler_pycutest``;
it does not install the upstream module named ``pycutest``.

PyCUTEst requires a separate upstream runtime installation.  Follow the
`upstream PyCUTEst installation guide
<https://jfowkes.github.io/pycutest/_build/html/install.html>`_ and verify that
the upstream ``pycutest`` module can load a problem before selecting the
OptiProfiler provider.  OptiProfiler deliberately does not automate this
runtime installation.

Discover and select a provider
------------------------------

Listing reads package metadata and does not import optional runtimes:

.. code-block:: python

    from optiprofiler import list_problem_libraries

    print(list_problem_libraries())

The result contains public names such as ``s2mpj`` or ``pycutest``.  A listed
provider may still report a missing compiler, runtime, executable, or data
directory when it is selected.  Runtime availability is checked lazily so that
one unavailable optional provider does not prevent importing OptiProfiler.

Select one or more providers with ``plibs``:

.. code-block:: python

    from optiprofiler import benchmark

    scores = benchmark(
        [solver1, solver2],
        plibs=['s2mpj', 'pycutest'],
    )

Library-specific experiment settings belong in ``plib_options``.  Runtime
paths, licenses, compiled executables, and caches remain installation concerns
owned by the provider and should not be placed in ``plib_options``.

Update an adapter
-----------------

During development, update the provider checkout to the commit selected by the
core ``problem_libraries.lock`` and reinstall it:

.. code-block:: bash

    git -C /path/to/pycutest fetch
    git -C /path/to/pycutest checkout <locked-commit>
    python -m pip install --upgrade /path/to/pycutest

Do not assume that the newest commit of every provider is compatible with an
older OptiProfiler core.  After the distributions are formally published, the
same ownership rule will apply to package-index upgrades: update the core and
each independently installed adapter explicitly.

Uninstall without deleting user data
------------------------------------

Python uninstall follows distribution ownership:

.. code-block:: bash

    python -m pip uninstall optiprofiler-pycutest

This removes only the adapter distribution and its entry point.  It does not
remove PyCUTEst, CUTEst, MASTSIF, compiled CUTEst problems, or user caches.
Use the corresponding distribution name to remove another adapter.

.. code-block:: bash

    python -m pip uninstall optiprofiler

Removing the core deletes the engine and bundled S2MPJ, but it does not delete
independently installed adapters, their runtimes, caches, or data.  Those
adapters will be unusable until a compatible core is installed again.

For SOLAR, uninstalling ``optiprofiler-solar`` also removes the slim runtime
source installed inside that distribution.  The default compiled cache under
``${XDG_CACHE_HOME:-~/.cache}/optiprofiler/solar`` is deliberately preserved.
Set ``SOLAR_CACHE_DIR`` when a different cache root is required.  Removing a
cache or an upstream runtime is a separate, user-directed operation and is not
performed by OptiProfiler.
