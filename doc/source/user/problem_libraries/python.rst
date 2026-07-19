.. _python_problem_libraries:

Python problem libraries
========================

Install a provider
------------------

Install the core package first.  S2MPJ is included and requires no additional
package:

.. code-block:: bash

    python -m pip install optiprofiler

Install only the optional providers that an experiment needs:

.. code-block:: bash

    python -m pip install optiprofiler-pycutest
    python -m pip install optiprofiler-solar
    python -m pip install optiprofiler-rs13

These distributions provide OptiProfiler adapters.  PyCUTEst is the exception
to a self-contained installation: its upstream PyCUTEst and CUTEst runtime must
already be installed and working.  The OptiProfiler adapter deliberately does
not install or remove that system-dependent runtime.  See the
`upstream PyCUTEst installation guide
<https://jfowkes.github.io/pycutest/_build/html/install.html>`_.

The SOLAR adapter ships its slim runtime source and builds a local executable
in a user cache when required.  The experimental RS13 adapter ships the
official Python problem definitions required for normal use.  Their detailed
ownership boundaries are summarized in :ref:`python_problem_library_providers`.

Discover and select providers
-----------------------------

List the public names registered in the current Python environment:

.. code-block:: python

    from optiprofiler import list_problem_libraries

    print(list_problem_libraries())

Listing reads package metadata without importing every optional runtime.  A
provider can therefore be listed even when a separately managed dependency,
such as CUTEst, is unavailable.  Availability is checked when the provider is
selected.

Pass one or more public names to ``plibs``:

.. code-block:: python

    from optiprofiler import benchmark

    scores = benchmark(
        [solver1, solver2],
        plibs=['s2mpj', 'pycutest'],
    )

Library-specific experiment settings belong in ``plib_options``.  Paths to
system runtimes, licenses, compiled executables, and caches remain installation
concerns owned by the provider and are not benchmark options.

Update providers
----------------

Update the core and each independently installed adapter explicitly:

.. code-block:: bash

    python -m pip install --upgrade optiprofiler
    python -m pip install --upgrade optiprofiler-pycutest
    python -m pip install --upgrade optiprofiler-solar
    python -m pip install --upgrade optiprofiler-rs13

Only run the commands for distributions that are installed.  Package metadata
prevents an adapter from being installed with an incompatible core version.
Upstream runtimes that are not owned by an adapter, such as CUTEst, follow
their own update procedure.

Uninstall without deleting user data
------------------------------------

Remove an optional adapter by its distribution name:

.. code-block:: bash

    python -m pip uninstall optiprofiler-pycutest

This removes only that adapter and its entry point.  It does not remove CUTEst,
compiled CUTEst problems, caches, or experiment output.  Use the corresponding
distribution name to remove SOLAR or RS13.

.. code-block:: bash

    python -m pip uninstall optiprofiler

Removing the core also removes bundled S2MPJ.  It does not delete independently
installed adapters, external runtimes, caches, or data.  Those adapters remain
installed but cannot run until a compatible core is installed again.

For SOLAR, uninstalling the adapter removes the slim source installed inside
the distribution but preserves the compiled cache under
``${XDG_CACHE_HOME:-~/.cache}/optiprofiler/solar`` by default.  Removing a
cache, an upstream runtime, or benchmark output is always a separate,
user-directed action.
