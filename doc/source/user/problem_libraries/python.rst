.. _python_problem_libraries:

Python problem libraries
========================

.. note::

   The external providers are separate distributions and are not yet published
   on PyPI. Install them from their repositories at the commits recorded in
   ``problem_libraries.lock``, as shown below; these are the revisions tested
   with this core release. The lock is a tested integration manifest, not an
   automatic pin applied by ordinary ``pip install``. Record the exact core and
   provider revisions, and separately installed runtime versions, with your
   results. Installing an adapter with ``--no-deps`` after preparing its
   prerequisites is a controlled development check, not certification of the
   public dependency chain.

Install a provider
------------------

Install the core package first.  S2MPJ is included and requires no additional
package:

.. code-block:: bash

    python -m pip install optiprofiler

Install only the optional providers that an experiment needs, at the exact
tested commits.  Each provider declares the compatible core version and pip
resolves its other dependencies normally:

.. code-block:: bash

    python -m pip install "optiprofiler-pycutest @ git+https://github.com/optiprofiler/pycutest@38119c384901f26e35eacca78ddd9d2ac3b47265"
    python -m pip install "optiprofiler-solar @ git+https://github.com/optiprofiler/solar_python@2478fe964ce4b697cc5ef2129ed5e7102dc7136f"
    python -m pip install "optiprofiler-rs13 @ git+https://github.com/optiprofiler/rs13@ddc11ff94f1d1df128a4f2c999c5267b3ad1818b"

Once a provider is published on PyPI, ``python -m pip install
optiprofiler-<name>`` installs its released version instead.

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

Update the core and each independently installed adapter explicitly.  After
upgrading the core, reinstall each adapter that is installed at the commit
listed in the new release's ``problem_libraries.lock`` with the command shown
in `Install a provider`_ (or, once published, ``python -m pip install
--upgrade optiprofiler-<name>``):

.. code-block:: bash

    python -m pip install --upgrade optiprofiler

With normal dependency resolution, package metadata declares the compatible
core versions; it does not verify a source commit or replace an installation
test.
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
