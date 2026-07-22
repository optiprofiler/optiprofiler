.. _python_custom_problem_libraries:

Custom Python problem libraries
===============================

Local filesystem adapter
------------------------

A local adapter can live outside the OptiProfiler source tree.  Create one
directory per public library name:

.. code-block:: text

    /path/to/my-libraries/
    └── myproblems/
        └── myproblems_tools.py

``myproblems_tools.py`` must define ``myproblems_select`` and
``myproblems_load``.  The selector returns problem-name strings and the loader
returns an :class:`~optiprofiler.Problem`.  A basic local adapter uses the
signatures ``myproblems_select(problem_options)`` and
``myproblems_load(problem_name)``.

To accept library-specific ``plib_options``, also define both
``myproblems_get_default_options()`` and
``myproblems_validate_options(options)``.  In that case the selector receives
``(problem_options, library_options)``, and the loader must accept
``library_options`` as a keyword argument.  Point the benchmark to either the
parent directory or the library directory:

.. code-block:: python

    scores = benchmark(
        [solver1, solver2],
        plibs=['s2mpj', 'myproblems'],
        custom_problem_libs_path='/path/to/my-libraries',
    )

An explicit filesystem adapter applies only to that benchmark run and takes
precedence over an installed provider with the same public name.

.. _problem_library_plugin_protocol:

Installable provider protocol
-----------------------------

An installable provider registers one entry point for each public name in the
``optiprofiler.problem_libraries`` group:

.. code-block:: toml

    [project.entry-points."optiprofiler.problem_libraries"]
    myproblems = "optiprofiler_myproblems.plugin:get_problem_library"

The entry point resolves to a zero-argument factory returning
:class:`~optiprofiler.ProblemLibraryPlugin`:

.. code-block:: python

    from optiprofiler import ProblemLibraryPlugin

    def get_problem_library():
        return ProblemLibraryPlugin(
            name='myproblems',
            api_version=1,
            select=select_problems,
            load=load_problem,
            check_available=check_available,
        )

``select(problem_options, library_options)`` returns an ordered sequence of
problem names.  ``load(problem_name, library_options)`` returns a
:class:`~optiprofiler.Problem`.  A configurable provider supplies both
``get_default_options`` and ``validate_options``; omitting both defines a
provider with no library-specific settings.

Discovery reads only entry-point metadata.  OptiProfiler loads the factory and
runs ``check_available`` when the provider is selected.  Factories must be
lightweight and safe in spawned workers, and callbacks must not depend on
mutable state stored in one process.  Provider distributions should declare a
core dependency that includes the released OptiProfiler version implementing
their protocol version.

The entry-point name, ``ProblemLibraryPlugin.name``, and the name used in
``plibs`` must be identical.  OptiProfiler rejects duplicate installed names
and conflicts with bundled providers.
