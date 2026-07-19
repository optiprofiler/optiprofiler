.. _matlab_problem_libraries:

MATLAB problem libraries
========================

Install optional providers
--------------------------

Run ``setup`` from the extracted MATLAB distribution.  S2MPJ is bundled and
requires no registry entry.  During an interactive installation, setup asks
whether optional providers should be installed.  For deterministic scripts,
state every choice explicitly:

.. code-block:: matlab

    setup(struct( ...
        'install_matcutest', false, ...
        'install_solar', true))

Optional adapters are installed under
``prefdir/optiprofiler/problem_libraries`` by default.  Choose another writable
parent directory with ``problem_library_root``:

.. code-block:: matlab

    setup(struct( ...
        'install_matcutest', false, ...
        'install_solar', true, ...
        'problem_library_root', '/path/to/problem-libraries'))

The setup flags mean "install or configure this provider during this run".
Setting a flag to ``false`` skips that provider; it does not unregister it or
delete existing files.

Select a provider
-----------------

Most users only need public names in the benchmark options:

.. code-block:: matlab

    options.plibs = {'s2mpj', 'solar'};
    scores = benchmark({@solver1, @solver2}, options);

Provider functions remain directly callable after setup:

.. code-block:: matlab

    names = solar_select(struct('ptype', 'n', 'maxdim', 20));
    problem = solar_load(names{1});

Generic tools can resolve the same provider without hard-coding function
names:

.. code-block:: matlab

    library = resolveProblemLibrary('solar');
    names = library.select(struct('ptype', 'n', 'maxdim', 20));
    problem = library.load(names{1});

Ordinary benchmark users do not need to call ``resolveProblemLibrary``.
OptiProfiler uses it internally when processing ``options.plibs``.

Storage and registry
--------------------

The default layout is:

.. code-block:: text

    prefdir/optiprofiler/
    ├── problem_libraries.mat
    └── problem_libraries/
        ├── matcutest/
        ├── solar/
        └── runtime/

``problem_libraries.mat`` is the persistent registry.  Adapter checkouts and
the MatCUTEst runtime live below the optional-library root.  A custom
``problem_library_root`` relocates optional checkouts and runtimes, not the
registry.  ``OPTIPROFILER_MATLAB_PROBLEM_LIBRARY_REGISTRY`` is intended only
for an isolated registry, such as a CI test.

Update optional providers
-------------------------

Rerun ``setup`` with the appropriate install flag.  Setup aligns a clean
setup-managed checkout with the provider version selected by the installed
OptiProfiler distribution:

.. code-block:: matlab

    setup(struct( ...
        'install_matcutest', false, ...
        'install_solar', true))

Setup does not overwrite a dirty checkout.  Review, commit, or move local
changes before updating.  This protects user work from a destructive reset.

Detach, uninstall, and remove files
-----------------------------------

.. list-table:: MATLAB removal operations
    :header-rows: 1
    :widths: 25 38 37

    * - Operation
      - Removes
      - Preserves
    * - ``unregisterProblemLibrary('solar')``
      - Registry entry plus current and persisted MATLAB path entries
      - Adapter source, runtime, cache, executable, and generated data
    * - ``setup uninstall``
      - OptiProfiler and setup-managed MATLAB path entries, including
        registered adapter roots
      - Registry entries, optional checkouts, runtimes, caches, executables,
        and generated data
    * - Manual removal after unregistering
      - Only the directory explicitly selected by the user
      - Everything outside that directory

Always unregister an optional or custom provider before manually deleting its
root.  OptiProfiler intentionally has no destructive ``removeProblemLibrary``
operation: an installation may use a custom root, contain local changes, share
a runtime, or have been created outside setup.
