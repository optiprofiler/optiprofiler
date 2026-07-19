.. _matlab_problem_libraries:

MATLAB problem-library lifecycle
================================

The MATLAB development branch keeps optional problem libraries outside the
OptiProfiler source tree.  ``setup`` installs or locates an adapter,
``registerProblemLibrary`` records its canonical functions, and
``resolveProblemLibrary`` returns validated function handles.  S2MPJ remains
the bundled default and requires no registry entry.

.. warning::

    This page describes the MATLAB development branch.  The current stable
    website continues to document the released package.  A MATLAB-only ZIP is
    a required deliverable for the next release, but is produced in Stage 4B;
    development installation still starts from a repository checkout.

What owns what
--------------

The core source, optional adapter checkout, runtime, registry, and generated
data are separate objects.  Detaching an adapter from MATLAB must not silently
delete the other objects.

.. list-table:: MATLAB provider ownership
    :header-rows: 1
    :widths: 12 25 28 35

    * - Public name
      - Adapter source
      - Runtime or data
      - Lifecycle boundary
    * - ``s2mpj``
      - Bundled under the OptiProfiler MATLAB source
      - Bundled S2MPJ subset
      - Installed and removed with the core; it cannot be registered or
        unregistered as an external provider.
    * - ``matcutest``
      - Optional ``optiprofiler/matcutest`` checkout
      - MatCUTEst runtime installed under the optional-library root
      - Linux only.  Registry removal and file removal are separate actions.
    * - ``solar``
      - Optional ``optiprofiler/solar_matlab`` checkout
      - Slim SOLAR runtime and locally compiled executable
      - Registry removal preserves the checkout, runtime, and executable.
    * - Custom name
      - User-owned directory
      - User-owned problem definitions and generated data
      - OptiProfiler records only the root and canonical function names; it
        never owns or deletes the directory.

Install the core and optional adapters
--------------------------------------

Clone the development repository, start MATLAB in its root directory, and run
``setup``.  For a deterministic noninteractive setup, state each optional
choice explicitly:

.. code-block:: matlab

    setup(struct( ...
        'install_matcutest', false, ...
        'install_solar', true))

S2MPJ is used from the bundled source.  Optional adapters are cloned under
``prefdir/optiprofiler/problem_libraries`` by default.  Choose another writable
parent directory with ``problem_library_root``:

.. code-block:: matlab

    setup(struct( ...
        'install_matcutest', false, ...
        'install_solar', true, ...
        'problem_library_root', '/path/to/problem-libraries'))

The setup flags mean "install or configure this provider during this run".
Setting ``install_solar`` or ``install_matcutest`` to ``false`` skips that
provider; it does not unregister it or delete files.

Storage and registry locations
------------------------------

With the default setup root, the relevant locations are:

.. code-block:: text

    prefdir/optiprofiler/
    |-- problem_libraries.mat
    `-- problem_libraries/
        |-- matcutest/
        |-- solar/
        `-- runtime/

``problem_libraries.mat`` is the persistent registry.  The adapter checkouts
and MatCUTEst runtime live below the optional-library root.  A custom
``problem_library_root`` relocates the optional checkouts and runtime, not the
registry file.  Set ``OPTIPROFILER_MATLAB_PROBLEM_LIBRARY_REGISTRY`` only when
an isolated or relocated registry is needed, for example in CI.

Select and call a provider
--------------------------

Most users only need the public name in the benchmark options:

.. code-block:: matlab

    options.plibs = {'s2mpj', 'solar'};
    scores = benchmark({@solver1, @solver2}, options);

The library-specific functions remain directly callable after setup:

.. code-block:: matlab

    names = solar_select(struct('ptype', 'n', 'maxdim', 20));
    problem = solar_load(names{1});

Generic tools can resolve the same provider without hard-coding its function
names:

.. code-block:: matlab

    library = resolveProblemLibrary('solar');
    names = library.select(struct('ptype', 'n', 'maxdim', 20));
    problem = library.load(names{1});

Ordinary benchmark users do not need to call ``resolveProblemLibrary``.
OptiProfiler uses the resolver internally when it processes ``options.plibs``.

Register a custom provider
--------------------------

A custom adapter may live in any user-owned writable directory.  Register the
root and canonical function names once:

.. code-block:: matlab

    registration = struct( ...
        'name', 'myproblems', ...
        'root', '/path/to/myproblems', ...
        'select_function', 'myproblems_select', ...
        'load_function', 'myproblems_load');
    registerProblemLibrary(registration);

An exactly identical registration is idempotent.  The same public name cannot
silently be redirected to a different root or different function metadata;
call ``unregisterProblemLibrary`` before registering a replacement.

Update a locked optional provider
---------------------------------

The repository-root ``matlab_problem_libraries.lock`` records the tested
adapter commits and canonical functions.  Rerun ``setup`` with the appropriate
install flag to align a clean optional checkout to the locked commit:

.. code-block:: matlab

    setup(struct( ...
        'install_matcutest', false, ...
        'install_solar', true))

Setup does not overwrite a dirty checkout.  If local modifications are found,
review, commit, or move them before updating.  This protects user work and
keeps a lock update from becoming a destructive reset.

Detach, uninstall, and remove files
-----------------------------------

These operations have deliberately different scopes:

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
    * - Manual file removal after unregistering
      - Only the directory explicitly chosen by the user
      - Everything outside that directory

Always unregister an optional or custom provider before manually deleting its
root.  OptiProfiler intentionally has no destructive ``removeProblemLibrary``
operation in this version because setup may reuse a custom root, a dirty
checkout, a shared MatCUTEst runtime, or an installation it did not create.
Likewise, ``setup uninstall`` accepts no install options and does not cascade
into external-data deletion.
