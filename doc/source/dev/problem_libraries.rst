Problem-library integration
===========================

OptiProfiler separates runtime discovery from repository integration.  Python
runtime code discovers independently installed providers through the
``optiprofiler.problem_libraries`` entry-point group and validates API version
1.  It does not read ``problem_libraries.lock``.

Repository lock
---------------

The repository-root ``problem_libraries.lock`` is the machine-readable source
for integration testing.  It records the exact repository and commit tested
for each library, its role, API version, and external package metadata.  Run
the validator after changing the protocol or any pinned library::

    python tools/check_problem_libraries_lock.py

The validator checks the lock schema, the core API version, and the bundled
S2MPJ gitlink.  The External Integration workflow derives its test matrix
directly from the lock, builds every external provider at the pinned commit,
installs it without its optional runtime, and executes its API-v1 protocol
tests.

Physical package boundary
-------------------------

S2MPJ is the one bundled default provider.  PyCUTEst and SOLAR are external
providers installed through their own distributions; they have no core
gitlinks and their sources are not included in core build artifacts.  RS13 is
external and experimental and likewise exists only in the lock-driven
integration matrix.

Each external library owns its build and runtime tests.  The core integration
workflow checks only the cross-repository boundary: locked source identity,
installed entry-point metadata, API compatibility, selection/loading adapter
behavior, and spawned-process safety.  This avoids making the core repository
responsible for downloading every optional runtime.

MATLAB parity
-------------

MATLAB uses its own lock because Python distribution and entry-point metadata
do not apply to MATLAB.  The repository-root
``matlab_problem_libraries.lock`` records the S2MPJ, MatCUTEst, and SOLAR
repositories, exact commits, canonical select/load functions, and supported
platforms.  Validate it with::

    python tools/check_matlab_problem_libraries_lock.py

At runtime, ``resolveProblemLibrary`` returns validated function handles and
``registerProblemLibrary`` persists the root and function names of an optional
or custom library.  The registry is stored under the MATLAB preferences
directory by default.  Set
``OPTIPROFILER_MATLAB_PROBLEM_LIBRARY_REGISTRY`` to isolate or relocate it.

``setup`` reads the MATLAB lock, checks out clean optional repositories at the
locked commits, and registers installed adapters without overwriting a dirty
checkout.  S2MPJ remains the bundled default.  A core-local folder fallback is
temporarily retained for compatibility, but new libraries should use the
explicit registry and may live outside the OptiProfiler source tree.
