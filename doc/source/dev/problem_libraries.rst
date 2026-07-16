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

The MATLAB implementation must ultimately provide the same engineering goals:
an explicit and testable problem-library contract, a clean boundary between
the default library and optional libraries, reproducible integration checks,
and preservation of default behavior.  It may use MATLAB-native setup,
registry, and path mechanisms rather than Python entry points.  Completion of
the Python physical split therefore does not by itself complete the overall
problem-library modernization.
