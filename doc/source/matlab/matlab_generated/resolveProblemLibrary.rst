.. _matresolveproblemlibrary:

resolveProblemLibrary
=====================

**library** = **resolveProblemLibrary(**\ *name*\ **)**
    Resolves and validates one MATLAB problem library by its public name.

-------------------------------------------------------------------------

The returned scalar structure includes canonical metadata and validated
``select`` and ``load`` function handles.  Optional ``collect_info`` and
``check_available`` handles are included when the provider declares them.

The resolver supports bundled providers and persistent external or custom
registrations.  It reports a clear error when a known provider is not
installed, the current platform is unsupported, or a registered function does
not resolve below the recorded root.

Most users select providers through ``options.plibs`` and do not need to call
this function directly.  It is useful for generic tools and direct adapter
inspection.  See :ref:`matlab_problem_libraries` for examples.
