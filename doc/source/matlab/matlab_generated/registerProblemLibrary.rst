.. _matregisterproblemlibrary:

registerProblemLibrary
======================

**library** = **registerProblemLibrary(**\ *registration*\ **)**
    Validates and persists one external or custom MATLAB problem-library
    registration and returns its resolved metadata and function handles.

-------------------------------------------------------------------------

``registration`` is a scalar structure with these required fields:

``name``
    Public provider name used in ``options.plibs``.

``root``
    Existing adapter root directory.

``select_function``
    Name of the canonical selector below ``root``.

``load_function``
    Name of the canonical loader below ``root``.

Optional fields are ``api_version``, ``role``, ``collect_info_function``,
``check_available_function``, and ``platforms``.  An identical registration is
idempotent.  If the public name already has different metadata, unregister it
before registering a replacement.

The registry stores only metadata and a root path.  This function does not copy
or delete adapter source, runtime files, caches, compiled artifacts, or user
data.  Bundled providers such as ``s2mpj`` cannot be registered externally.

See :ref:`matlab_problem_libraries` for the complete lifecycle.
