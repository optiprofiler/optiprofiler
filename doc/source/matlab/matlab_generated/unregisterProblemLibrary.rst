.. _matunregisterproblemlibrary:

unregisterProblemLibrary
========================

**registration** = **unregisterProblemLibrary(**\ *name*\ **)**
    Removes one persistent external or custom problem-library registration.

-------------------------------------------------------------------------

The function removes the registry entry and the provider root from the current
and persisted MATLAB path.  It returns the removed registration metadata.

It does not delete adapter source, runtime files, caches, compiled artifacts,
or generated data.  Bundled S2MPJ and the bundled custom example cannot be
unregistered.  See :ref:`matlab_problem_libraries` for the distinction between
unregistering a provider, uninstalling OptiProfiler, and manually removing
files.
