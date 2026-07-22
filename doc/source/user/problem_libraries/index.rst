.. _problem_libraries:

Problem libraries
=================

OptiProfiler separates the benchmarking engine from the collections of
optimization problems used by an experiment.  The Python and MATLAB
distributions include S2MPJ as the default library.  Other providers can be
installed independently without copying their source into the core package.

For a typical experiment, the workflow is:

1. install OptiProfiler and any optional provider;
2. inspect or select providers by their stable public names;
3. pass those names to Python ``plibs`` or MATLAB ``options.plibs``; and
4. update, detach, uninstall, or remove files according to each provider's
   ownership boundary.

The provider pages distinguish the OptiProfiler engine, the adapter or
distribution, an upstream runtime or data collection, the MATLAB registry, and
generated caches or experiment output.  Uninstall and detach operations affect
only their documented ownership boundary; they do not silently delete
independently managed source, runtimes, or user data.

.. toctree::
    :maxdepth: 1

    python
    providers
    custom_python
    matlab
    matlab_providers
    custom_matlab
