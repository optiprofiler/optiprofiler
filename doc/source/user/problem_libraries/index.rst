.. _problem_libraries:

Problem libraries
=================

OptiProfiler separates the benchmarking engine from the collections of
optimization problems used by an experiment.  The core distribution includes
S2MPJ as the default library.  Other providers can be installed independently
without copying their source into the core package.

For a typical experiment, the workflow is:

1. install OptiProfiler and any optional provider;
2. inspect the available public library names;
3. pass the selected names to ``plibs``; and
4. update or uninstall each independently owned distribution when needed.

The provider pages distinguish four objects that are easy to confuse: the
OptiProfiler engine, the adapter, an upstream runtime or data collection, and
generated caches or experiment output.  An uninstall operation removes only
the files owned by the selected distribution.

.. toctree::
    :maxdepth: 1

    python
    providers
    custom_python
