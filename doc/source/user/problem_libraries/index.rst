.. _problem_libraries:

Problem libraries
=================

OptiProfiler separates the benchmarking engine from the collections of
optimization problems used by an experiment.  The MATLAB package includes
S2MPJ as the default library.  Other providers are installed outside the core
source and registered under stable public names.

For a typical experiment, the workflow is:

1. install OptiProfiler and any optional provider;
2. select providers by public name in ``options.plibs``;
3. run the benchmark without hard-coding provider paths; and
4. update, detach, or remove each independently owned installation when needed.

The provider pages distinguish the MATLAB engine, the adapter checkout, a
runtime or data installation, the persistent registry, and generated files.
Detaching a provider from MATLAB does not silently delete its source or data.

.. toctree::
    :maxdepth: 1

    matlab
    providers
    custom_matlab
