# Third-party notices

OptiProfiler's root `LICENSE` applies to OptiProfiler-authored code. It does
not change the license of third-party problem definitions or supporting files.

## S2MPJ MATLAB subset

The MATLAB distribution bundles a pinned subset from these repositories:

- OptiProfiler adapter and synchronized subset:
  <https://github.com/optiprofiler/s2mpj_matlab>
- Original S2MPJ source:
  <https://github.com/GrattonToint/S2MPJ>

The exact adapter revision is recorded in `matlab_problem_libraries.lock`.
The bundled files include MATLAB problem definitions, `s2mpjlib.m`, the
problem list, selection metadata, and the OptiProfiler adapter.

Please cite:

> S. Gratton and Ph. L. Toint, S2MPJ: An automatic converter from
> CUTEst-SIF optimization problems to Matlab, Python and Julia,
> Optimization Methods and Software, 2025.
> <https://doi.org/10.1080/10556788.2025.2490640>

As of 2026-07-20, the original S2MPJ repository does not contain a standalone
`LICENSE`, `COPYING`, or `NOTICE` file. Before a public release artifact that
contains this subset is distributed, the maintainers must obtain and include
the applicable upstream license or written redistribution permission. This
notice records provenance; it does not create redistribution rights.
