# Third-party notices

OptiProfiler's root `LICENSE` applies to OptiProfiler-authored code. It does
not change the license of third-party problem definitions or supporting files.

## S2MPJ Python and MATLAB subsets

The Python and MATLAB distributions bundle pinned subsets from these repositories:

- OptiProfiler Python adapter and synchronized subset:
  <https://github.com/optiprofiler/s2mpj_python>
- OptiProfiler MATLAB adapter and synchronized subset:
  <https://github.com/optiprofiler/s2mpj_matlab>
- Original S2MPJ source:
  <https://github.com/GrattonToint/S2MPJ>

The exact adapter revisions are recorded in `problem_libraries.lock` and
`matlab_problem_libraries.lock`, respectively. The bundled files include
problem definitions, supporting routines, selection metadata, and the
OptiProfiler adapters. Each subset's `THIRD_PARTY_NOTICES.md` records its
source provenance and retained local differences.

Please cite:

> S. Gratton and Ph. L. Toint, S2MPJ: An automatic converter from
> CUTEst-SIF optimization problems to Matlab, Python and Julia,
> Optimization Methods and Software, 2025.
> <https://doi.org/10.1080/10556788.2025.2490640>

S2MPJ is distributed under the BSD 3-Clause terms in upstream `LICENCE.txt`,
copyright (c) 2026, S. Gratton and Ph. L. Toint. The unmodified text is included
in each bundled subset and in `licenses/S2MPJ-LICENSE.txt`.

License provenance:
<https://github.com/GrattonToint/S2MPJ/blob/fea6a70048eaad28b13a08703ddbfdbf65cd9c30/LICENCE.txt>.
The SHA-256 of the license file is
`a8636fc42ac474fc85fbf451c6a0316f6cbd9efa9031d549797dec6b43e9e5b4`.
This revision identifies the license text, not a wholesale update of the
bundled problem sources. The upstream license does not assign ownership of
the independently maintained OptiProfiler adapters.
