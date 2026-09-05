# Third-party notices

OptiProfiler's root `LICENSE` applies to OptiProfiler-authored code. It does
not change the license of third-party problem definitions or supporting files.

## S2MPJ Python and MATLAB subsets

This paper maintenance version retains its pinned problem sources and adapter
interfaces. Exact Python and MATLAB subset revisions are recorded by the Git
submodules; `setup.m` also pins the MATLAB installation revision.

- Python subset: <https://github.com/optiprofiler/s2mpj_python>
- MATLAB subset: <https://github.com/optiprofiler/s2mpj_matlab>
- Original source: <https://github.com/GrattonToint/S2MPJ>

Please cite S. Gratton and Ph. L. Toint, *S2MPJ: An automatic converter from
CUTEst-SIF optimization problems to Matlab, Python and Julia*, Optimization
Methods and Software, 2025, <https://doi.org/10.1080/10556788.2025.2490640>.

The unmodified upstream BSD 3-Clause license, copyright (c) 2026, S. Gratton
and Ph. L. Toint, is included as `LICENCE.txt` in each subset and as
`licenses/S2MPJ-LICENSE.txt` here. Its source is
<https://github.com/GrattonToint/S2MPJ/blob/fea6a70048eaad28b13a08703ddbfdbf65cd9c30/LICENCE.txt>,
with SHA-256
`a8636fc42ac474fc85fbf451c6a0316f6cbd9efa9031d549797dec6b43e9e5b4`.
That revision identifies the license text, not a wholesale update of the
problem corpus. Each subset's notice describes retained source provenance.
The upstream license does not assign ownership of the independently maintained
OptiProfiler adapters.
