# Release Checklist

## Required for the Next Release

- [ ] Build a MATLAB-only ZIP from the same tagged commit used for the Python
  release.
- [ ] Include `setup.m`, the MATLAB engine, `matlab_problem_libraries.lock`,
  required licenses and metadata, and the bundled S2MPJ contents.
- [ ] Obtain and include the upstream S2MPJ license or written redistribution
  permission. `THIRD_PARTY_NOTICES.md` records provenance but does not replace
  that authorization.
- [ ] Exclude Python sources, documentation build inputs, workflows, local
  development files, and Git history.
- [ ] Test installation and the default S2MPJ benchmark from a clean extraction
  of the ZIP on every supported MATLAB platform.
- [ ] Attach the verified ZIP to the release. A MATLAB toolbox package may be
  added separately, but it does not replace the required ZIP.
- [ ] Run `python tools/check_matlab_zip.py --release-ready <archive>` before
  attaching the asset.

The OptiProfiler repository remains the single source of truth. Do not maintain
a second hand-edited MATLAB core repository; generate any distribution mirror
from the tagged source commit.
