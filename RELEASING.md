# OptiProfiler release acceptance

Changing package metadata is preparation, not publication. Do not create a
release tag until the exact integrated commit and its artifacts are accepted.
The paper maintenance branches keep their own versions and provider pins;
V2 features and release metadata must not be backported to them.

## Required core gates

1. Integrate the reviewed runtime fixes. Require terminal, successful checks
   for the exact final SHA, including Windows spawned-provider tests.
2. Validate both integration locks and their bundled S2MPJ gitlinks. A lock
   schema pass proves consistency, not that newer provider fixes were adopted.
3. Build wheel, sdist and MATLAB ZIP from a clean checkout with the exact
   bundled providers. Preserve schemas, licenses and provenance notices.
4. Run `python tools/check_release.py --dist-dir dist --tag v2.0.0` for a
   V2.0.0 release. The source version, tag, artifact filenames, metadata and
   packaged version must agree. The checker does not import the core package.
5. Install each Python artifact through normal dependency resolution, run
   `pip check`, and exercise `tools/smoke_python_distribution.py` outside the
   source checkout. CI tests wheel and sdist on Linux, macOS and Windows.
   It checks actual solver evaluations, composed features, reference scoring,
   parallel execution, readable PDFs and solver-free save/load with unchanged
   source data. These are headless Agg tests, not native GUI-backend tests.
6. Require clean-extraction MATLAB ZIP tests on all supported platforms and
   S2MPJ coexistence, fresh-process restoration and real process-worker tests.
   A same-session load is not evidence of successful native-state restoration.
7. Run the wider numerical and whole-library acceptance on syu-ubuntu using
   the final SHA, record versions and budgets, and retain reports and figures.
   Earlier candidate results must not be relabeled as final release results.

`build.yml` blocks both publishers until the built-artifact acceptance matrix
passes. A tag whose version differs from the package fails before publication.
An already-existing PyPI artifact is an error rather than a silently skipped
upload. Tag creation and uploading still require an explicit release decision.

## External providers

Python PyCUTEst, SOLAR and RS13 are independent distributions, not core
submodules. Their lock entries identify reviewed source snapshots for
integration testing; the runtime discovers installed entry points and does not
read the lock as an automatic installer. Provision upstream native prerequisites
separately. Do not use `--no-deps` as end-user installation acceptance.

Each provider needs its own source/artifact acceptance, compatible declared
core requirement, license/provenance inventory and publication decision.
Passing core CI does not publish those packages. Documentation must distinguish
source installation from a package that is actually available on PyPI.

Outstanding owner decisions are kept in the release audit, not silently settled
by a technical test: licenses for independently owned adapter/glue files,
MatCUTEst native redistribution materials, and the provider publication scope.
Preserve upstream licenses unchanged.

Before tagging, also settle these prerequisites:

- Every locked provider commit (both lock files and both S2MPJ gitlinks) must
  be reachable from a protected branch or a tag of its provider repository.
  `setup.m` and `git submodule update` fetch these exact commits; a deleted or
  rewritten branch would make the tagged release uninstallable from Git.
- Until a provider is actually published on PyPI, the documentation must give
  the exact-SHA source installation (see the Python problem-library guide) and
  must not present `pip install optiprofiler-<provider>` as available.
- The wheel and sdist ship the S2MPJ licence and notice with the bundled
  subset, and the sdist also ships `THIRD_PARTY_NOTICES.md` and `licenses/`
  (enforced by `tools/check_release.py`). The core package metadata lists only
  the OptiProfiler `LICENSE`; PEP 639 `license-files` needs setuptools 77 or
  later (Python 3.9 or later), so adopt it when Python 3.8 support ends.

## Reference scoring boundary

V2 acceptance requires an explicit provider/user reference to define standard
profile convergence targets, while an omitted reference keeps the V1
solver-history minimum (including the optional plain run). Test both fresh
execution and solver-free save/load; metadata propagation alone is not
acceptance of this behavior.

The record stores the raw `fun` and total `maxcv` with one of two kinds,
`feasible` (`maxcv` is zero) and `merit` (nonnegative `maxcv`), rather than a
cached merit; an optional `x_ref` is metadata only. Legacy kind aliases remain
readable. Targets are converted with the current merit function and each run's
initial violation. Built-in coordinate features keep the record, including a
nonzero violation, because truth is evaluated in the original coordinates;
custom stages and `quantized` with `ground_truth=True` fall back to the
empirical target and persist the first fallback reason.

A load converts targets with the merit of that load. A load without
`merit_fun` of an archive with explicit references recomputes every cached
merit (histories, initial values and the plain baseline) with the default
merit and warns when they change; supply the original `merit_fun` to reproduce
a custom-merit experiment. Archives without references keep reusing their
cached merits, as in V1. A record that the selected merit maps to a non-finite
target is reported before the problem is solved. Require both language
regression suites and final installed-artifact tests to pass before closing
this gate; candidate code alone is not release acceptance.
