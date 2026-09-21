# Frozen 1.x FeaturedProblem fixture

`b0-featured-problem.mat` holds two live `FeaturedProblem` objects saved by the
MATLAB implementation at source `1c4b7d9170ad40e1112e33e5a5638cecc9e644eb`
(the `main` baseline before Feature composition and EvalReport), MATLAB R2026a
on syu-ubuntu:

- `fp_noisy`: `FeaturedProblem(problem, Feature('noisy', 'noise_level', 0.5,
  'n_runs', 2), 6, 17)` after two objective evaluations (`[1; 2]`, `[0.5; 1]`)
  and one inequality-constraint evaluation (`[1; 2]`);
- `fp_plain`: `FeaturedProblem(problem, Feature('plain'), 6, 17)` after one
  objective evaluation (`[1; 2]`);

where `problem = Problem(struct('fun', fun, 'x0', [1; 2], 'cub', cub, 'name',
'b0fixture'))` with `fun = str2func('@(x) sum(x.^2)')` and
`cub = str2func('@(x) x(1) - 0.5')`. Their `feature` properties are the
1.x name/options layout, which the current `Feature.loadobj` returns as a
`LegacyFeatureEnvelope`; `FeaturedProblem.loadobj` rebuilds the single-feature
runtime through `importLegacyFeature`.

`b0-featured-problem-continuation.mat` records what the same live objects
returned next under the 1.x implementation (`fun([0.25; 0.5])`,
`cub([0.25; 0.5])`, `maxcv([0.25; 0.5])`, the histories and counters), which a
restored object must reproduce.

SHA-256:

- `b0-featured-problem.mat`
  `f878e3a38ae1b8e2406ef86ea9f82c603737e5508f2f954d3463e49fbf2b2181`
- `b0-featured-problem-continuation.mat`
  `ab677e7560098a78464cde3976aee5e7eca79c6c22cb17409518861d7c49aa6a`
- `gen_b0_fixture.m` (the generator, in this folder)
  `b5c864c71fbdb229588de325d8d399fbc3f7b2c16078af190a08a59e014e8990`

## Provenance

Two captures, both by the unmodified `1c4b7d9` tree with an isolated path
(`restoredefaultpath`, then the `src` folder of that tree only), MATLAB R2026a
(26.1.0.3203278) on syu-ubuntu.

1. 2026-09-17, `b0-featured-problem.mat` with SHA-256
   `5bfb0da47db2af835e351aad12d5510bb5a44d68de92aaa07b9b8c285a532513`, and the
   continuation file above. The two callbacks were anonymous functions written
   in the generator. MATLAB records the absolute path of the file an anonymous
   function is written in, here
   `~/audits/op-feature-report-stabilization-20260917/evidence/gen_b0_fixture.m`
   under the home directory of the account that ran it,
   and warns `MATLAB:dispatcher:UnresolvedFunctionHandle` when the file is
   loaded on a machine without that path. The tests that require a silent load
   therefore passed on syu-ubuntu and failed on every clean machine (GitHub run
   35536565688). That file is superseded; it is in the history of this
   repository (commit `82013d8`).
2. 2026-09-21, the `b0-featured-problem.mat` of this folder. Same tree, same
   MATLAB, same commands, by `gen_b0_fixture.m` of this folder, whose only
   difference from the first generator is that the callbacks are built with
   `str2func`, which records no file: `functions(h).file` is empty for all
   eight handles in the file, so nothing is looked up when it is loaded.
   The continuation that this capture wrote is identical to the one of the
   first capture in all eleven fields (`isequaln`), so the continuation file
   of 2026-09-17 is kept as it is.

Never regenerate this fixture with a newer implementation: it is evidence of
what 1.x saved. Native MAT loading is a trusted-input boundary.

## Regenerating it (with the 1.x implementation only)

From a checkout of this repository, with `OUT` an empty directory:

```
git worktree add --detach /tmp/op-b0-tree 1c4b7d9170ad40e1112e33e5a5638cecc9e644eb
matlab -batch "restoredefaultpath; addpath('/tmp/op-b0-tree/matlab/optiprofiler/src'); addpath(fullfile(pwd, 'matlab', 'optiprofiler', 'tests', 'fixtures', 'feature-v2')); gen_b0_fixture('OUT')"
```

The generator prints the `Feature.m` it ran with, which has to be the one of
`/tmp/op-b0-tree`, and stops if a callback records a file. Check the result
before it replaces anything:

```
matlab -batch "a = load('matlab/optiprofiler/tests/fixtures/feature-v2/b0-featured-problem-continuation.mat'); b = load('OUT/b0-featured-problem-continuation.mat'); assert(isequaln(rmfield(a.continuation, 'matlab'), rmfield(b.continuation, 'matlab')))"
```

(`matlab` is the version string of the release that wrote the file; everything
else has to be identical.) Then run
`TestFeatureReviewRegressions/legacyFeaturedProblemRestoresAndContinuesItsStream`
and `TestProblemReference/objectsFromBeforeTheRecordLoadAsUnknown` with the new
file in place.
