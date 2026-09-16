# Frozen 1.x FeaturedProblem fixture

`b0-featured-problem.mat` holds two live `FeaturedProblem` objects saved by the
MATLAB implementation at source `1c4b7d9170ad40e1112e33e5a5638cecc9e644eb`
(the `main` baseline before Feature composition and EvalReport), MATLAB R2026a
on syu-ubuntu, 2026-09-17:

- `fp_noisy`: `FeaturedProblem(problem, Feature('noisy', 'noise_level', 0.5,
  'n_runs', 2), 6, 17)` after two objective evaluations (`[1; 2]`, `[0.5; 1]`)
  and one inequality-constraint evaluation (`[1; 2]`);
- `fp_plain`: `FeaturedProblem(problem, Feature('plain'), 6, 17)` after one
  objective evaluation (`[1; 2]`);

where `problem = Problem(struct('fun', @(x) sum(x.^2), 'x0', [1; 2], 'cub',
@(x) x(1) - 0.5, 'name', 'b0fixture'))`. Their `feature` properties are the
1.x name/options layout, which the current `Feature.loadobj` returns as a
`LegacyFeatureEnvelope`; `FeaturedProblem.loadobj` rebuilds the single-feature
runtime through `importLegacyFeature`.

`b0-featured-problem-continuation.mat` records what the same live objects
returned next under the 1.x implementation (`fun([0.25; 0.5])`,
`cub([0.25; 0.5])`, `maxcv([0.25; 0.5])`, the histories and counters), which a
restored object must reproduce.

SHA-256:

- `b0-featured-problem.mat`
  `5bfb0da47db2af835e351aad12d5510bb5a44d68de92aaa07b9b8c285a532513`
- `b0-featured-problem-continuation.mat`
  `ab677e7560098a78464cde3976aee5e7eca79c6c22cb17409518861d7c49aa6a`

Generator: the script `gen_b0_fixture.m` recorded in the stabilization
evidence (`~/audits/op-feature-report-stabilization-20260917/evidence` on
syu-ubuntu), run against the `1c4b7d9` tree with an isolated path. Native MAT
loading is a trusted-input boundary.
