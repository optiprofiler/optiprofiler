# Development Plan

This document records the development roadmap, not a list of released APIs.
The provider split is implemented on the development line. The running-record,
feature-composition, and reference/scoring work below is still planned; names
and file formats are provisional until implementation and consumer tests agree.
The paper maintenance line receives applicable fixes, not these new features.

## Separate Problem Libraries from the Engine

Keep S2MPJ as the bundled default problem library so that `benchmark(solvers)`
continues to work out of the box. Optional problem libraries are separated
from the OptiProfiler engine and loaded through explicit library contracts.

For Python, the implemented design is a plugin-style protocol: external packages
such as `optiprofiler-solar` can expose problem-library metadata and `select`
/ `load` entry points, while OptiProfiler discovers and validates them without
hard-coding their implementation details.

For MATLAB, the implemented design uses an explicit registry and setup contract
for optional problem libraries. The engine depends on the agreed
library interface, not on optional library source trees being present inside the
package.

The remaining distribution work is to verify compatible released package
metadata, normal dependency resolution in a fresh environment, and the actual
wheel, sdist, and MATLAB ZIP contents. The repository locks record tested
integration revisions; they do not prove that an arbitrary installed provider
contains that revision. Reproduction instructions must record the exact core,
provider, and separately managed runtime versions without forcing every daily
installation to use one immutable integration set.

## Problem-Library Lifecycle and Diagnostics

The implemented lifecycle lets MATLAB users unregister an optional
or custom provider without deleting its source, runtime, cache, or other user
data. Core uninstall and provider removal have separate, explicit
semantics.

A later release should add richer, outcome-equivalent inspection in both
languages. Python may extend its entry-point-based listing, while MATLAB may
inspect its registry, but both should report comparable fields: public name,
role, API version, source and installed version or locked commit, location,
runtime availability, cache or build location, duplicate or stale state, and a
suggested repair. This future work may expose `list`, `status`, and `doctor`
operations; it is not part of the current lifecycle closeout.

A future destructive MATLAB removal API must remain separate from `setup`
and `unregisterProblemLibrary`. It should operate only on files recorded in a
setup-generated ownership manifest, show a dry-run removal plan, and treat the
adapter checkout, runtime, and cache as separate scopes. It must refuse bundled
providers, custom or unknown roots, and dirty checkouts; it must not expose a
blanket force option that recursively deletes an arbitrary registered root.

## MATLAB Release Artifact (Required for the Next Release)

CI already builds and tests a MATLAB-only ZIP. The next release should publish
the accepted archive from the same tagged monorepo commit. The archive should
contain `setup.m`, the MATLAB engine, the MATLAB problem-library lock,
licenses and required metadata,
and the bundled S2MPJ contents, without Python sources, documentation build
inputs, workflows, or Git history. A MATLAB toolbox (`.mltbx`) may be added in
the same or a later release.

The monorepo remains the single source of truth. Do not maintain a second
hand-edited MATLAB core repository; if a Git-based distribution mirror becomes
necessary, generate it automatically from the release commit.

## Machine-Readable Experiment Records and Evolve Feedback

Start with a small, provisional `run_record` consumed by a real evaluator.
Core records experiment facts: run identity, actual problem/provider identity,
feature and randomness policies, budgets, solver identities, status and errors,
score definitions, and references to saved numerical data with axes and units.
Keep numerical arrays in the existing data files; JSON uses standard types and
an explicit representation for unavailable or non-finite values.

Distinguish experiment completion from solver failures and rendering failures.
Publish required data before a completed record; a surviving `running` record
alone does not prove that its owner has stopped. Preserve existing `score_only`
and legacy load semantics. Decide any new output mode from consumer needs rather
than adding a broad experiment-configuration API in advance.

Evolve owns compact, stage-specific feedback for its agents, including permitted
scores, error categories, history/profile summaries, and hidden-set isolation.
Core must not need an `agent_mode` or know whether the consumer is an agent,
CLI, or web service. The withdrawn `resolve_problem_selection` / `agent_report`
candidate is not the contract for this work.

## Ordered Feature Composition

Design `FeaturePipeline` as an ordered list of stages with explicit identities
and per-stage options. A spelling such as `noisy+perturbed_x0` may be shorthand,
not the canonical saved identity. Compile the stages into one FeaturedProblem;
do not recursively wrap counters, budgets, truth evaluation, or histories.

Start with a small supported whitelist: `noisy` with `perturbed_x0`, followed
by order-sensitive `noisy` / `truncated` combinations. Distinguish conflicting
combinations from valid but not yet supported ones. Spatial transformations
require equivalent bounds, constraints, truth, and derivative semantics before
their combinations are enabled. Python and MATLAB must agree on the contract
and each support language-local replay; matching random samples across the two
languages is not promised.

## Reference Facts and an Independent Arena Scorer

Begin with an offline reference catalog and a frozen scoring manifest, not a
new public scalar `Problem` property. Distinguish a rigorous lower bound, exact
optimum, best-known feasible value, and chosen scoring target. Record provenance
and applicability; classify how each feature preserves, transforms, or invalidates
those facts. Missing reference facts remain unknown, not a claimed lower bound
at the initial point.

The manifest fixes the target, per-problem feasibility tolerance, budget axis,
weights, and policy version. Feasibility-only tasks need a separate interpretation.
Provider updates must not silently change a signed-off scoring standard.
Validate with a small, identifiable problem set and actual saved truth histories.

Ordinary profiles remain relative comparisons. Do not use the minimum of a
reference and the current solver cohort's best value as a fixed Arena target:
the result would still depend on the competing solvers. The independent scorer
must pass this test: hold one solver's history fixed, add/remove/reorder other
solvers, and its absolute score stays unchanged. Legacy data need verifiable
provenance before Arena use; ordinary trusted-data reload remains supported.

## Separate Follow-Ups and Release Boundaries

- Keep common-oracle noise coupling as a separate experimental policy decision;
  do not silently change existing single-feature seeds or run indexing.
- Incremental checkpoints/restarts, the MATLAB constraint-count callback
  contract, and richer provider diagnostics are separate bounded designs.
  Final atomic saving is not a solver checkpoint or a power-loss guarantee.
- Public upload services must not treat current pickle-bearing H5/options as
  safe numerical input. Trusted-archive warnings are already documented; safe
  formats and isolation are still service prerequisites.
- Before releasing, verify the tag against the artifact version, distinguish
  paper/Python/MATLAB publication, and inspect final licenses, notices, and
  reproducibility metadata. S2MPJ's existing BSD-3-Clause text must be retained;
  it is no longer an upstream-missing-license item.

Implement in small accepted slices: close remaining documentation/provider
issues; prepare installation/release checks alongside the run-record design;
validate the record with Evolve; add feature combinations; then run the offline
reference/scoring pilot. A future 2.0.0 release need not wait for every research
item. Do not publish or change the version as a side effect of this roadmap.
