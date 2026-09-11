# Development Plan

This document records the development roadmap, not a list of released APIs.
The provider split, the machine-readable evaluation report and the ordered
feature composition are implemented on the development line (their current
contracts are documented in the user guide); the reference/scoring work below
is still planned, and its names and file formats are provisional until
implementation and consumer tests agree. The paper maintenance line receives
applicable fixes, not these new features.

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

The record is the opt-in `eval_report` (`report_path=`): a versioned main
report (`optiprofiler.eval_report/2`, emitted by both languages; version 1
stays immutable as the contract of reports written before it) with a numeric companion
(`optiprofiler.plot_data/1`), both pinned by packaged JSON Schemas and selected
by the document's schema identifier. Core records experiment facts: run
identity, actual problem/provider identity, the canonical feature specification
and the experiment plans (run counts by role), randomness policies, budgets,
solver identities, status and errors, score definitions, and references to
saved numerical data with axes and units. Numerical arrays stay in the existing
data files; JSON uses standard types and an explicit representation for
unavailable or non-finite values.

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

`Feature` is the canonical pipeline specification: an ordered list of stages
with explicit identities (`name#occurrence`, frozen literal seed codes) and
stage-local options, built once from the `feature_name` shorthand or the
structured `feature` entries and never reparsed. The identity pipeline has
zero effective stages and the name `plain`; empty or malformed input is an
error. Experiment settings such as `n_runs` are not feature options: the
experiment layer resolves one plan per role (primary, plain reference) from the
specification and the solver metadata. Every built-in stage composes with every
other in any order and any number of times; a single effective stage keeps the
established single-feature execution and seeds, and a genuine composition is
executed by one recorder over lazily composed views with per-stage, per-channel
seeds derived by a language-local policy over the same stage identities
(Python `seedsequence-v2`; MATLAB `matlab-stage-horner32-v1`, an exact 32-bit
Horner fold feeding the unchanged legacy stream kernel; neither language
promises the other's samples). Counters, budgets, truth evaluation and histories
live in the recorder only; stages are never wrapped recursively. Spatial
transformations carry bounds, constraints, truth and derivative semantics
through the same view protocol.

Provenance is explicit (`feature_pipeline-v3`: feature block plus experiment
block; `options_refined-v2` for native replay), historical payloads are kept
verbatim, and historical serialized configurations enter only through the
trusted compatibility boundary (`optiprofiler.legacy_compat`). Python and
MATLAB agree on the contract and each support language-local replay and seed
policies; matching random samples across the two languages is not promised.

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

Implemented in small accepted slices on isolated candidate branches (the
provider split is integrated on the development line; the versioned
evaluation report `optiprofiler.eval_report/2` with the immutable version 1
and the `plot_data/1` companion, the ordered feature composition with the 2.0
`Feature` contract, the experiment plans and provenance, and the trusted
boundary for configurations written by earlier versions are candidate
implementations awaiting integration into the maintained branch after
independent review). Still future work: installation/release checks, validation of the
record by Evolve as an external integration gate, and the offline
reference/scoring pilot. A future 2.0.0 release need not wait for every research
item. Do not publish or change the version as a side effect of this roadmap.
