# Development Plan

This document records the development roadmap, not a list of released APIs.
Implemented on this development line: the provider split, the 2.0 `Feature`
contract with ordered feature composition, experiment plans and provenance,
the trusted boundary for earlier configurations, version 2 of the
machine-readable evaluation report, and the optional feasible reference fact
on `Problem` and `FeaturedProblem`. The offline reference catalog and Arena
scorer remain future work; storing a reference does not change profile scores.
These are development features, not a published 2.0.0 release. The paper
maintenance line receives applicable fixes, not these new features.

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
(Python `seedsequence-v2`; MATLAB `matlab-stage-horner32-v2`, an exact 32-bit
Horner fold over the stage identity and, per query, over the IEEE-754 words of
the observed payload, feeding the unchanged legacy stream kernel; version 1
handed the payload to the legacy product mixer, which collapsed whenever any
payload element was zero; neither language promises the other's samples). One outer recorder owns solver budgets, public
call counters and histories; stage views own local served-query state and
transform observation and reference paths. Public recorders are not nested.
Spatial views transport bounds and constraints; composed derivative methods
explicitly reject unsupported calls (Python raises `NotImplementedError`,
MATLAB errors).

Provenance is explicit (`feature_pipeline-v3`: feature block plus experiment
block; `options_refined-v2` for native replay). Source historical archives are
never modified and their payloads are read without upgrading; the metadata a
new report retains from a loaded archive is a sanitized, bounded copy, never a
verbatim copy. Historical serialized configurations enter only through the
trusted compatibility boundary (`optiprofiler.legacy_compat`). Python and
MATLAB agree on the contract and each support language-local replay and seed
policies; matching random samples across the two languages is not promised.

## Reference Facts and an Independent Arena Scorer

Build the scorer around an offline reference catalog and a frozen scoring
manifest, not a naked scalar `Problem` property. Distinguish a rigorous lower bound, exact
optimum, best-known feasible value, and chosen scoring target. Record provenance
and applicability; classify how each feature preserves, transforms, or invalidates
those facts. Missing reference facts remain unknown, not a claimed lower bound
at the initial point.

The in-memory carrier of one such fact is implemented on this development
line (see `REFERENCE_CONTRACT.md`; it supersedes the point-carrying record of
the earlier `codex/reference-facts` candidate). `Problem.reference` is
an optional record of exactly four fields: a finite real scalar `merit`, a
`kind` (`lower_bound`, `optimum`, `best_known`, `target`; every kind is a claim
over feasible points), a non-empty `source` and a `mapping` token of a closed
registry whose only member is `feasible_objective/1`. It is still not a public
scalar property: a naked scalar is rejected, and the scalar is never read
without its kind, provenance and mapping. There is no reference point, no
callback and no user-defined mapping. The record is validated fail-closed
without evaluating the objective; a serialized record that the running version
does not accept reads as unknown or raises and is never reinterpreted.
`FeaturedProblem` either retains the record unchanged or reports it as unknown:
observation-only stages, `perturbed_x0`, `permuted` and `linearly_transformed`
retain it, `custom` only when its options are a subset of `mod_x0` and
`mod_affine`, `quantized` with `ground_truth=true` never, and a composition only
if every stage does; the rule is the same for every kind.

The fact is a feasible reference, not a run-history minimum and not the dynamic
cohort minimum of the profiles, and nothing derives it from solver output. It
is not a floor either: on a constrained problem a run merit may be below it,
because the merit function tolerates or penalizes small violations, and such
merits are never clamped. No profile formula, archive or report reads it yet.
A future consumer must first establish that the merit function in use
preserves the feasible identity `merit_fun(f, 0, maxcv_init) == f` (the default
one does, for every `maxcv_init` including NaN); otherwise run merits and the
record are not in the same space. The offline catalog, the scoring manifest and
the independent scorer remain future work; the catalog can fill the record at
load time through a provider.

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

The provider split, evaluation report `optiprofiler.eval_report/2` with the
immutable version 1 and the `plot_data/1` companion, ordered feature composition,
experiment plans and provenance, trusted historical-configuration boundary,
and feasible reference fact are implemented here. Still future work:
installation/release checks, validation of the current report contract by
Evolve as an external integration gate, and the offline reference/scoring
pilot. A future 2.0.0 release need not wait for every research item. Do not
publish or change the version as a side effect of this roadmap.
