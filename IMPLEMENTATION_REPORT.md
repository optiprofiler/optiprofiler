# Implementation report: revised reference contract (`codex/reference-facts-scalar`)

## 1. Identity

- Branch `codex/reference-facts-scalar`, created from
  `54d42bd550d0493db661711715f868518d433924`.
- Tested SHA (all code, tests and documentation):
  `249cda5e4d52185f46f7f4fb43c82c658ecb2b4f`. The last commit on the branch
  adds only this report and `EVIDENCE_INDEX.md`, and excludes these two files
  from the spelling workflow; section 10.4 says how that commit is checked.
- The superseded candidate `codex/reference-facts` stays at
  `0aea026f52b9d96603054178a2044b88acce27b3`, locally and on the audit remote.
  `codex/feature-composition-python` stays at `458ea2d`. No other branch, tag,
  paper reference, version, provider, lock or gitlink was touched: the diff of
  `.gitmodules`, both lock files, `pyproject.toml` and both provider
  directories against the base is empty (`preserved_files_diff.txt`).
- Push target: the configured audit remote only, that is `origin`, which in
  this clone is the local audit repository
  `/tmp/op-feature-opus5-audit-20260915`. GitHub is not a target. The push
  happens after this commit, so its receipt and the listings of the remote refs
  before and after it are in the evidence directory and not in this file.

Commits on the branch:

1. `006766a` `Problem.reference`: four-field feasible reference fact in Python
   and MATLAB (implementation, tests, fixtures).
2. `8440690` Documentation, `REFERENCE_CONTRACT.md` and spelling hygiene.
3. `aa433b7` Validate the record first; say what a reference costs.
4. `d836c35` Constructor docstring.
5. `69eb3c2` Close the registry against subclasses and text subclasses.
6. `249cda5` MATLAB usage page.
7. This report, the evidence index and their two spelling excludes.

## 2. The contract as implemented

`Problem.reference` is optional author/provider metadata with exactly four
fields: a finite real scalar `merit`; a `kind` in `lower_bound`, `optimum`,
`best_known`, `target`, each a claim over the feasible points of the problem; a
non-empty `source`; and a `mapping` token of a closed registry whose only
member is `feasible_objective/1`. There is no reference point, no callable and
no user-defined mapping. All four fields are required on input; an omitted
mapping is not defaulted. `REFERENCE_CONTRACT.md` is the full statement.

## 3. Files

Python
- `python/optiprofiler/opclasses.py`: `ProblemReference` (final, immutable),
  `_reference_merit`, `_reference_text`, `_restore_problem_reference`,
  `_propagated_reference`, `Problem(..., reference=None)` and
  `Problem.reference`, the single-stage propagation in `FeaturedProblem`.
- `python/optiprofiler/feature_definitions.py`: `retains_reference` and
  `REFERENCE_SAFE_CUSTOM_OPTIONS`, next to the other pure stage definitions.
- `python/optiprofiler/composition.py`: every stage view inherits its
  predecessor's record only if its stage retains it; the recorder takes the
  record of the final view.
- `python/optiprofiler/__init__.py`: exports `ProblemReference`.
- `python/optiprofiler/tests/test_problem_reference.py`: new.

MATLAB
- `matlab/optiprofiler/src/Problem.m`: dependent, validated `reference` over
  the stored `reference_`; constructor field; help text.
- `matlab/optiprofiler/src/private/normalizeProblemReference.m`: new validator.
- `matlab/optiprofiler/src/+optiprofiler_internal/propagateProblemReference.m`:
  new retained-or-unknown rule.
- `matlab/optiprofiler/src/FeaturedProblem.m`,
  `matlab/optiprofiler/src/+optiprofiler_internal/FeatureProblemView.m`:
  propagation; the internal scoring method of the views is renamed
  `referenceValue` because `reference` is now the `Problem` property.
- `matlab/optiprofiler/tests/unit_tests/TestProblemReference.m`: new.
- `matlab/optiprofiler/tests/fixtures/reference/`: `legacy-point-record.mat`
  (objects saved by the superseded branch at `0aea026`),
  `RawReferenceProblem.m`, the `reference_invariance` provider, `README.md`.

Documentation
- `REFERENCE_CONTRACT.md` (new), `doc/DEVELOPMENT.md`,
  `doc/source/user/usage_python.rst`, `doc/source/user/usage.rst`,
  `doc/source/user/problem_libraries/custom_python.rst`,
  `doc/source/user/problem_libraries/custom_matlab.rst`,
  `doc/source/matlab/matlab_generated/Problem.rst`,
  `doc/source/dev/problems.rst`.
- `.github/actions/spelling/expect.txt`: the MATLAB API name `intmin`, next to
  the `intmax` it already listed.

## 4. API

Added: Python `optiprofiler.ProblemReference(merit, kind, source, mapping)`
with the read-only properties `merit`, `kind`, `source`, `mapping`, the
informational class attributes `KINDS`, `MAPPINGS`, `FIELDS`, `as_dict()` and
`from_record()`; the trailing optional argument `reference` of `Problem`; the
read-only property `reference` of `Problem` and `FeaturedProblem`. MATLAB: the
optional input struct field `reference` and the property `reference` of
`Problem` and `FeaturedProblem` (`[]` when unknown). The public `benchmark`
signature is unchanged in both languages, and no option reads a reference.

Absent, as required: `point`, `with_point`, `fun`, `maxcv` on the record, the
`affine_inverse` of the MATLAB stage views and `transportProblemReference`.
This branch starts from the base, so these never existed on it; tests assert
their absence so that they cannot come back unnoticed.

Kept: the `reference` property of `Problem` and `FeaturedProblem` and the
load/save plumbing (Python pickling, MATLAB `save`/`load`, provider load).

## 5. Validation and deserialization

Construction fails closed on: a non-finite, non-real or non-scalar merit
(logical values, text and one-element containers included); an integer merit
that a double cannot hold exactly; an unknown or non-text kind; an empty, blank
or non-text source; an unknown, non-text or callable mapping; a missing or
unknown field; a naked scalar; and any record with `fun`, `maxcv` or `point`.
The record is validated before any other constructor argument, so a malformed
record is rejected before any callback of the problem is touched.

| Case | Outcome |
| --- | --- |
| Object serialized before the record existed | unknown, no warning |
| Python pickle of the current layout whose record is rejected | unknown (`None`) with a `RuntimeWarning` |
| Python pickle that restores instance state directly (the superseded layout) | `TypeError` |
| Python instance whose stored value is not a validated record | unknown |
| MATLAB file whose stored record is rejected, including files of the superseded branch | the object loads; `reference` reads as `[]` with `MATLAB:Problem:reference_StoredRecordRejected` |

Nothing is reinterpreted in any of these cases: no field is read as the merit
of another layout and no mapping is guessed.

## 6. Propagation

Retained unchanged or unknown; nothing is transported; the same for every kind.

| Stage | Reference |
| --- | --- |
| `plain`, `noisy`, `truncated`, `random_nan`, `nonquantifiable_constraints`, `unrelaxable_constraints`, `quantized` with `ground_truth=false` | retained (observation only) |
| `perturbed_x0` | retained |
| `permuted`, `linearly_transformed` | retained (coordinate permutation, valid affine change of variables) |
| `custom` with option names within `{mod_x0, mod_affine}` | retained |
| `custom` with any other option | unknown |
| `quantized` with `ground_truth=true` (the default) | unknown |
| any other stage name | unknown |

A composition retains the record only if every stage does. The rule reads
stage names and option names only, so it calls nothing.

## 7. Where the required statements are documented

Feasible reference fact, not a run-history minimum and not the dynamic cohort
minimum; a constrained run merit may be below the reference because of
violation handling and is never clamped; a custom merit function must preserve
the feasible identity before a consumer uses the record:
`ProblemReference` docstring, `Problem` docstring, `Problem.m` help,
`normalizeProblemReference.m` help, `FeaturedProblem.m` help,
`usage_python.rst`, `usage.rst`, `Problem.rst`, `doc/DEVELOPMENT.md` and
`REFERENCE_CONTRACT.md` sections 1.2 and 2.

## 8. Tests

`test_problem_reference.py` has 60 test functions (961 cases);
`TestProblemReference.m` has 37 test methods. Every test carries a comment that
says what it pins and why.

| Requirement | Python | MATLAB |
| --- | --- | --- |
| Four-field shape | `TestRecordShape` (12 functions) | `recordHasExactlyFourFields`, `everyKindIsAcceptedWithTheSameShape`, `omittedReferenceIsUnknown`, `referenceIsSetAtConstructionOnly`, `everyFieldIsRequired`, `unknownFieldsAndNakedValuesAreRejected`, `kindAndSourceAreValidated` |
| Rejection of legacy fields | `test_superseded_fields_are_rejected_not_reinterpreted`, `test_superseded_record_and_transport_api_is_gone` | `supersededFieldsAreRejectedNotReinterpreted`, `supersededRecordAndTransportApiIsGone` |
| Unknown mappings, closed registry | `TestMappingRegistry` (6 functions) | `mappingRegistryIsClosed` |
| Finite validation | `TestFiniteValidation` (6 functions) | `nonFiniteMeritsAreRejected`, `nonScalarAndNonRealMeritsAreRejected`, `integersThatADoubleCannotHoldAreRejected`, `finiteRealScalarsAreStoredExactlyAsDouble`, `validationNeverEvaluatesTheProblem`, `malformedRecordIsRejectedBeforeAnyCallbackIsTouched` |
| Propagation table, single | `TestPropagationRule`, `test_single_stage` (25 stage variants times 4 kinds) | `ruleMatchesTheTable`, `quantizedIsRetainedOnlyForAnExplicitFalse`, `singleStageTableIsUniformOverTheKinds` |
| Propagation table, composed | `test_every_ordered_pair_is_safe_only_if_both_stages_are` (all 576 ordered pairs), `test_longer_compositions` | `everyOrderedPairIsSafeOnlyIfBothStagesAre` (576 pairs), `longerCompositions`, `shorthandOptionsAndPlainTokensReachTheRule` |
| All custom key subsets | `TestCustomKeySubsets` (all 256 subsets: single, first and last in a composition, pure rule) | `allCustomKeySubsets` (256 subsets times three placements), `customWhitelistIsExactlyX0AndAffine` |
| Callback-count invariance | `test_a_reference_costs_no_callback_call` (13 features, user callbacks included) | `aReferenceCostsNoCallbackCall` |
| Safe affine and permutation, test-local optimum without point | `TestSafeCoordinateChanges` (9 pipelines times 3 seeds) | `retainedClaimIsTrueForTheFeaturedProblem`, `invalidAffineNeverYieldsAFeaturedProblem` |
| Default merit feasible identity for `maxcv_init` 0, 0.2, 50, NaN | `TestDefaultMeritFeasibleIdentity` | `defaultMeritIsTheObjectiveAtFeasiblePoints`, `customMeritMustPreserveTheIdentityBeforeAnyConsumerUsesTheReference` |
| Constrained merit below reference, never clamped | `TestConstrainedMeritBelowReference` | `runMeritsFallBelowTheReferenceAndAreNotClamped` |
| Benchmark invariance with and without reference | `TestBenchmarkInvariance` (scores, curves and archived histories, three features) | `benchmarkIsIdenticalWithAndWithoutReferences` |
| Legacy and provider load | `TestLegacyAndProviderLoad` (8 functions) | `saveAndLoadRoundTrip`, `objectsFromBeforeTheRecordLoadAsUnknown`, `supersededLayoutFileLoadsAsUnknown`, `rejectedStoredRecordsReadAsUnknown`, `providerLoadPreservesTheReferenceWithoutExecutingAnything`, `projectX0KeepsTheReference` |

The safe-coordinate test knows the optimum only locally: minimize
`(x1 - 1)^2 + (x2 - 2)^2 + 3` subject to `x1 + x2 <= 2.5` and a box, whose
solution is the projection `(0.75, 1.75)` with value `3.125`. It maps that
point into solver coordinates itself, checks that the featured truth there is
`(3.125, 0)`, and checks on 300 sampled points per case that the featured truth
equals the original objective and violation, so no feasible point is below the
reference. The constrained-merit test uses `minimize 1e7 * x` subject to
`-x <= 0`: a tolerated point has merit `-5e-4` and a penalized point about
`-990`, both below the optimum `0`, and the histories and the archive keep
them as they are.

## 9. Mutation checks

Passing tests show little unless they fail for the right reason. The harness
`evidence/refscalar/mutate.py` exports the tree of the tested SHA, seeds one
contract violation at a time into a fresh copy and runs the reference tests.
A violation that cannot be seeded is reported as a harness error, never as a
detected one. Result at `249cda5`: the unmutated control is
green in both languages (961 Python cases, 37 MATLAB tests), and every one of
the 15 Python and 10 MATLAB mutations turns the tests red, with no harness
error.

| Seeded violation | Python | MATLAB |
| --- | --- | --- |
| `custom` always retained (the superseded rule) | P1 | M1 |
| `quantized` with `ground_truth=true` retained | P2 | M2 |
| an unknown mapping token accepted | P3 | M3 |
| a non-finite merit accepted | P4 | M4 |
| a legacy record reinterpreted (`fun` read as `merit`) | P5 | M5 |
| recorded truth clamped to the reference | P6 | M6 |
| a reference costs one more affine callback call | P7 | M7 |
| a later safe stage restores the record in a composition | P8 | |
| the stored record handed out without validation | | M8 |
| an omitted mapping defaulted | P9 | M9 |
| deserialization maps an unknown token to `feasible_objective/1` | P10 | |
| a callable accepted as mapping | P11 | |
| the record class can be subclassed | P12 | |
| validation reads the assignable class attribute | P13 | |
| text of a `str` subclass trusted as it is | P14 | |
| the record validated last, after the constraint probe | P15 | |
| a sparse merit stored sparse | | M10 |

## 10. Gates on syu-ubuntu

All gates run from the audit root
`~/audits/op-feature-report-stabilization-20260917` through
`refscalar_gates.sh <sha> <tag> <bundle> all`. The driver verifies the bundle,
fetches the branch, requires the bundle head to equal the SHA, checks the SHA
out in the clean clone and then uses the existing isolation wrappers
(`opmatlab`, `remote_mltests.sh`, `remote_full_python.sh`,
`remote_ml_gates.sh`). MATLAB R2026a; Python 3.12, 3.11, 3.10 and 3.8.

### 10.1 Exact SHA and clean tree

The driver printed `HEAD 249cda5e4d52185f46f7f4fb43c82c658ecb2b4f dirty=0`
after the checkout and recorded zero dirty paths again after the focused gates.
Each full gate runs in its own disposable working tree; every `tree_sha` file
(four Python suites, the MATLAB full run, the eval-report gate and the ZIP
gate) holds `249cda5e4d52185f46f7f4fb43c82c658ecb2b4f` and `0`. The diff of the
preserved files against the base is empty, and the branch changes 25 paths (15
modified, 10 added).

### 10.2 Focused gates at `249cda5`

| Gate | Result |
| --- | --- |
| MATLAB focused, JVM (15 classes: `TestProblemReference`, `TestProblem`, `TestFeaturedProblem`, `TestFeatureCompositionV2`, `TestFeatureCustomV2`, `TestFeatureReviewRegressions`, `TestHistoricalNativeFeatures`, `TestProblemLibraryRegistry`, `TestFeatureNativeV2`, `TestConstraintValidity`, `TestFeatureComposedStreamsV2`, `TestTruthRobustness`, `TestConstraintCountInvariant`, `TestMeritFunCompute`, `TestPlainReferenceIdentity`) | Tests=170 Failed=0 Incomplete=0; junit 170 tests, 0 failures, 0 errors; EXIT=0 |
| MATLAB focused, `-nojvm` (6 classes: `TestProblemReference`, `TestProblem`, `TestFeaturedProblem`, `TestFeatureCompositionV2`, `TestFeatureCustomV2`, `TestTruthRobustness`) | Tests=62 Failed=0 Incomplete=0; junit 62 tests, 0 failures, 0 errors; EXIT=0 |
| Python focused, 3.12 / 3.11 / 3.10 / 3.8 (10 modules: `test_problem_reference`, `test_problems`, `test_features`, `test_feature_specification`, `test_feature_composition`, `test_composition_custom_boundary`, `test_quantized_truth`, `test_problem_libraries`, `test_legacy_compat`, `test_profile_utils`) | 1497 passed on each interpreter; EXIT=0 |
| checkcode of the 9 changed `.m` files, at head and at base | 9 of 9 analyzed; 4 messages at head, the same 4 at base (`Problem.m`, `project_x0`, untouched); 0 new |

### 10.3 Full gates at `249cda5`

| Gate | Result |
| --- | --- |
| Python full suite, 3.12.14 | 2855 passed, 3 skipped; junit 2858 tests, 0 failures, 0 errors; EXIT=0 |
| Python full suite, 3.11.16 | 2855 passed, 3 skipped; junit 2858 tests, 0 failures, 0 errors; EXIT=0 |
| Python full suite, 3.10.12 | 2855 passed, 3 skipped; junit 2858 tests, 0 failures, 0 errors; EXIT=0 |
| Python full suite, 3.8.20 | 2855 passed, 3 skipped; junit 2858 tests, 0 failures, 0 errors; EXIT=0 |
| Documentation tests, 3.12 | Ran 14 tests, OK; EXIT=0 |
| MATLAB full unit run (`tools/run_matlab_full_unit_ci`, pool capped at 4) | Tests=498 Failed=0 Incomplete=2; junit 498 tests, 0 failures, 0 errors, 2 skipped; "Full unit suite completed without failed tests"; EXIT=0 (10:44:06Z to 11:38:22Z) |
| MATLAB eval-report gate (`real-workers`) | `TestEvalReport` 12 of 12 cases, Failed=0, Incomplete=0; `TestProfileBands` 6 of 6, Failed=0, Incomplete=0; EXIT=0 |
| MATLAB ZIP gate | lock validated; ZIP validated with 1279 files; `SMOKE-PASSED`; MATLAB-EXIT=0; user startup and pathdef unchanged |

The suite grew from 1894 to 2855 passing Python cases, which is the 961 cases
of the new module, and from 461 to 498 MATLAB tests, which is the 37 tests of
the new class. The three Python skips are the suite's usual ones (LaTeX not
installed, a Windows-only junction case, one case that needs an environment
variable). The two incomplete MATLAB tests are the Windows-only
`TestSetupPathOwnership/testWindowsRelativePathsRejected` and
`testWindowsAbsolutePathForms`, filtered by their platform assumption on Linux.
No gate at `249cda5` was retried.

### 10.4 The report commit

A commit cannot contain the results of its own test run. The gates above ran at
`249cda5`, before this report was written. The last commit adds only
`IMPLEMENTATION_REPORT.md` and `EVIDENCE_INDEX.md`, which no test, package or
archive reads, and two lines in `.github/actions/spelling/excludes.txt` that
keep these two files out of the spelling workflow (they list tool names and log
paths, not prose). The same driver is run once more at the pushed head after this
commit. Those receipts cannot be part of the commit; they are written to
`FINAL_HEAD_GATES.md` in the evidence directory named in `EVIDENCE_INDEX.md`.

### 10.5 Earlier heads

The full gates also ran at two earlier heads of this branch, as an early
warning and not as acceptance evidence: `006766a` (the implementation commit)
and `d836c35`. At `006766a` everything finished green: Python
2852 passed and 3 skipped on each interpreter, MATLAB Tests=497 Failed=0
Incomplete=2, eval-report and ZIP gates EXIT=0. At `d836c35` the focused gates,
the four Python suites (2853 passed, 3 skipped each), the eval-report gate and
the ZIP gate finished green; its full MATLAB run was stopped when the head moved
to `249cda5`, with no failure logged in the 346 tests it had completed.

## 11. Findings that a reviewer should know

1. **A checkcode count of zero was an artifact.** The first version of the
   checkcode step counted console lines with a pattern that never matched
   (the lines start with one blank, the pattern expected two), so it printed
   `messages=0` while the log held eleven messages. The driver of the
   superseded branch has the same pattern: its report states zero messages,
   but its log (`logs/reffacts-focused-reffacts2/checkcode.log`) holds four.
   The step now lets MATLAB write a TSV, analyzes the base version of every
   changed file as well, and reports head, base and new messages. Seven
   messages of the new test class were real and are fixed. The remaining four
   are in `Problem.m` in `project_x0`, which this branch does not touch; they
   are identical at the base.
2. **"Executes nothing" was wrong for Python.** The Python `Problem`
   constructor probes `cub` and `ceq` at `x0` for their dimensions, with or
   without a reference. The accurate statement is that a reference adds no
   callback call. The wording is corrected everywhere, the record is now
   validated before any other argument in both languages, and tests pin that a
   malformed record is rejected before any callback runs.
3. **The closed registry had a hole in Python.** Validation read
   `self.MAPPINGS`, so a subclass could bring its own registry, and a `str`
   subclass with an overridden `__eq__` could pass the membership test. The
   class is now final, validation reads module-private literals, and text is
   copied to exact `str` before it is judged.
4. **MATLAB load semantics decided the storage design.** On R2026a, when a
   property set method raises during `load`, a class without `loadobj` keeps
   the default silently, and a class with `loadobj` (`FeaturedProblem`)
   receives a struct and MATLAB warns that the constructor must preserve the
   class. A dependent property without a set method behaves the same way for a
   value saved under the superseded layout. Hence `reference` is a dependent,
   validated view whose protected set method stores raw: loading never raises,
   and the getter decides what a caller sees. The probes are in
   `evidence/refscalar/probes/`.
5. **One infrastructure retry.** At `d836c35` the Python 3.11 worker failed
   before any test ran: the PyPI mirror answered 403 while four installs
   resolved at once. The failed attempt is kept in
   `logs/py-refscalar2-311/attempt1-mirror-403/`; the retry passed. The driver
   now staggers the installs, and the run at `249cda5` needed no retry.

## 12. Out of scope

- No consumer: no benchmark option, score, profile formula, archive or report
  reads the reference.
- No offline catalog, scoring manifest or independent scorer.
- **The bounds-versus-linear diagonal-check tolerance issue is recorded, not
  fixed.** For `linearly_transformed` and for `custom` with `mod_affine`, the
  bounds modifier tests the *inverse* for exact diagonality (MATLAB
  `isdiag(inv)`) while the linear modifiers test the *matrix* (`isdiag(A)`).
  With `A = diag(2, 4)` and an inverse carrying off-diagonal entries of size
  `1e-17`, the bounds become infinite and no bound rows are added, so the
  bounds vanish from the structure handed to the solver
  (`evidence/refscalar/probes/probe_diag.py`). The scoring truth still measures
  the violation with the original problem, so the reference rule does not
  depend on this defect. It needs its own fix and tests.
  *Note added 2026-09-20: fixed by the affine safeguard commit that follows
  this report on the branch, and hardened by the follow-up commit after it
  (exact structural decision, one transformation per problem, two-sided
  inverse, overflow guards); see `REFERENCE_CONTRACT.md`, section 7.*

## 13. Platform limitations

- The local sandbox lacks `pandas` and `pypdf` and has no provider checkout;
  all evidence comes from syu-ubuntu (Linux, MATLAB R2026a only). macOS,
  Windows and older MATLAB releases are not covered by these runs.
- syu-ubuntu was shared and saturated during the runs (load average about 125
  on 120 cores, mostly another user's jobs), so the durations are not
  representative.
- The Sphinx documentation was not built: syu-ubuntu has no Sphinx toolchain and
  a documentation build is not one of the gates. The six changed `.rst` pages
  were parsed with docutils, with the Sphinx roles and directives replaced by
  inert stand-ins: no structural warning. This does not check cross-references
  or the numpydoc rendering of the new class page.
- The check-spelling workflow runs on GitHub only. It was approximated offline
  with the workflow's own dictionaries (`spell_added_lines.py`): no unknown
  word and no forbidden pattern in the lines this branch adds. This is an
  approximation and not the workflow.
- The known harness quirks apply unchanged: `qpdf`/`gs` wrappers on the PDF
  tool path, the ZIP smoke session executes the user's `startup.m` before
  `restoredefaultpath` (its checksum is verified unchanged), and the default
  parallel pool is capped at four workers.

## 14. Judgment calls to confirm

1. All four fields are required on input, `mapping` included. A default would
   be friendlier, but an omitted mapping would then be interpreted.
2. A rejected record of the current pickle layout reads as unknown with a
   warning, while a stream of the superseded layout raises. Both satisfy
   "unknown or error, never reinterpret"; the first keeps a problem written by
   a later version loadable.
3. In MATLAB a rejected stored record warns on every read rather than once at
   load, because `Problem` has no `loadobj` and adding one would change how
   unrelated load failures surface.
4. An integer merit that a double cannot hold exactly is rejected rather than
   rounded.
5. `intmin` was added to the check-spelling expect list; coined fixture names
   were replaced by dictionary words instead of being added. This report and
   the evidence index are excluded from the spelling workflow instead of
   adding the names of audit tools and log directories to the expect list. If
   the two documents are not wanted on the development line, they and the two
   exclude lines can be dropped together.
