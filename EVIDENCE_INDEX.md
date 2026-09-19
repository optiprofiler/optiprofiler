# Evidence index: `codex/reference-facts-scalar`

Tested SHA `249cda5e4d52185f46f7f4fb43c82c658ecb2b4f`, base
`54d42bd550d0493db661711715f868518d433924`. See `IMPLEMENTATION_REPORT.md` for
what the evidence shows; this file says where it is and how to reproduce it.

- `$D` is the audit root on syu-ubuntu,
  `~/audits/op-feature-report-stabilization-20260917`. Older evidence trees in
  it were not modified; everything of this branch has `refscalar` in its name.
- `$H` is the local handoff directory,
  `/tmp/op-feature-report-stabilization-20260917/reference-facts-scalar`. It
  holds a copy of everything below (`evidence/refscalar-evidence.tgz`, unpacked
  next to it) with a SHA-256 manifest. Nothing in `$H` is part of the
  repository.

## 1. Reproduce

```
# on syu-ubuntu; the bundle is incremental on the base, which the clean clone has
cd $D
./refscalar_gates.sh 249cda5e4d52185f46f7f4fb43c82c658ecb2b4f <new-tag> $D/refscalar-249cda5.bundle all
evidence/refscalar/refscalar_mutation.sh 249cda5e4d52185f46f7f4fb43c82c658ecb2b4f
```

Use a new tag for every run: a tag names the log directories and the
disposable working trees, and the helper scripts reuse a working tree that exists.

## 2. Gates at the tested SHA (tag `refscalar3`)

| Evidence | Path under `$D` | Shows |
| --- | --- | --- |
| Driver console | `logs/refscalar-driver-refscalar3.console` | the checked-out SHA and `dirty=0` |
| Bundle verification | `logs/refscalar-bundle-verify-refscalar3.txt` | the bundle is valid and requires only the base |
| Focused summary | `logs/refscalar-focused-refscalar3/summary.txt` | every focused result on one page |
| Exact SHA, clean tree | `logs/refscalar-focused-refscalar3/head_sha`, `dirty_count` | SHA and `0` |
| Changed files | `logs/refscalar-focused-refscalar3/changed_files.txt` | the 25 paths of the branch (15 modified, 10 added) |
| Preserved files | `logs/refscalar-focused-refscalar3/preserved_files_diff.txt` | empty: no gitlink, lock, version or provider change |
| MATLAB focused, JVM | `logs/ml-refscalar3-focused-jvm/console.log`, `junit.xml` | 15 classes |
| MATLAB focused, `-nojvm` | `logs/ml-refscalar3-focused-nojvm/console.log`, `junit.xml` | 6 classes |
| Python focused | `logs/refscalar-focused-refscalar3/pytest_venv{312,311,310,38}.log` | 10 modules per interpreter |
| checkcode | `logs/refscalar-focused-refscalar3/checkcode.tsv`, `checkcode_{head,base,new}_messages.txt` | messages at head and at base per changed file; the new ones |
| Python full suites | `logs/py-refscalar3-{312,311,310,38}/pytest.log`, `junit.xml`, `tree_sha`, `freeze.txt`, `python_version` | one CI-faithful worktree and venv per interpreter |
| Documentation tests | `logs/py-refscalar3-312/doc_tests.log` | `unittest discover -s doc/tests` |
| MATLAB full unit run | `logs/ml-full-refscalar3full/ci.log`, `junit.xml`, `runner.log`, `coverage.xml`, `tree_sha` | `tools/run_matlab_full_unit_ci` |
| MATLAB eval-report gate | `logs/ml-eval-refscalar3/console.log`, `receipt.json`, `junit.xml`, `profile-bands-receipt.json` | `tools/run_matlab_eval_report_ci(..., 'real-workers')` |
| MATLAB ZIP gate | `logs/ml-zip-refscalar3/console.log`, `python_steps.log`, `zip_sha256.txt` | lock check, deterministic ZIP, clean extraction, smoke |

Every `tree_sha` file holds two lines: the SHA of the worktree and the number of
dirty paths. The gate is an exact-SHA gate only if they read
`249cda5e4d52185f46f7f4fb43c82c658ecb2b4f` and `0`.

## 3. Earlier heads of this branch (early warning, not acceptance evidence)

| Head | Tag | Note |
| --- | --- | --- |
| `006766a` | `refscalar1` | focused and all full gates finished. Its focused `summary.txt` says `messages=0` for checkcode; that count is the artifact described in the report, and `checkcode.log` next to it holds the eleven messages |
| `d836c35` | `refscalar2` | focused, four Python suites, eval and ZIP gates finished. Python 3.11 needed one retry after a mirror 403 (`logs/py-refscalar2-311/attempt1-mirror-403/`). The full MATLAB run was stopped when the head moved (`logs/ml-full-refscalar2full/STOPPED_SUPERSEDED.txt`) |

## 4. Mutation checks

| Evidence | Path under `$D/evidence/refscalar` | Shows |
| --- | --- | --- |
| Harness | `mutate.py`, `refscalar_mutation.sh` | the seeded violations and how the tree of the SHA is exported |
| Result at the tested SHA | `mutation-249cda5e4d52185f46f7f4fb43c82c658ecb2b4f/mutation-python.json`, `mutation-matlab.json`, `*.console`, `source_sha`, `source_tree_digest` | control green; every mutation turns the tests red |
| Result at `d836c35` | `mutation-d836c353cec339414a9414967dd18c73fe788a2f/` | the same with the smaller mutation set of that head |

## 5. Probes

All under `$D/evidence/refscalar/probes`.

| Probe | Shows |
| --- | --- |
| `mlprobe/` (`LoadProbe`, `probe.sh`) | R2026a `load` when a set method raises: silent default without `loadobj`, a struct and a warning with `loadobj` |
| `mlprobe2/` (`RefProbe`, `probe2.sh`) | a dependent validated view over a stored record loads without raising, for the current and the superseded layout |
| `make_legacy_fixture.m`, `make_legacy_fixture.sh`, `load_legacy_fixture.m` | how `legacy-point-record.mat` was written with the sources of `0aea026`, and what loading it does |
| `probe_diag.py` | the bounds-versus-linear diagonal-check issue, reproduced and not fixed |
| `probe_compose.py` | all ordered stage pairs build; an empty `custom` stage is legal |
| `smoke.py`, `ml_smoke.m` | first end-to-end checks of both implementations |

## 6. Spelling approximation

`$H/spelling/spell_added_lines.py` (also `$D/evidence/refscalar/`), with the
digests of the dictionaries it used in `$H/spelling/dictionary_digests.sha256`
and its output in `$H/spelling/result_54d42bd_to_249cda5.txt`. It is an offline
approximation of the GitHub workflow, not the workflow.

## 7. Refs and push

The push follows the commit that contains this file, so the last two rows are
written after it.

| Evidence | Path under `$H/refs` | Shows |
| --- | --- | --- |
| Before | `audit_remote_refs_before.txt` | the 21 refs of the audit remote before the push |
| After | `audit_remote_refs_after.txt`, `refs_diff.txt` | the difference between the two listings, which must be the new branch only |
| Receipt | `push_receipt.txt` | the exact push command and its output |
| Final head | `$H/FINAL_HEAD_GATES.md` | the gates run once more at the pushed head |

## 8. Integrity

`$D/evidence/refscalar/MANIFEST.sha256` lists the SHA-256 of every collected
file; `$H/evidence/refscalar-evidence.tgz.sha256` is the digest of the tarball.
The bundle of the tested SHA is `$D/refscalar-249cda5.bundle`, SHA-256
`5c9105a500897a70158eac2c40a0259b770f5e71c53fcd0a66ed3d906aaf7937`.
