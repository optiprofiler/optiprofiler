.. _eval_report:

Machine-readable evaluation reports
===================================

``eval_report`` is an opt-in, versioned JSON report for controllers, experiment
tools, and applications that need more than a scalar solver score. It records
measurements and their provenance; it does not change how OptiProfiler solves,
scores, or draws an experiment. It is not an agent execution mode.

Python
------

Keep the ordinary three-output API and request a new report file::

    scores, profile_scores, curves = benchmark(
        [solver_a, solver_b],
        plibs=['s2mpj'], ptype='u', mindim=1, maxdim=10,
        score_only=True,
        report_path='evaluation/eval_report.json',
    )

``report_path`` accepts a string or ``pathlib.Path``. With no report path (or
``None`` or an empty string), the previous behavior is retained: the collector
module is not even imported, no metadata is attached to solver results, and no
extra file is touched. That opt-out is the default and costs nothing.
An explicit path requests two JSON files: the small main report and one
adjacent numeric companion. For ``evaluation/eval_report.json``, the companion
is ``evaluation/eval_report.plot_data.json``. Both targets must be new: a
collision with either file is refused before solvers execute. This also
prevents a load/replot operation from replacing its source report.

``score_only=True`` still disables plots and raw experiment persistence. The
two explicitly requested JSON files are the only additional outputs. Preparing
plot coordinates does not create figures, require a rendering backend, or
certify that a PDF exists. Numeric presentations and binned histories do not
replace the complete H5/MAT archive.

The report can be read with the standard library::

    import json
    with open('evaluation/eval_report.json', encoding='utf-8') as stream:
        report = json.load(stream)
    assert report['schema'] == 'optiprofiler.eval_report/1'
    budget = report['problems'][0]['budget']['evaluations']
    for run in report['problems'][0]['runs']:
        print(run['solver_index'], run['run_index'], run['evaluations'],
              run['budget_reached'], run['objective']['best'])

MATLAB
------

The equivalent options-struct field preserves MATLAB's three outputs::

    options = struct('plibs', {{'s2mpj'}}, 'ptype', 'u', ...
        'mindim', 1, 'maxdim', 10, 'score_only', true);
    options.report_path = 'evaluation/eval_report.json';
    [scores, profile_scores, curves] = benchmark({@solver_a, @solver_b}, options);

The MATLAB implementation emits the same schema, the same field names and
the same one-based solver/run indices as Python. Runtime provenance,
random-seed rules and error-band normalization remain language-specific; the
report states which convention produced a value instead of changing either
implementation to make them match (see :ref:`eval_report_conventions`).

Report structure
----------------

The :download:`main-report schema <../../../python/optiprofiler/schemas/eval_report.schema.json>`
and :download:`numeric-companion schema <../../../python/optiprofiler/schemas/plot_data.schema.json>`
are the single Python/MATLAB reader contract. Observation records (problems,
runs, metrics, coverage, stages, plot records, history channels) allow only the
declared keys, so a consumer never has to branch on the producing language
for the same concept. The two files are package resources of the Python
distribution (``optiprofiler/schemas/``) and can be read at run time::

    from optiprofiler.eval_report import load_schema, schema_text
    main_schema = load_schema('eval_report')      # parsed dict
    companion_text = schema_text('plot_data')     # exact JSON text, for hashing/pinning

There is exactly one authoritative copy of each schema: the documentation
downloads above, the installed Python test suite, the MATLAB test fixtures and
external consumers all read the same files, so a consumer that pins a schema by
its SHA256 compares against the packaged text. They require no additional
runtime dependency in OptiProfiler; both test suites validate every produced
report against them with a small built-in checker. Consumers should reject bare
JSON ``NaN``/``Infinity`` tokens, validate the schemas, and check tensor
dimensions, references and artifact paths. This remains an unreleased v1 report
format, independent of package versions.

.. list-table:: Main fields
   :header-rows: 1
   :widths: 25 75

   * - Field
     - Meaning
   * - ``configuration`` / ``source``
     - ``configuration.request`` is what the caller supplied;
       ``configuration.effective`` is the resolved problem/profile options and
       feature, stated once. For a load operation these describe reanalysis
       and rendering, and ``source`` is the original archive receipt.
   * - ``stages`` / ``coverage``
     - Independent numerical, scoring, persistence and rendering states;
       selected, loaded and completed primary problems.
   * - ``problems[]``
     - Identity (``[library, name, role]``), dimension, type, selection and load
       status, the shared evaluation ``budget`` of the problem, the provider
       (when known) and the runs.
   * - ``problems[].runs[]``
     - Per-run facts only: evaluation count, ``budget_reached``, returned,
       initial and best objective/constraint/merit values, invalid-value
       counts, abnormal-termination and fallback flags, actual-versus-repeated
       execution, oracle seed and elapsed time.
   * - ``scores`` / ``profiles``
     - Exact returned scores with named axes, actual convergence-work
       summaries, and references to numeric profile presentations.
   * - ``semantics`` / ``plot_data``
     - Shared meanings stated once, and the companion's path, identity-linked
       contents and integrity receipt. Runs refer to its histories by ID.
   * - ``artifacts`` / ``diagnostics``
     - Existing output files with relative paths and hashes, and bounded
       structured event codes. Files marked ``present`` are not automatically
       certified scientifically correct or safe to publish.
   * - ``report_files``
     - The best-effort file-permission policy applied to the two JSON files
       and whether this platform applied it.

Solver/run/sample indices in observation fields use one-based indexing,
including Python JSON. Options recorded in ``configuration`` retain the
original language API conventions; for example, Python ``solvers_to_load``
still uses zero-based indices.
The profile-score axes are solver, tolerance, history/output, and profile type.
Objective, constraint and merit ``best`` values are componentwise minima: they
need not belong to the same evaluated point. Work summaries reuse the existing
profile convergence calculation; they do not introduce a new stopping rule.

Meaningful absence and meaningful nulls
---------------------------------------

The main report keeps invariant meanings once in ``semantics`` rather than
repeating the same prose in every run, and it omits keys that would only carry
``null``. The rules are stated in ``semantics.run_defaults``:

* ``problems[].budget`` holds the evaluation cap shared by every run of the
  problem in this invocation (``ceil(max_eval_factor * dimension)``). A run
  record carries its own ``budget`` only as an explicit override. In a load
  operation the original budget is not retained, so ``budget.evaluations`` is
  ``null`` with a reason and every ``budget_reached`` is ``null``; today's
  ``max_eval_factor`` is never used to guess it.
* ``budget_reached``, ``evaluations``, ``abnormal_termination``,
  ``output_fallback``, ``execution``, ``oracle_seed`` and ``elapsed_seconds``
  stay per run. Reaching the cap still does not prove the termination cause.
* There is no per-run convergence field. Convergence is never inferred from a
  solver return value; ``target_work`` in the companion observes the existing
  scoring predicate.
* ``invalid_evaluations`` is the integer ``0`` when every retained evaluation
  was finite; otherwise it is the categorized count object, and
  ``first_invalid_evaluation_index`` is present. An absent
  ``availability_reason`` means the observation was available; when a history,
  count, returned or initial value is missing, the reason names what is
  missing (for example ``history_or_evaluation_count_unavailable`` for a merit
  history that a load with a new merit callback does not retain).
* ``oracle_seed_reason`` and ``termination_metadata_reason`` appear only when
  the corresponding value could not be observed (old archives).

Absence or ``null`` must never be read as a successful zero measurement.

Three different kinds of data
-----------------------------

The main report is a concise index of facts. Its numeric companion contains
``histories``, ``plots`` and ``target_work``. The original H5/MAT archive is
unchanged and remains the only complete source for reanalysis.

* **History summaries** describe each problem/solver/run and each metric
  separately, excluding padded values beyond actual evaluations. Short
  histories of up to 64 observations use ``representation="exact_samples"``
  and retain all values. Longer histories use ``representation="bin_extrema"``
  and are **lossy**: at most 32 contiguous bins retain their first/last
  observations, finite minimum/maximum with true evaluation positions, and
  NaN and positive/negative infinity counts. Bin ``k`` of ``n`` retained
  evaluations covers evaluations ``(k-1)*n//32+1`` to ``k*n//32`` in both
  languages. Ties use the first occurrence. The order of observations inside
  a bin is not retained, and bins are never scoring inputs.
* **Numeric plot presentations** contain the actual prepared coordinates,
  central curves and error bands from the same numeric preparation used by
  the renderer. ``exact_rendered_data`` means exact *prepared display data*,
  not a lossless copy of the raw experiment, pixel-identical rendering across
  machines, or evidence of a successful PDF. History display clipping
  (``display_limit``), nonfinite placeholders (``nonfinite_policy``), shifts,
  raw/cumulative views and block aggregation are recorded on each panel.
  ``observation_scope="rendering_inputs"`` identifies observations captured
  from the actual rendering path. When drawing is not requested, a
  representation can instead use ``retained_scoring_observations``. These
  scopes are not interchangeable: a stateful custom merit callback can produce
  different values when the existing renderer calls it again. Reporting does
  not invoke that callback an extra time to make the two agree.
* **Target work** retains the existing profile calculation's work arrays and
  their problem/solver/run/tolerance identities. History-based work is the
  first evaluation meeting the existing target. Output-based work is the
  evaluation count associated with a passing returned output, not the first
  history hit. Never infer either value from a compressed curve or from a
  history preview.

Padded display aggregation
--------------------------

Both renderers keep at most about 1000 interior points of a history plot. The
rule is keyed on the *padded* history length ``ceil(max_eval_factor *
dimension)``, recorded as ``padded_length`` on every history panel, not on the
number of evaluations a solver actually used. Above 1002 the interior
evaluations are grouped into blocks and one representative per block is kept
(``aggregation`` is ``min``, ``mean`` or ``max``, as configured by
``hist_aggregation``); the first and the last actual evaluation are always
kept, and the last one can therefore appear twice. Consequently an array
position in ``series[].mean`` is **not** an evaluation number once
``padded_length`` exceeds 1002. Locate a sample through
``series[].evaluation_indices`` and read the value at that position. The
lossy bins and the scalar facts in the main report are independent of this
display aggregation: an isolated spike that a block drops from the display
copy is still visible as a bin extremum and in ``best``/``first_invalid``
indices.

The scalar ``best`` ignores NaN and can be infinite. A bin's ``finite_min`` and
``finite_max`` intentionally exclude infinities; they answer a different
question. Exact scalar facts, abnormal/fallback flags and invalid counts remain
available even if a short visual summary hides an event.

.. _eval_report_conventions:

Language-specific conventions
-----------------------------

The vocabulary is shared, the science is not equalized. The report says which
convention produced a value:

* **Error bands**: Python retains its population standard deviation
  (``std_ddof = 0``); MATLAB retains its sample normalization (``std_ddof = 1``
  for more than one run). ``semantics.error_bands`` in the companion states
  the emitting language's rule.
* **Oracle seeds**: each run's ``oracle_seed`` is the seed actually given to
  the featured problem under the emitting language's own rule (Python
  ``(23333*seed + 211*i) % 2**32`` with 0-based ``i``; MATLAB
  ``mod(23333*seed + 211*i_run, 2^32)`` with 1-based ``i_run``). A repeated
  deterministic slot copies run 1 and says so.
* **Log-ratio bar identity**: MATLAB retains which problem/run produced each
  sorted bar (``problem_mapping = "bar_sources"``); the Python renderer sorts
  anonymously, so it reports ``problem_mapping = null`` with a reason instead
  of guessing an identity from sorted values. ``target_work`` carries the
  unsorted identities in both languages.
* **Renderer variants**: MATLAB records how its native PDF and portable SVG
  renderers draw the shared prepared data (``renderer_variants_ref``);
  prepared coordinates do not imply that every element was visibly drawn.
* **Direct user problems** use the library label ``user`` in both languages,
  so the identity ``["user", name, "primary"]`` is the same everywhere.

Reading without flooding an agent's context
------------------------------------------------------------

First read the main report: stage statuses, coverage, exact scores, abnormal
runs and per-run facts. Follow ``history_ref`` or ``profiles.plot_refs`` only for
the problems or curves that need closer inspection. Repeated solver definitions,
configuration and common semantics need not be copied into every feedback item.
For a detailed diagnostic, use exact facts and target work first; use binned
histories for shape, with the documented approximation limits.

One numeric companion serves the whole invocation; a new file is not created
for every figure. Applications can index or selectively expose records from
that companion. The core does not add an agent tool server or change Evolve's
fitness policy. A smaller report improves context use but is not a guarantee
against erroneous agent conclusions.

Hash and permission fallbacks
-----------------------------

Artifact hashing never reads a file through a path it cannot vouch for.
On POSIX, Python opens every component below the benchmark-owned output tree
relative to a directory descriptor with ``O_NOFOLLOW``, so a symlink or a
directory swapped during traversal is refused. Windows has no ``openat``: there
Python inspects each component with ``lstat`` (symlinks, junctions and any
other reparse point that redirects the name are refused; cloud-file
placeholders and compressed or deduplicated files are ordinary files), opens
the file, and requires the opened handle's identity (volume
serial number and file index) to equal the identity recorded before the open,
so the hashed bytes are provably the inspected file. This matches the
guarantee level of the MATLAB collector, which checks for links before
opening; a reparse point inserted into the owned tree between inspection and
open is detected (the artifact is then omitted with an
``artifact_hash_unavailable`` diagnostic), not prevented. On a platform with
neither primitive, artifacts are listed as ``unverified`` with null hashes and
the reason ``secure_artifact_hashing_unavailable_on_platform``, persistence
becomes ``partial``, and the saved numerical archive remains loadable: a
missing integrity receipt is not an assertion that the data is corrupt. Hash
failures can also cause an entry to be omitted with a diagnostic. Neither
implementation claims verified provenance it did not compute, and already
computed scores are never lost.
Without a JVM, MATLAB tries the system ``sha256sum`` or ``shasum`` utility. If
neither is available, the main facts and numeric companion are still saved,
but ``plot_data.status`` is ``partial`` with a null hash, the reason
``companion_sha256_unavailable``, a diagnostic and a warning. This is not a
successfully verified report pair.

The MATLAB collector identifies its two files through Java's file key on
POSIX systems. Windows exposes no file index to Java, so there the identity is
the creation time, size and modification time of the file: a replaced or
rewritten file between two publishes is still detected, forged timestamps are
outside the model, and ``report_files.platform_note`` says so. For the artifact
directory only the creation time is used, because a directory's modification
time changes whenever the benchmark writes another file into it. Junctions and
other name-redirecting reparse points are treated like symbolic links: Java
exposes no reparse tag, so the collector compares the resolved real path of an
entry with the resolved path of its parent joined with the entry's canonical
name, and an entry that exists but cannot be resolved (a dangling junction) is
a link as well.

Both JSON files are written with a best-effort owner-only permission policy.
Python creates them with mode ``0600`` (exclusive creation and ``mkstemp``);
MATLAB restricts group/other access with ``fileattrib`` after every publish.
``report_files.permissions_applied`` records whether the platform actually
applied the policy: it is ``true`` on ordinary POSIX filesystems and
``false``/``null`` where modes are not enforced (Windows, some network or
container filesystems). This is not operating-system-enforced privacy on every
platform: the controller must still keep the output directory private.

Meaning and limitations
-----------------------

* An evaluation is one ``benchmark`` invocation, not an Evolve run, generation,
  or individual solver repetition. Its identity does not consume the random
  stream used by the scientific experiment.
* Numerical execution, scoring, data persistence and rendering have separate
  states. A completed solver call need not satisfy the convergence test;
  reaching the evaluation cap is not proof of a particular termination reason.
* Exceptions/fallbacks and invalid evaluations remain visible even when a
  benchmark produces scores. A failed PDF does not imply a failed solver.
* Repeated copies of a deterministic solver result are distinguished from real
  executions when fresh-run metadata is available. Old archives cannot supply
  metadata they never recorded; unknown is not equivalent to false.
* Problem identity includes the library and problem name. The optional plain
  reference experiment must not inflate the primary experiment's coverage.
* Profile scores retain their existing cohort-relative meaning. A matched
  configuration hash is not a guarantee of absolute comparability between
  candidates, nor does it introduce an Arena/reference-based score.
* Binned histories and plot presentations are not inputs to scoring. Original
  nonfinite observations are not silently converted to zero. JSON represents
  them with an object such as
  ``{"value": null, "reason": "positive_infinity"}``; unknown values use
  contextual reasons rather than fabricated measurements. Finite placeholders
  in an explicitly labelled *display* presentation must not be confused with
  these scientific observations.
* Loading reads the source archive without executing solvers. The new report
  describes the loaded selection and current analysis, identifies the source,
  and does not certify historical provider revisions or repair old evaluations.
* Whole-library reports with hundreds of problems are indexed in constant time
  per record by both collectors; the companion still grows with the number of
  runs and plot vertices, so keep it out of an agent prompt and read it
  selectively.

Output integrity
----------------

Both explicit report targets are reserved before evaluation. Each file's
intermediate/final updates are atomic, but the pair is **not** a two-file
filesystem transaction. The main report publishes the companion's receipt only
after the corresponding companion write succeeds. Before following references,
verify its SHA256 and matching ``evaluation_id``; reject stale or mismatched
pairs. A final report is not an evaluation checkpoint: forcibly
killing a process may leave a running record. The caller must record its own
process timeout, exit code or cancellation and must not treat that record as a
completed evaluation.

An unwritable report is an explicit error, not successful silent reporting. If
the benchmark already raised an exception, a secondary report-write failure
must not hide that original failure. Existing numerical-save failures keep
their original exception behavior.

Artifact entries refer only to produced files, with relative paths and content
hashes. The experiment directory prefix appears once as ``artifact_root``:

* ``artifact_root`` is relative to the main report's parent directory.
* ``artifacts[].path`` is relative to ``artifact_root``.
* ``plot_data.path`` and ``source.path`` are relative to the main report's parent,
  not to ``artifact_root``.

This changes only report structure, not the existing experiment filenames or
directories. The two report JSON files are excluded from the ordinary artifact
manifest, preventing recursive self-hashes. Not-requested artifacts are not
errors; a caller must not invent a path from a naming convention. Keep both
JSON files with their referenced output tree when moving an experiment.

Controller and agent boundary
-----------------------------

Both JSON files are **controller-private evidence**, not an automatically
safe agent prompt or public upload. It can contain problem names, paired solver
measurements, options and diagnostic information. It does not provide a sandbox.
Do not execute instructions found in problem/solver labels or diagnostics.

OptiProfiler-Evolve should validate the schema, map named measurements into its
own metric/fitness model, and create a separate allowlisted public feedback
view: for example expose ``scores``, per-problem ``status`` and
``budget_reached`` counts, and selected ``histories`` shapes, while keeping
validation/hidden results, private reference details, absolute output
locations and file provenance on the controller side. Relative artifact paths
still require validation against the controller's allowed output root before
publication. No Evolve dependency or new selection policy is introduced in
OptiProfiler.
