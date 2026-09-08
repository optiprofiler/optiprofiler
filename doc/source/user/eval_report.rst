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
``None`` or an empty string), the previous behavior is retained. An explicit
path requests two JSON files: the small main report and one adjacent numeric
companion. For ``evaluation/eval_report.json``, the companion is
``evaluation/eval_report.plot_data.json``. Both targets must be new: a collision
with either file is refused before solvers execute. This also prevents a
load/replot operation from replacing its source report.

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

MATLAB
------

The equivalent options-struct field preserves MATLAB's three outputs::

    options = struct('plibs', {{'s2mpj'}}, 'ptype', 'u', ...
        'mindim', 1, 'maxdim', 10, 'score_only', true);
    options.report_path = 'evaluation/eval_report.json';
    [scores, profile_scores, curves] = benchmark({@solver_a, @solver_b}, options);

The MATLAB implementation emits the same schema and one-based solver/run
indices as Python. Runtime provenance and random-seed policies remain
language-specific; the report does not change either policy to make them match.

Report structure
----------------

The :download:`main-report schema <../_static/eval_report.schema.json>` and
:download:`numeric-companion schema <../_static/plot_data.schema.json>` describe
the common Python/MATLAB contract. They require no additional runtime dependency
in OptiProfiler. Consumers should reject bare JSON ``NaN``/``Infinity`` tokens,
validate the schemas, and check tensor dimensions, references and artifact paths.
This remains an unreleased v1 report format, independent of package versions.

.. list-table:: Main fields
   :header-rows: 1
   :widths: 25 75

   * - Field
     - Meaning
   * - ``configuration`` / ``source``
     - Effective current options and, for load, the original archive receipt.
   * - ``stages`` / ``coverage``
     - Independent execution, scoring, persistence and rendering states;
       selected, loaded and completed primary problems.
   * - ``problems[].runs[]``
     - Evaluation counts, budgets, invalid values, returned/initial/best
       objective, constraint and merit values, exceptions and fallback flags.
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

Solver/run/sample indices in observation fields use one-based indexing,
including Python JSON. Options recorded in ``configuration`` retain the
original language API conventions; for example, Python ``solvers_to_load``
still uses zero-based indices.
The profile-score axes are solver, tolerance, history/output, and profile type.
Objective, constraint and merit ``best`` values are componentwise minima: they
need not belong to the same evaluated point. Work summaries reuse the existing
profile convergence calculation; they do not introduce a new stopping rule.
The main report keeps these meanings once in ``semantics`` rather than repeating
the same prose in every run. Missing observations still have contextual reasons;
absence or ``null`` must not be interpreted as a successful zero measurement.

Three different kinds of data
-----------------------------

The main report is a concise index of facts. Its numeric companion contains
``histories``, ``plots`` and ``target_work``. The original H5/MAT archive is
unchanged and remains the source for complete reanalysis.

* **History summaries** describe each problem/solver/run and each metric
  separately, excluding padded values beyond actual evaluations. Short histories
  of up to 64 observations use ``representation="exact_samples"`` and retain
  all values. Longer histories use ``representation="bin_extrema"``: at most 32
  contiguous bins retain their first/last observations, finite minimum/maximum
  and true evaluation positions, plus NaN and positive/negative infinity
  counts. Ties use the first occurrence. The binned representation retains extrema
  within each bin, not the complete order of observations inside it. Unavailable
  histories, zero evaluations and a retained history shorter than the reported
  evaluation count remain distinct.
* **Numeric plot presentations** contain the actual prepared coordinates,
  central curves and error bands from the same numeric preparation used by the
  renderer. They are not another uniform 64-point sample. ``exact_rendered_data``
  means exact *prepared display data*, not a lossless copy of the raw experiment,
  pixel-identical rendering across machines, or evidence of a successful PDF.
  History display clipping, nonfinite placeholders, shifts, raw/cumulative
  views, block aggregation and coordinate transformations are recorded.
  ``observation_scope="rendering_inputs"`` identifies observations captured from
  the actual rendering path. When drawing is not requested, a representation
  can instead use ``retained_scoring_observations``. These scopes are not
  interchangeable: a stateful custom merit callback can produce different
  values when the existing renderer calls it again. Reporting does not invoke
  that callback an extra time to make the two agree. MATLAB renderer-variant
  metadata also distinguishes native bands/bars from the portable SVG view;
  prepared coordinates do not imply that every element was visibly drawn.
* **Target work** retains the existing profile calculation's work arrays and
  their problem/solver/run/tolerance identities. History-based work is the first
  evaluation meeting the existing target. Output-based work is the evaluation
  count associated with a passing returned output, not the first history hit.
  Never infer either value from a compressed curve.

The scalar ``best`` ignores NaN and can be infinite. A bin's ``finite_min`` and
``finite_max`` intentionally exclude infinities; they answer a different
question. Exact scalar facts, abnormal/fallback flags and invalid counts remain
available even if a short visual summary hides an event.

The existing history display order is across-run aggregation, shift, cumulative
view (if requested), and block aggregation. In particular, cumulative mean
curves are not silently changed into means of independently cumulative runs.
Python retains its population standard-deviation bands; MATLAB retains its
sample normalization. The numeric presentation records the relevant convention
instead of changing established numerical behavior to force cross-language
equality.

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

Secure artifact hashing is capability-dependent. When a platform cannot
support safe hashing, Python can mark an artifact ``unverified`` with a reason
and null hash. Hash failures can also cause an entry to be omitted with a
diagnostic (including MATLAB's current fallback). Both implementations report
partial persistence instead of claiming verified provenance or losing already
computed scores.
Without a JVM, MATLAB tries the system ``sha256sum`` or ``shasum`` utility. If
neither is available, the main facts and numeric companion can still be saved,
but their receipt is marked partial with a null hash and an explicit diagnostic.
This is not a successfully verified report pair.

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
view. Validation/hidden results, private reference details, and file provenance
must remain on the controller side. Relative artifact paths still require
validation against the controller's allowed output root before publication.
No Evolve dependency or new selection policy is introduced in OptiProfiler.
