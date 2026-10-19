Reference targets
=================

Without a reference, profiles retain the V1 rule: use the least merit
observed in the selected solver histories, initial values and optional plain
run. This is an empirical comparison value, not a known optimum or lower bound.

A problem author or provider can instead supply raw reference facts:

.. code-block:: python

    reference = ProblemReference(fun=0.0, maxcv=0.0, kind='feasible',
                                 source='analytic optimum')
    problem = Problem(objective, x0, reference=reference)

MATLAB accepts the equivalent struct with fields ``fun``, ``maxcv``, ``kind``
and ``source``. The two kinds are ``feasible`` (zero violation) and ``merit``
(finite nonnegative violation). Optional ``x_ref`` is original-coordinate
source metadata and is not transformed by features or queried for metadata.
Built-in coordinate changes retain even nonzero reference violation because
scoring uses original-coordinate residuals. ``custom`` and ground-truth
``quantized`` invalidate the reference and use the default empirical target.
The first fallback reason is retained across compositions and save/load.

Standard profiles convert the effective featured reference as
``merit_fun(fun, maxcv, maxcv_init)`` separately for each run. This target
replaces the empirical minimum, including when it is beaten; measured merits
are never clamped. A custom merit need not preserve the ordering implied by
a raw objective lower bound. Its converted value is a scoring target, not an
automatically certified lower bound for that callback.
A converted target must be a finite real scalar. For example, the default
merit may discard a highly infeasible raw record by returning infinity; the
benchmark then reports an error instead of silently using a different target.

Saved experiments retain the raw records. Solver-free load uses the current
merit function, not a cached converted value, and problem/solver selection
keeps each reference attached to its original problem. Old archives without
reference records keep the V1 rule. To reproduce the same comparison, retain
the same scoring options, including ``run_plain``.

The optional EvalReport includes one effective raw ``reference_records`` entry
per stored run. A null entry means the empirical V1 target was used; an absent
field means the archive did not retain these facts. Reporting observes these
records without another merit-function call. Its provenance text is a bounded
metadata projection; the numerical archive retains the full raw records.

Feature propagation follows scoring truth, not the constraints as represented
to a solver. Built-in coordinate changes retain raw values even for nonzero
violation. No eager point transport or custom reference callback is required.
Providers choose which available reference to load through their own options;
the core adds no reference selection option to ``benchmark``.

Compatibility with earlier reference records
--------------------------------------------

New code uses only the two kinds above and does not supply ``mapping``.
For existing providers and saved records, ``lower_bound`` and ``optimum``
normalize to ``feasible``; ``best_known`` and ``target`` normalize to ``merit``.
The old ``feasible_objective/1`` and ``objective_violation/1`` mapping tokens
remain accepted in both languages but are removed from normalized records.
The former requires zero violation. Unknown tokens and callbacks are rejected.
The exact older record with fields ``merit``, ``kind``, ``source`` and
``mapping='feasible_objective/1'`` is read as ``fun=merit, maxcv=0``. Other
precomputed-merit layouts are not guessed or converted.
The Python constructor retains its optional fourth positional ``mapping``
argument solely for source compatibility; it is not an additional feature
mapping API. A ``custom`` stage still falls back to the default reference.
