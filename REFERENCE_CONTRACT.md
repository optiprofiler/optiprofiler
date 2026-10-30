# Reference contract: convergence targets and source points

This is the current V2 development contract. Earlier scalar-only and
point-transport candidates are superseded; section 7 retains the existing
numerical safeguards for the feature transformations themselves.

## 1. The record

`Problem.reference` is optional. A stored (serialized) record has four
required fields; the Python `ProblemReference` constructor may omit `maxcv`,
which then defaults to zero and is stored explicitly:

| Field | Meaning |
| --- | --- |
| `fun` | Finite raw objective reference, not a precomputed merit. |
| `maxcv` | Finite nonnegative total violation paired with `fun`. |
| `kind` | `feasible` requires `maxcv=0`; `merit` allows nonzero violation. |
| `source` | Non-empty provenance text. |
| `x_ref` (optional) | Finite reference point in the original problem's coordinates. |

`merit` users are responsible for the selected merit function and reference.
Neither kind asks the engine to prove optimality. A feasible reference may
also be a value-only bound; absence of a point does not invalidate its value.
No user mapping token, callback, or point-transport option is required.

## 2. Convergence targets

An effective explicit record replaces the empirical target with
`merit_fun(fun, maxcv, maxcv_init)` separately for each problem/run.
Missing or invalidated references retain the V1 minimum over selected
solver histories, initial values and the optional plain-reference result.
`default` is this selection rule, not a third kind and not a value written
back into `Problem.reference`.

Initial violation may vary by run, so the converted merit is not cached as
a problem fact. An invalid converted target raises; it is not silently
replaced. The record itself is valid; it is incompatible with the selected
merit at that run's `maxcv_init`. Because both are known once the featured
problem of a primary run is built, a benchmark reports this before that run's
first solver call (plain-baseline records are never targets); a load with a
different merit reports it when scoring. Measurements below a supplied target
are retained, never clamped. A feasible-objective reference need not bound
merits at infeasible points.

A load uses one merit for every cached quantity and for the targets. With
`merit_fun`, all cached merits are recomputed. Without it, an archive with
explicit references is rescored entirely with the default merit (a warning
and the report diagnostic `cached_merits_recomputed` state when the cached
values change), and an archive without references keeps its cached merits
(V1 behaviour). An archive saved after a merit failure has no cached merits
and is scored with the merit of the load.

## 3. Features and original-coordinate points

| Stage | Explicit reference |
| --- | --- |
| `plain`, `noisy`, `truncated`, `random_nan`, `nonquantifiable_constraints`, `unrelaxable_constraints`, `perturbed_x0` | retained unchanged |
| `permuted`, `linearly_transformed` | retained unchanged, including nonzero `maxcv` |
| `quantized`, `ground_truth=false` | retained unchanged |
| `quantized`, `ground_truth=true` | default fallback, even with `x_ref` |
| `custom` | default fallback for every option set; no reference callback in this release |
| unknown stage | default fallback |

Scoring truth evaluates residuals in predecessor/original coordinates, not
the transformed constraints passed to a solver. Coordinate changes therefore
do not require transporting reference `fun/maxcv`.

`x_ref` is source metadata, not a point in featured coordinates. It stays
unchanged even after a permutation or affine stage. Maintaining the record
performs no point matrix operation, oracle evaluation, or random-stream draw.
Do not pass this original-coordinate point directly to a transformed solver.

A composition keeps a reference only if every stage keeps it. Once lost,
later stages cannot restore it. The first invalidating stage is retained in
`reference_fallback_reasons` as `custom_feature`, `quantized_ground_truth`,
or `unsupported_feature`. A null/empty reason means no feature-triggered
fallback, not a claim that a reference was supplied.

## 4. Validation

Validation evaluates no objective or constraint. Raw scalars must be real,
finite, non-logical, exactly representable as doubles, with nonnegative
`maxcv`; feasible records require zero violation. Kind and provenance must
be valid text. A supplied point must be a finite vector matching the original
problem dimension. Python owns a copy and exposes the immutable record;
MATLAB exposes a validated struct. A naked scalar and the ambiguous old
`point` field are rejected.

## 5. Compatibility and persistence

Earlier `lower_bound`/`optimum` inputs normalize to `feasible`;
`best_known`/`target` inputs normalize to `merit`. This also supports existing
provider packages. Earlier `objective_violation/1` and zero-violation
`feasible_objective/1` tokens are accepted as compatibility inputs, not
emitted in new records. Unknown tokens are rejected. Only the exact earlier
`merit/kind/source/mapping` feasible-objective draft is read as
`fun=merit,maxcv=0`; arbitrary cached merits are not reinterpreted.

Archives retain full raw records, original-coordinate points and fallback
reasons. Reports embed points of up to 256 coordinates; larger points have
an explicit dimension/omission record, with the full point in the numerical
archive when that archive exists. Reporting performs no extra merit calls.
Load recomputes targets with the current merit function without executing
solvers or reference-point oracles. Source archives remain unchanged.

Old archives without records/reasons remain readable. Canonicalizing their
record vocabulary does not rerun or repair their recorded experiments.
Rejected serialized records read as unknown with a warning or raise,
rather than guessing their meaning. Native featured-problem pickles preserve
the saved runtime maps and state without calling modifiers again.

## 6. Scope

No new benchmark option or solver signature is introduced. Providers choose
their references using existing library-owned loading options. No custom
reference-conversion hook, automatic reference optimization, or new Arena
scorer is introduced here. The V1 report schema remains immutable.

## 7. Affine validation and transport

An affine stage uses the change of variables `x = A * y + b`. Its structure,
truth evaluations and supported derivatives must use the same validated
triple `(A, b, inv)`, where `inv` is a candidate inverse. A failed check raises
during construction rather than silently dropping a constraint or retaining a
reference for a failed construction. These are floating-point safeguards, not
exact rank certificates, proofs of the author's reference or guarantees of exact
reversibility at every representable point. Here `eps` is double-precision
machine epsilon, approximately `2.22e-16`.

### 7.1 Validate the matrix independently of its claimed inverse

The triple is converted to double precision, then checked for real finite
entries and the required dimensions. Callback outputs, including integers
beyond `2^53`, use the rounded double values consistently; unlike a recorded
reference merit, they need not be exactly representable as doubles.

The established pair scale policy remains
`norm(abs(inv) * abs(A), inf) < 1 / eps`. For a user-supplied inverse, both
products must also satisfy

    norm(max(abs(P - I) - 64 * n * eps * T, 0), 'fro') <= 1e-8 * n,

where `(P, T)` is `(A * inv, abs(A) * abs(inv))` or
`(inv * A, abs(inv) * abs(A))`. Non-finite test quantities fail. The allowance
accounts for rounding of large cancelling terms; it does not independently
establish that a claimed inverse is valid.

User-supplied matrices therefore have an additional check that uses `A`
alone. Its rows and then columns are equilibrated by powers of two, using
their largest absolute entries. Reversing each scaling must recover every
entry exactly; a zero row or column, or any loss during scaling, refuses the
transform. Python requires `cond(E, 1) < 1 / eps` for the equilibrated matrix;
MATLAB requires `rcond(E) > eps`. These are independent numerical conditioning
estimates, not mathematical certificates, and may differ near the cutoff.
There is no extra factor `n` in this cutoff. The empty `0 x 0` matrix is
accepted. Internally constructed rotations and nonzero diagonal scalings
retain their existing checks without an additional cubic-cost factorization.

### 7.2 Preserve the box and explicit linear constraints

The diagonal shortcut is selected only when `A` is exactly diagonal and
`abs(inv(i, i) * A(i, i) - 1) <= 8 * eps` in every component. Bounds are then
scaled by `diag(inv)`, with endpoints exchanged for negative scales.
Off-diagonal entries of a supplied inverse do not select the representation.
Every other accepted map uses the generic representation: finite original
bounds become linear rows of `A`, and fixed variables become equalities.
Small nonzero entries of `A` are never discarded as a structural shortcut.

`mod_bounds` replaces the logical original box, independently of whether an
affine stage would represent it as bounds or generated rows. The replaced box
is not added again; original explicit linear constraints are still
transported. Without `mod_bounds`, a custom linear modifier is rejected when
its verbatim replacement would discard generated bound rows of that type.
Any `mod_bounds` makes the reference unknown, as in section 4.

### 7.3 Reject unrepresentable transport, allow ordinary rounding

Finite shifted bounds, transported linear rows and right-hand sides must stay
finite. A nonzero diagonally scaled bound must remain at least the smallest
normal magnitude. A transported row or right-hand-side entry is refused when
it has a product of nonzero factors but the sum of its absolute terms is below
that magnitude; cancellation of normal terms is not classified as underflow.
Both overflow and underflow raise rather than silently removing a constraint.

Bounds are checked as a pair after translation and after diagonal scaling.
If a strict original interval acquires equal floating-point endpoints,
construction raises. The translation check also applies when the box becomes
linear rows. A genuinely fixed original variable remains allowed.

For finite original `x0`, retaining any affine default candidate
`y = inv * (x0 - b)` requires finite `y`, `A * y + b` and rounding allowance,
and the established componentwise rounding test:

    alpha = 64 * n * eps
    abs(A * y + b - x0) <= alpha * (abs(A) * abs(y) + abs(b) + abs(x0)).

Only a custom user-supplied inverse has the additional test below before its
candidate can be retained unchanged:

    abs(A * y + b - x0) <= sqrt(eps) * min(abs(x0), abs(x0 - b)).

This extra threshold is about
`1.49e-8 * min(abs(x0(i)), abs(x0(i) - b(i)))` per component. It
determines whether to keep the inverse candidate unchanged or solve
`A * y = x0 - b` independently. It is independent of the claimed inverse and
of other coordinates: an inflated inverse or an unrelated large coordinate
cannot relax it. This prevents retaining an inaccurate inverse candidate
merely because its own large entries inflated the rounding allowance.
The smaller of the original and centered coordinate magnitudes prevents a
large coordinate origin from hiding error in a small `x0 - b`, or a large
shift from excusing error in a small original coordinate. There is no absolute
floor: any nonzero recovery error at a zero coordinate of either `x0` or
`x0 - b` triggers an independent solve. An omitted shift is zero.

The independently solved point is checked for finiteness and the established
rounding allowance only, not against the extra local threshold. Failure
raises an initial-point error. The threshold is a solve trigger, not a bound
promised for the final forward error or a proof of exact feasibility. Even an
honestly computed custom inverse can trigger a solve, and a well-conditioned rotation
that mixes small and large coordinates may not recover every small coordinate
within that threshold. Nor is small forward error guaranteed for arbitrary
ill-conditioned maps. This policy removes reliance on the supplied inverse
for a suspect start without imposing a new all-map forward-accuracy contract.
Intrinsic conditioning or mixed coordinate scales can cause forward drift even
after solving. Initial truth evaluation uses the posed start; the existing
MATLAB single-stage fallback still retries the original initial point when
an evaluation returns an empty result. Thus `fun_init` and `maxcv_init` do
not certify recovery of an originally feasible point.
Internally generated `linearly_transformed` maps retain their existing
finite/rounding-only policy; this custom-inverse solve trigger does not change
their starts. Forward drift for extreme conditioning or mixed coordinate
scales remains a limitation of that existing policy.

The allowance's nonnegative terms are scaled before addition to avoid
overflow in a valid large identity case. An overflowing absolute matrix
product is still a failed check; `Inf <= Inf` is never accepted as evidence.
Custom inverse candidates passing both tests are kept bitwise; no such promise
is made for all honestly computed custom inverses. Explicit `mod_x0` still
overrides the affine default initial point and is not required to map back to
the predecessor's `x0`. The existing policy for an originally non-finite user
`x0` is unchanged: it is transported
as given. These tests permit ordinary rounding, not exact equality.

Transport checks are local to each stage. The strict-interval test protects
the original box of that stage, not arbitrary intervals encoded as linear
rows or a proof of lossless transport across a composition.

MATLAB stores ordinary problem data (`x0`, bounds and linear constraints) as
doubles, as Python does, so integer arithmetic cannot change transported
bounds or truth evaluations.

### 7.4 Keep one runtime map and the established derivative semantics

Each stage runtime caches one triple, keyed by problem identity and seed.
Every structure, truth and derivative reader uses that triple; Python arrays
are read-only and MATLAB saves the triple with the kernel. Specifications
carry no sampled map. A problem must not be mutated in place while a featured
problem uses it. Deprecated `Feature.modifier_*` conveniences create a new
runtime per call and do not promise a shared map across separate calls.
Python restoration is runtime-local as described in section 3.

For a single `permuted`, `linearly_transformed` or custom affine stage, the
existing derivative wrappers apply the chain rule at `x = A * y + b`:
`A' * grad(x)`, `A' * hess(x) * A`, `J(x) * A` and
`A' * H_i(x) * A`. They use original callbacks, not noisy or otherwise
observed values; they consume no evaluation budget and record no history.
Absent derivatives remain absent. Other single stages retain their existing
passthrough behavior, and compositions still do not provide derivatives.
The numerical safeguards do not change these derivative or cache contracts.
