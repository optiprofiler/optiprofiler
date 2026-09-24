# Reference contract: the feasible reference fact of a problem

This contract is implemented on the development line. It originated on
`codex/reference-facts-scalar`, based on
`54d42bd550d0493db661711715f868518d433924`, and was hardened on
`codex/reference-facts-scalar-followup`. It supersedes the
point-carrying record of the candidate branch `codex/reference-facts`
(`0aea026f52b9d96603054178a2044b88acce27b3`), which is left untouched.

## 1. The record

`Problem.reference` is optional author/provider metadata. When present it is a
record of exactly four fields.

| Field | Content |
| --- | --- |
| `merit` | A finite real scalar. |
| `kind` | `lower_bound`, `optimum`, `best_known` or `target`. |
| `source` | Non-empty provenance text, for example `author` or a citation. |
| `mapping` | A token of a closed registry. The only token is `feasible_objective/1`. |

There is no reference point, no objective/violation pair, no callable and no
user-defined mapping. An omitted reference means unknown; it is never a bound
at the initial point.

Python: `Problem(..., reference=None)` accepts a `ProblemReference` or a
mapping with exactly the four fields; `problem.reference` is the immutable
record or `None`. MATLAB: the optional input struct field `reference` is a
struct with exactly the four fields; `problem.reference` is that struct in the
canonical field order or `[]`, and it cannot be assigned.

### 1.1 Every kind is a claim over feasible points

- `lower_bound`: `f(x) >= merit` for every feasible `x`.
- `optimum`: `merit` is the exact optimal value of the objective over the
  feasible set.
- `best_known`: `merit` is the objective value of a known feasible point, hence
  an upper bound on the optimum.
- `target`: a level of the objective that the author chose as a target for
  feasible points. It carries no mathematical claim.

No kind is derived from solver output, and no kind is inferred when the field
is omitted.

### 1.2 The mapping token fixes how the scalar is read

`feasible_objective/1` means: `merit` is an objective value over feasible
points. It therefore equals the merit of a feasible point under every merit
function with the *feasible identity*

    merit_fun(f, 0, maxcv_init) == f    for every maxcv_init.

The default merit function has this identity for every initial violation,
including `NaN`. The tests pin it for `maxcv_init` in `{0, 0.2, 50, NaN}`.

The registry is closed. It is a literal list in each language with no
registration entry point. In Python the record class cannot be subclassed and
validation reads a module-private literal, so neither a subclass nor an
assignment to the informational attribute `ProblemReference.MAPPINGS` extends
it; in MATLAB the list is a literal inside the private validator. A mapping is
never a callable or a function handle.
An unknown token is rejected; it is matched exactly, so another version
(`feasible_objective/2`), another case or surrounding blanks are different
tokens and not spellings of the known one. A changed reading of the scalar
gets a new token, never a new meaning for an old one.

### 1.3 Why the record is not "a public scalar property"

`doc/DEVELOPMENT.md` asks for reference facts that are not a new public scalar
`Problem` property. The record honors this: a naked scalar is rejected, and the
scalar is never read without its kind, its provenance and its mapping. All four
fields are required on input. In particular an omitted `mapping` is not read as
`feasible_objective/1`: the author states how the scalar is to be read.

## 2. What the record is not

**Not a run-history minimum and not the dynamic cohort minimum.** The profile
baseline (`merit_mins`) is the least merit observed over the selected solvers,
runs and evaluations, including the initial point and, with `run_plain`, the
plain run. It changes when solver membership changes, is recomputed on every
load, and is never stored in a problem. The reference is an absolute claim
about the problem that does not depend on which solvers ran. The two are never
merged, and nothing writes a cohort minimum into `Problem.reference`.

**Not a floor for run merits; never a clamp.** On a constrained problem the
merit of a run may be *below* the reference. A merit function tolerates or
penalizes small violations, so an infeasible point can have a lower merit than
every feasible point. With the default merit function, a point whose violation
is at most `v1 = min(0.01, 1e-10 * max(1, maxcv_init))` has merit `f`, and a
point with violation in `(v1, v2]` has merit `f + 1e5 * (v - v1)`, which is
below the feasible optimum whenever the objective falls faster than the penalty
grows. Such merits are legitimate. They are recorded, archived and scored as
they are, and no code path clamps a history, a merit or a score to the
reference. This holds for every kind: not even a rigorous `lower_bound` bounds
the merits of infeasible points.

**Not comparable with run merits under an arbitrary merit function.** The
scalar is stated for feasible points. Before a consumer compares run merits
with the record under a custom `merit_fun`, that function must be known to
preserve the feasible identity of section 1.2. A merit function without it (for
example `f + 1`, or `f * (1 + maxcv_init)`) puts run merits in another space
than the record. No consumer exists on this branch; the requirement is
documented for the first one, and the tests contain the probe it has to make.

## 3. Validation and deserialization

Validation is structural and evaluates nothing, so a reference adds no callback
call to building or loading a problem. (The Python constructor probes the
constraint callbacks at `x0` for their dimensions, with or without a reference;
the reference is validated before that, so a malformed record is rejected
before any callback of the problem is touched.)

Rejected at construction, in both languages: a non-finite, non-real or
non-scalar `merit` (logical values, text and one-element containers included);
an integer `merit` that a double cannot hold exactly; an unknown or non-text
`kind`; an empty, blank or non-text `source`; an unknown, non-text or callable
`mapping`; a missing field; an unknown field; a naked scalar; and any record
with a field of the superseded layout (`fun`, `maxcv`, `point`), even when the
four current fields are present as well. Nothing is repaired, defaulted or
guessed.

A serialized record that the running version does not accept is read as
unknown or raises. It is never reinterpreted: no field is read as the merit of
another layout, and no mapping is guessed.

| Case | Outcome |
| --- | --- |
| Object serialized before the record existed | unknown, no warning |
| Python pickle of the current layout whose record is rejected (unknown mapping of a later registry, non-finite merit, foreign field) | unknown (`None`) with a `RuntimeWarning`; the problem stays usable |
| Python pickle that restores instance state directly (the superseded layout) | `TypeError`; the stream is refused as a whole |
| Python instance whose stored value is not a validated record | unknown |
| MATLAB file whose stored record is rejected (the superseded five-field layout, an unknown mapping, a non-finite merit) | the object loads; `reference` reads as `[]` with the warning `MATLAB:Problem:reference_StoredRecordRejected` |

MATLAB detail. `reference` is a dependent, validated view of the stored record.
Its set method is protected and stores the value as is, so loading a file never
raises in that method; without it MATLAB cannot assign the value saved by the
superseded layout, hands a struct to `loadobj` and warns that the constructor
must preserve the class. Every read goes through the validator, so a caller
sees a valid four-field record or `[]`, whatever a file left in the object.

Python featured-problem pickles. A built `FeaturedProblem` is restored from
its saved instance state, not through its constructor. Single stages and
compositions therefore keep their sampled affine map without calling feature
callbacks again. This is a compatibility path for trusted pickle input, not a
safe interchange format. New construction safeguards do not recompute stored
initial points or histories in existing pickles or experiment archives; loading
does not migrate or repair previously recorded results.

The reducer stores state in its third item, after the instance is registered,
so callback-owner cycles refer to that same restored instance. Instance
dictionaries and initialized subclass slots (including inherited and private
slots) are preserved; uninitialized slots remain uninitialized. Each stage
runtime restores the read-only flags of its own cached affine arrays. No
general traversal of callback objects, closures or other user-owned state is
performed.

The first reducer's state-as-argument layout remains readable: its state
dictionary is kept by identity because pickle may still be filling it when a
callback-owner cycle calls the reconstructor. Compositions written before
there was a reducer are also read without resampling. The oldest single-stage
pickles, whose reconstruction lacks the required constructor arguments and a
kept transformation, still raise `TypeError`; reading them must not silently
call a stateful modifier again.

## 4. Propagation (`FeaturedProblem.reference`)

The record is retained unchanged or reported as unknown. Nothing is
transported, because the record holds no point, and nothing is derived. The
rule is the same for every kind: each is a claim over the feasible points, so
each survives exactly the stages that leave the feasible points and the
objective values over them unchanged up to a change of variables.

| Stage | Reference | Reason |
| --- | --- | --- |
| `plain` | retained | identity |
| `noisy`, `truncated`, `random_nan`, `nonquantifiable_constraints`, `unrelaxable_constraints` | retained | observation only: the truth recorded in the histories is the original problem's |
| `quantized`, `ground_truth=false` | retained | observation only: the mesh changes what the solver sees, not the truth |
| `perturbed_x0` | retained | only the initial point moves |
| `permuted` | retained | coordinate permutation |
| `linearly_transformed` | retained | invertible affine change of variables built by the framework |
| `custom`, option names within `{mod_x0, mod_affine}` | retained | `mod_x0` moves the initial point; `mod_affine` must pass the numerical validation and transport checks of section 7 before the featured problem is built |
| `custom`, any other option (`mod_fun`, `mod_cub`, `mod_ceq`, `mod_bounds`, `mod_linear_ub`, `mod_linear_eq`, or an option added later) | unknown | arbitrary user code that may change values, constraints or bounds; the framework cannot prove that the feasible reference survives |
| `quantized`, `ground_truth=true` (the default) | unknown | the truth becomes the mesh problem; no proof that the same feasible reference holds is available |
| any other stage name | unknown | fail closed |

A composition is safe only if every stage is safe. Each stage view inherits its
predecessor's record only if its own stage retains it, so once a stage reports
unknown no later stage can restore the record.

The custom rule is a whitelist on the option *names*. It does not look at what
a callback does: a `mod_fun` that returns the predecessor's value unchanged
still makes the reference unknown. Only names are read and no callback is
called, so a reference costs no callback call; the tests compare the call
counts of every problem and custom callback with and without a reference.

`quantized` with `ground_truth=true` stays unknown for every kind. A proof that
the mesh problem has the same feasible reference would need more than the
record holds, so no such proof is attempted and the rule has no kind-specific
exception.

## 5. Removed relative to the superseded branch

`point`, `with_point`, `fun` and `maxcv` no longer exist on the record, and the
affine inverse that the MATLAB stage views kept for transporting a point
(`affine_inverse`) and the transport function (`transportProblemReference`) are
gone. The `reference` property of `Problem` and `FeaturedProblem` and the
save/load plumbing are kept in both languages. The internal MATLAB scoring
method of the stage views is named `referenceValue`, because `reference` is the
`Problem` property.

## 6. Out of scope

- No consumer: no benchmark option, no score, no profile formula, archive or
  report reads the reference. The public `benchmark` signature is unchanged.
- No offline catalog, scoring manifest or independent scorer.
- No change to solver interfaces, seeds, random streams, history or archive
  layouts, profile formulas, the report schema, providers, locks, gitlinks,
  versions or paper references. Section 7 describes the correction or refusal
  of previously accepted affine initial points, not a new experiment format.

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
