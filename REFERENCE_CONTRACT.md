# Reference contract: the feasible reference fact of a problem

Branch `codex/reference-facts-scalar`, based on
`54d42bd550d0493db661711715f868518d433924`. This contract supersedes the
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
| `custom`, option names within `{mod_x0, mod_affine}` | retained | `mod_x0` moves the initial point; `mod_affine` is validated by `A * inv == I` when the problem is built, so it is a valid affine coordinate change |
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
- No change to solver calls, seeds, random streams, histories, profile
  formulas, archives, the report schema, providers, locks, gitlinks, versions
  or paper references.

## 7. The posed problem is the scored problem (affine safeguard)

This section was "out of scope" when the contract was written and is fixed by
the affine safeguard commit that follows it on this branch.

The defect. For `linearly_transformed` and for `custom` with `mod_affine`, the
bounds modifier decided with an exact test on the *inverse* whether the map is
diagonal (MATLAB `isdiag(inv)`), while the linear modifiers decided with an
exact test on the *matrix* (`isdiag(A)`). With `A = diag(2, 4)` and a supplied
inverse carrying off-diagonal entries of size `1e-17`, the bounds became
infinite and no bound row was added: the bounds left the problem handed to the
solver, while the truth went on scoring them (it reported a violation of 19 at
a point mapped far outside the box). A retained reference then described a
problem that the solver had not been given.

The rule now, the same in both languages. The pair is validated before any
simplification: real, finite arrays of matching sizes;
`norm(abs(inv) * abs(A), inf) < 1 / eps`, a condition number that scaling the
rows of `A` cannot change, so that an exact but badly scaled transformation is
not refused; and, for a supplied inverse,
`norm(A * inv - I, 'fro') <= 1e-8 * n`. One decision is then made from both
matrices and read by the bounds and by both kinds of linear constraints: an
off-diagonal entry is negligible if
`abs(M(i, j)) <= n * eps * min(abs(M(i, i)), abs(M(j, j)))`, for `M = A` and for
`M = inv`. If both are diagonal in this sense the bounds stay bounds, scaled by
`diag(inv)`; otherwise every finite bound is posed as a linear row of `A`,
which needs `A` only and is valid for every invertible `A`. The tolerance
therefore chooses a representation and can never decide whether a bound is
posed. A finite bound that would overflow when scaled raises instead of
becoming infinite. So does a supplied `mod_linear_ub` or `mod_linear_eq` under a
`mod_affine` that is not diagonal, without `mod_bounds`, if the problem has
bounds that the framework would pose as those rows: a supplied linear modifier
replaces the rows verbatim, so the bounds would have nowhere left to go.

Follow-up (the hardening commit after the safeguard). An independent audit of
the rule above found four ways in which the posed problem could still differ
from the scored one, and the rule is now as follows; where this paragraph and
the previous one differ, this one holds.

- The decision is exact. `abs(M(i, j)) <= n * eps * min(...)` called an entry of
  `4e-16` next to a diagonal of 1 negligible. With bounds of `1e16` that entry
  moves the feasible set by 4: `y = (-4, 1e16)` is mapped to the vertex
  `(0, 1e16)` of the original box, and the box posed in the new variables
  rejected it. What an off-diagonal entry of `A` moves depends on the size of
  the other variable, which no tolerance on the entry knows. The bounds stay
  bounds only if `A` is exactly diagonal (then the set is a box) and
  `abs(inv(i, i) * A(i, i) - 1) <= 8 * eps` (then `diag(inv)` times a bound is
  that bound to roundoff; a supplied inverse is otherwise only held to `1e-8`).
  Off-diagonal entries of `inv` change no feasible set and are not read. The
  pair that `linearly_transformed` builds without rotation is exactly diagonal
  and reciprocal to `1.5 * eps`, so its bounds are bitwise what they were.
- One transformation per problem and seed. The modifiers and every evaluation
  of the single-feature path each called `mod_affine` again. User code with a
  state (a counter, the global random stream instead of the stream it is
  handed) then gave the bounds one map and the linear rows another: with a
  diagonal and a dense answer in turn, the bounds became infinite and no bound
  row was added. The runtime now produces and validates the triple once, keeps
  it (read-only in Python; saved with the kernel in MATLAB, so that a loaded
  featured problem goes on with the map its structure was built with), and
  every reader gets that one. A specification keeps nothing.
- `inv` has to invert `A` from both sides, and neither set of units matters.
  With `A = diag(1e-8, 1e8)` and `inv = [1e8, 1e-8; 0, 1e-8]`, `A * inv` is the
  identity to `1e-16` while `inv * A` misses it by 1. Requiring
  `norm(inv * A - I, 'fro') <= 1e-8 * n` as well would refuse a rotation with
  one new variable in units of `1e13`, which is consistent to roundoff. The
  plain rule on `A * inv` already refused the same rotation with one original
  variable in those units; and for a rotation whose products happen to cancel
  exactly in plain arithmetic (entries `t` and `2 * t`) NumPy, whose kernels
  are fused, refused it while MATLAB accepted it. Each entry is therefore
  measured against the terms
  it is summed from: `norm(E, 'fro') <= 1e-8 * n` for
  `E = (A * inv - I) ./ max(1, abs(A) * abs(inv))` and for
  `E = (inv * A - I) ./ max(1, abs(inv) * abs(A))`. Where the terms are below 1
  this is the plain rule; it is never stricter than the plain rule on
  `A * inv`, and an error of the size of the terms is refused at every scale.
- Nothing finite leaves the floating-point range unnoticed. `xu - b` with
  `1e308 + 1e308` is infinite before anything is scaled, so the scale guard saw
  "no bound"; as a row it had an infinite right-hand side, which counts as no
  constraint; `bub - aub * b` and `aub * A` gave `-Inf`, `NaN` and `Inf`. Every
  shifted bound, right-hand side, composed row and the pulled-back initial
  point is now checked, for `linearly_transformed` as well, and raises.

Two conversions were found on the way. MATLAB kept the data of a `Problem`
(`x0`, `xl`, `xu`, `aub`, `bub`, `aeq`, `beq`) in the class they were given in,
and combines an integer with a scalar double in integer arithmetic: with `int32`
data the truth called `x = 4.3` feasible for `xu = 4`, and a change of variables
posed the bounds `-1` and `2` for `-0.25` and `1.75`. The data are now double
precision numbers, as in Python. In a Python composition, an integer beyond the
range of a float returned by a custom callback raised `OverflowError` instead of
the `ValueError` naming the stage.

Consequence for this contract. `linearly_transformed` and `custom` within
`{mod_x0, mod_affine}` retain the reference (section 4). That is a claim about
the problem handed to the solver, and it holds because that problem is now
always the scored one in new coordinates: a transform that cannot be
represented raises, so no featured problem, and no reference, exists for it.
