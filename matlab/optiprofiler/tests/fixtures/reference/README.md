# Reference-fact test fixtures

## `legacy-point-record.mat`: objects of the superseded record layout

The file holds two live objects saved by the MATLAB implementation of the
superseded candidate branch `codex/reference-facts` at source
`0aea026f52b9d96603054178a2044b88acce27b3`, MATLAB R2026a on syu-ubuntu,
2026-09-19, format `-v7`:

- `problem`: `Problem(struct('name', 'LEGACY', 'fun', @(x) sum((x(:) - [1;
  2]).^2), 'x0', [0; 0], 'reference', record))`;
- `featured`: `FeaturedProblem(problem, Feature('plain'), 5, 0)`;

where `record = struct('fun', 0, 'kind', 'optimum', 'source',
'codex/reference-facts@0aea026', 'point', [1; 2])`. In that layout `reference`
was a stored property holding the five fields `fun`, `maxcv`, `kind`, `source`
and `point`.

The current contract stores exactly `merit`, `kind`, `source` and `mapping`. The
fixture pins what happens to a record of the superseded layout: both objects
load as objects of their classes, their objectives still evaluate, and their
`reference` reads as `[]` (unknown) with the warning
`MATLAB:Problem:reference_StoredRecordRejected`. The legacy `fun` is never read
as a `merit`.

SHA-256:

- `legacy-point-record.mat`
  `ef4d3da28f1b1a655113b56d4f5ea7cf6b25900df24c0b0dd10619d41e778d6a`

## `RawReferenceProblem.m`: a stored record that skips validation

A `Problem` subclass whose constructor stores a raw value as the reference
record through the protected setter. It emulates, without a binary file, what
a file written under another mapping registry or a damaged file leaves in an
object (for example the token `feasible_objective/2` or a non-finite merit),
so that the tests can pin that such a record reads as unknown.

## `reference_invariance_select.m`, `reference_invariance_load.m`: the reference-invariance library

A two-problem provider whose loader states reference facts only when the
environment variable `OPTIPROFILER_TEST_STATE_REFERENCES` is `1`. The library
name, root, problem names and callbacks are the same in both variants, so a
benchmark with references can be compared bit for bit with one without.
