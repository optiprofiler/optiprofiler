# Reference-fact test fixtures

## `legacy-point-record.mat`: objects of the superseded record layout

The file holds two live objects saved by the MATLAB implementation of the
superseded candidate branch `codex/reference-facts` at source
`0aea026f52b9d96603054178a2044b88acce27b3`, MATLAB R2026a on syu-ubuntu,
format `-v7`, regenerated on 2026-09-20:

- `problem`: `Problem(struct('name', 'LEGACY', 'fun', objective, 'x0', [0; 0],
  'reference', record))` with `objective = str2func('@(x) sum((x(:) - [1;
  2]).^2)')`;
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

The objective must come from `str2func`. The first version of this file
(2026-09-19) wrote the anonymous function inside a generator function file, so
MATLAB stored that file's path, in a scratch directory under `/tmp`, in the MAT
file. Loading it warned
`MATLAB:dispatcher:UnresolvedFunctionHandle` wherever that scratch file did
not exist, which is every machine but the one that wrote it, and that one only
until the scratch directory was cleaned. The test that loads this fixture
requires a load without warnings, so it passed only there and then. A handle
from `str2func` has no parent file; the regenerated file was verified to load
without a warning after its generator had been deleted.

SHA-256:

- `legacy-point-record.mat`
  `9425143f238085f34aacc57b0c29003f654085043656fdb90cb079bdb1c05089`

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
