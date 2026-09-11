.. _use_python:

Usage for Python
================

OptiProfiler provides a :func:`~optiprofiler.benchmark` function. This is the main entry point to the package. It benchmarks given solvers on the selected test suite.

We provide below simple examples on how to use OptiProfiler in Python. For more details on the signature of the :func:`~optiprofiler.benchmark` function, please refer to the :ref:`Python API documentation <pythonapi>`.

Examples
--------

.. _py_example1:

Example 1: first example to try out
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

Let us first try to benchmark two callable optimization solvers **solver1** and **solver2** on the default test suite.
(Note that each **solver** must accept signatures mentioned in the `Cautions` part of the :func:`~optiprofiler.benchmark` function according to the type of problems you want to solve.)

To do this, run:

.. code-block:: python

    from optiprofiler import benchmark

    scores = benchmark([solver1, solver2])

This will benchmark the two solvers under the default test setting, which means ``'plain'`` feature (see :class:`~optiprofiler.Feature`) and unconstrained problems from the default problem library whose dimension is smaller or equal to 2. It will also return the scores of the two solvers based on the profiles.

There will be a new folder named ``out`` in the current working directory, which contains a subfolder named ``plain_<timestamp>`` with all the detailed results.

.. figure:: images/py_subfolder_structure.png
   :width: 80%
   :align: center
   :alt: Structure of the plain_<timestamp> subfolder

   Figure 1: Screenshot of the subfolder containing detailed results of the benchmarking run.

Additionally, a PDF file named ``summary.pdf`` is generated, summarizing all the performance profiles and data profiles.

.. figure:: images/py_summary_pdf.jpg
   :width: 90%
   :align: center
   :alt: Summary PDF preview

   Figure 2: Screenshot of the summary PDF file summarizing all the performance profiles and data profiles.

The subfolder ``test_log`` contains diagnostic files for the experiment. In
particular, ``test_log/report.txt`` records selected problem names, timing
information, and special cases detected while building the profiles: problems
where ``merit_init = phi(x_0) = inf`` (all solvers are declared to pass that
problem/run), solver runs that terminated abnormally, and solver outputs that
were replaced by the initial point as an output-based penalty. The file
``test_log/log.txt`` contains the messages printed during the run.

If history plots are enabled, each problem library has a subfolder under
``history_plots``. The top-level ``PROBLEM.pdf`` files keep the combined view
with both raw and cumulative-minimum histories. The ``raw`` and ``cummin``
subfolders contain the corresponding single-view PDFs. The merged
``*_history_plots_summary.pdf`` file uses the top-level combined PDFs.

.. include:: history_display.rst

.. _py_example2:

Example 2: one step further by adding options
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

You can also add options to the benchmark function. For example, if you want to benchmark three solvers **solver1**, **solver2**, and **solver3** on the test suite with the ``'noisy'`` feature and all the unconstrained and bound-constrained problems with dimension between 6 and 10 from the default problem set, you can run:

.. code-block:: python

    from optiprofiler import benchmark

    scores = benchmark(
        [solver1, solver2, solver3],
        ptype='ub',
        mindim=6,
        maxdim=10,
        feature_name='noisy',
    )

This will create the corresponding folders ``out/noisy_<timestamp>`` and files as in :ref:`Example 1 <py_example1>`. More details on the options can be found in the :func:`~optiprofiler.benchmark` function documentation.

For the deterministic noisy variant from Moré and Wild's benchmarking model,
set ``noise_mode='deterministic'``. If ``n_runs`` is not provided, OptiProfiler
uses one run for this deterministic feature unless ``solver_isrand`` marks at
least one solver as randomized, in which case OptiProfiler uses five runs as
usual.

.. code-block:: python

    scores = benchmark(
        [solver1, solver2, solver3],
        feature_name='noisy',
        noise_mode='deterministic',
        noise_map='chebyshev',
    )

By default, ``n_jobs`` is set conservatively to about half of the available
workers instead of all workers. For the most reproducible timing experiments,
set ``n_jobs`` explicitly, for example ``n_jobs=1`` for sequential runs.

.. _py_composing_features:

Composing features
^^^^^^^^^^^^^^^^^^

Several features can be applied in order by joining their names with ``+``.
The first name is applied to the original problem first: ``'noisy+truncated'``
adds noise to every objective value and then truncates the noisy value to
``significant_digits`` digits, whereas ``'truncated+noisy'`` truncates first and
adds noise afterwards. Any number of stages is allowed, names may repeat, and
``'plain'`` is an identity that can appear anywhere.

.. code-block:: python

    scores = benchmark(
        [solver1, solver2],
        feature_name='perturbed_x0+noisy+truncated',
        perturbation_level=1e-2,
        noise_level=1e-4,
        significant_digits=4,
    )

Order matters for features that move the query point as well. With
``feature_name='quantized+linearly_transformed+quantized'`` a solver point ``x``
is first snapped to the mesh by the last stage, then mapped through the linear
transformation, and snapped again by the first stage before the original
objective is evaluated. For this chain, an objective query reads the original
objective once at ``snap(A @ snap(x))`` for the value the solver observes and
once more for the scoring reference recorded in the history; each constraint
channel behaves the same way. That count is a property of this chain of lazy
point maps and linear transport, not a general promise. The general rule is
only that a composition performs no unused or exponentially repeated reads:
each stage reads its predecessor exactly as its operation requires, so chains
of point maps stay linear in length. Additional reads of the original problem
are legitimate wherever a stage asks for them: a ``custom`` callback that
probes its predecessor, the ``unrelaxable_constraints`` gate reading the
predecessor's constraints during an objective query, and the reference reads
for the initial point, the histories and the final scoring.

The rules are:

- Every supplied feature option is passed to every stage that accepts it, and
  each stage validates it with its own rules. ``feature_name='noisy+perturbed_x0'``
  therefore uses Gaussian noise and a spherical perturbation by default,
  ``distribution='gaussian'`` configures both stages, and
  ``distribution='uniform'`` is rejected with an error naming the
  ``perturbed_x0`` stage. Options accepted by no stage are rejected. Repeated
  stages share the supplied options; configuring two occurrences of the same
  feature differently is not expressible through the flat option list.
- ``n_runs`` is an experiment option, never a feature option:
  ``benchmark(..., n_runs=N)`` sets it, the experiment plan records it once,
  and no stage carries a run count (``Feature(..., n_runs=N)`` is an error
  naming the benchmark option). Unless given explicitly, it is 5 when
  ``solver_isrand`` marks a randomized solver (not in a load) and otherwise
  the largest established default of the effective stages. Those defaults
  are the ones of the single features: five
  for ``perturbed_x0``, ``noisy`` (one with ``noise_mode='deterministic'``),
  ``permuted``, ``linearly_transformed`` (also with ``rotated=False``,
  although that variant is deterministic), ``random_nan`` and ``truncated``
  with ``perturbed_trailing_digits=True``; one for ``plain``, ``truncated``,
  ``unrelaxable_constraints``, ``nonquantifiable_constraints``,
  ``quantized`` and ``custom`` (although ``custom`` counts as stochastic).
  So ``'noisy+truncated'`` uses five runs and ``'truncated+quantized'`` uses
  one. The composition counts as stochastic when any stage is stochastic.
- Value changes (``noisy``, ``truncated``, ``random_nan``,
  ``nonquantifiable_constraints``) affect what solvers observe, never the
  scoring reference recorded in the histories. ``permuted`` and
  ``linearly_transformed`` transport both. ``perturbed_x0`` only moves the
  initial point. ``quantized`` with ``ground_truth=False`` snaps observations
  only; with ``ground_truth=True`` the inherited reference objective and
  nonlinear constraints are read at the snapped point as well, while bounds and
  linear constraints are always checked at the unsnapped point.
- ``unrelaxable_constraints`` makes the objective infinite where the
  constraints of the problem it wraps are violated, as that problem observes
  them: after ``'noisy+unrelaxable_constraints'`` the gate reads noisy
  constraint samples, after ``'unrelaxable_constraints+noisy'`` it reads the
  original constraints. After a rotation, former bounds are linear constraints
  of the wrapped problem and belong to the linear category. Violations that
  contain NaN never close the gate, as for the single feature.
- Custom callbacks of a ``'custom'`` stage receive the problem produced by the
  preceding stages, so ``problem.fun(x)`` inside a callback is a genuine query
  of that problem. Each call is a separately served query: a predecessor whose
  observations are randomized per query (such as ``noisy`` or ``random_nan``)
  draws its own sample for it, while randomness fixed when the problem is
  built (``perturbed_x0``, ``permuted``, ``linearly_transformed``) and a
  deterministic predecessor return the same value again. Callback outputs are
  validated at the custom stage
  (objective values follow the ``Problem.fun`` scalar policy and are recorded
  as NaN with a logged warning when they are not real scalars; constraint
  outputs must be real one-dimensional arrays of the predecessor's size,
  otherwise an error names the stage and callback).
- A name whose only effective stage is a single feature, such as
  ``'plain+noisy'``, behaves exactly like that feature, including its random
  streams and its output folder name. A genuine composition seeds every
  channel of every stage separately from the run seed (policy
  ``seedsequence-v2``): the seed is the first 32-bit word of
  ``numpy.random.SeedSequence(run_seed, spawn_key=(code, occurrence, tag))``
  with a frozen numeric code per feature, the occurrence index among stages of
  the same name, and the tags ``fun=0``, ``cub=1``, ``ceq=2`` and
  ``construction=3`` (draws made when the problem is built). Inserting
  ``'plain'`` therefore does not change the derivation identity of any other
  stage, and distinct identities remove the structural alias that omitting the
  channel tag would cause between the objective and constraint channels. This
  is not an independence proof of the generator, and 32-bit seeds can still
  coincide. Archives record the policy identifier; archives written under an
  earlier policy keep their own identifier when loaded.
- The output folder name joins the stamps of the stages with ``__`` and is
  shortened with a digest for long chains. The archive and the structured
  report record the declared specification, the effective stages with their
  local options and the seed policy, and separately the experiment plan (see
  :ref:`py_structured_feature` for the payload). Derivatives of a composed
  problem are not provided.

.. _py_structured_feature:

Structured feature specification
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

The ``feature_name`` shorthand broadcasts every supplied option to all stages
that accept it, so two occurrences of the same feature, or a shared key such
as ``distribution`` for ``noisy`` and ``perturbed_x0``, cannot be configured
differently. The ``feature`` option gives each stage its own options. It is
one stage entry ``{'name': ..., 'options': {...}}`` or an ordered list (or
tuple) of entries; a bare name is a stage with default options, ``options``
may be omitted, and ``'plain'`` entries are identities.

.. code-block:: python

    scores, profile_scores, curves = benchmark(
        [solver1, solver2],
        n_runs=3,
        feature=[
            {'name': 'noisy', 'options': {'distribution': 'uniform', 'noise_level': 1e-2}},
            {'name': 'perturbed_x0', 'options': {'distribution': 'gaussian'}},
            {'name': 'noisy', 'options': {'noise_level': 1e-4, 'noise_type': 'absolute'}},
        ],
    )

    # The same experiment as feature_name='noisy', noise_level=1e-2:
    scores, _, _ = benchmark([solver1, solver2], feature={'name': 'noisy', 'options': {'noise_level': 1e-2}})

The rules are:

- Option ownership is explicit. ``n_runs`` belongs to the experiment and is
  given at the top level; inside an entry it is rejected. Every other option
  belongs to the stage named in its entry and is validated with that stage's
  rules, so ``{'name': 'perturbed_x0', 'options': {'noise_level': 1e-2}}`` is
  an error naming the entry and the stage. Profile options such as ``seed``
  or ``max_eval_factor`` are never stage options.
- One route per call. ``feature`` and ``feature_name`` cannot both be given
  (an explicit ``feature_name='plain'`` counts as given), ``feature=None``
  and an empty specification are errors, a string belongs to
  ``feature_name``, and with ``feature`` every flat stage option is rejected:
  there is no override hierarchy. These errors are raised before any output
  directory is created; a requested ``report_path`` still records the failure.
- Every entry, ``'plain'`` ones included, is validated before ``'plain'``
  entries are removed, so an invalid option on an identity stage is never
  dropped silently. Names are stripped and lowercased, and stage names are
  atomic (``{'name': 'noisy+truncated'}`` is an error). A specification with
  one effective stage runs the established single-feature path, with its
  seeds and folder name.
- Equal settings give equal numbers. ``feature=[{'name': 'noisy', 'options':
  {'noise_level': 1e-2}}, 'truncated']`` and ``feature_name='noisy+truncated',
  noise_level=1e-2`` produce the same histories, outputs, folder stamp and
  scores for the same seed. Stage identities stay name plus occurrence, so
  differently configured repeats keep their seeds and differ only by the
  configured amounts.
- ``feature`` cannot be combined with ``load``. A specification means
  transformations to execute, and loading keeps the archived pipeline of the
  saved experiment. Experiments created with ``feature`` load like any other,
  with their archived pipeline retained; ``feature_name`` keeps labelling a
  load as before.

Feature objects. ``benchmark(feature=...)`` also accepts an already built
:class:`~optiprofiler.Feature`, and both routes build one. A ``Feature`` is an
immutable pipeline specification: ``feature.stages`` is the tuple of effective
stage records, each with ``name``, ``occurrence``, ``identity``
(``'noisy#0'``), ``code`` (the frozen seed code), ``position`` and ``options``
(a read-only view of the validated stage-local options after defaults:
arrays come back as isolated read-only copies and lists/tuples as fresh
copies, while callables and other opaque values keep their identity);
``feature.name`` is the effective name (``'plain'`` for the identity, which
has no stage), ``feature.declared`` the declaration as given (route and entries
before defaults, with the same ownership rule) and ``feature.is_stochastic``
whether any stage draws random numbers. Ownership of option values is
explicit: exact built-in data containers (``list``, ``tuple``, ``dict``) and
NumPy arrays are copied when the specification is built and the views return
isolated copies, so changing such a value afterwards, or what a view returned,
does not change the specification. Opaque values are kept by identity and are
not copied: callables, container subclasses, other objects and the elements
of object-dtype arrays; the specification never mutates them, but their own
state remains yours to keep unchanged. No copy, iteration or reduction hook
of your objects is ever invoked. The declaration route (``'feature_name'`` for the shorthand string,
``'feature'`` for structured entries) is a property of the specification: a
``Feature`` declared by shorthand and passed as ``benchmark(feature=obj)``
keeps its declaration and its receipt. The provenance records that route as
``declaration_route`` and, separately, the benchmark keyword that carried the
specification as ``route``. A specification imported from a historical object
keeps the declaration that object recorded (``plain`` entries included); one
imported from an object that recorded no declaration has declaration route
``None``, no declared entries and no declared name. In both cases its
effective stages are complete. Each trial builds fresh
runtime state from the records, so reusing
one ``Feature`` object across benchmarks gives the same numbers as building it
again. Pickling a ``Feature`` writes its versioned native form: the effective
stages with every validated local option explicit (so the defaults of the
version that loads it cannot change the numbers) and, separately, the
declaration as provenance; loading validates the effective options again,
callables included as the native objects, and rejects a native form version
it does not know. Objects pickled by OptiProfiler 1.x are read only through
the trusted decoder described below; pickles of other layouts or from other
versions are not a supported input. ``Feature.options``
and the ``modifier_*`` methods remain as deprecated one-stage conveniences
(``DeprecationWarning``; an error on a composition) and are used by nothing
in the engine, the archives or the reports.

Provenance and replay. The archive entry ``feature_pipeline`` (payload
``feature_pipeline-v3``) has two blocks. The ``feature`` block records the
benchmark keyword that carried the specification (``route``), the declaration
route of the specification (``declaration_route``), the declared name, the
declared entries with the options as supplied, the effective name, the seed
policy, the stamps with their origin (explicit or generated) and, for every effective
stage, its position, name, code, occurrence, identity and local options
(callables are described by name, never executed). The ``experiment`` block
records the resolved plan of the primary role: ``n_runs`` and where it came
from (``explicit``, ``randomized_solvers`` or ``stage_hints``), the run
policy, the execution strategy and the language-local runtime identifier; run
counts appear there and nowhere inside a stage. The structured report shows the same feature block under
``configuration.effective.feature`` and the plans of every executing role
under ``configuration.effective.experiment``. Archives written by earlier
versions (``feature_pipeline-v1`` and ``-v2``, or none at all) keep their own
payload verbatim when loaded; nothing is migrated or relabelled. The file
``test_log/options_refined.pkl`` (``schema`` ``options_refined-v2``) stores
``feature_route``, ``feature_name`` (the declared name), ``n_runs`` and
``feature_specification``: the ordered effective stages with their validated
options as native Python values, callables included, and no flat stage keys.
That entry is valid ``feature`` input, so
``scores, _, _ = benchmark(feature=refined['feature_specification'], n_runs=refined['n_runs'], ...)``
reproduces the effective experiment for both routes; forwarding the whole
refined dictionary is not a supported call. ``optiprofiler.legacy_compat``
is the trusted boundary for files written by earlier versions:
``load_options`` decodes historical enumeration members (including the former
``n_runs`` feature option), ``replay_arguments`` maps the supported layouts to
``feature`` and ``n_runs`` (a flat 1.x file records no feature identity and
needs an explicit ``feature_name``), and ``import_legacy_feature`` converts a
1.x pickled ``Feature`` into a canonical ``Feature`` plus its retained run
count, which is passed to ``benchmark`` separately. Only load pickle and H5
files from trusted sources: unpickling can execute code, and these readers
offer no safety for untrusted payloads. The callback descriptions in the
archive and the report are one-way informational metadata (a class or function
name), never a recipe for reconstructing an executable callable.

MATLAB mapping. The MATLAB implementation follows the same contract:
``options.feature`` as a cell array of structs
(``struct('name', 'noisy', 'options', struct('noise_level', 1e-3))``) or
names, ``options.n_runs`` at the top level, the same ``feature_pipeline-v3``
fields and report version, and a language-local seed policy for compositions
(``matlab-stage-horner32-v1``); random samples are not matched across the two
languages. ``options.feature_name`` keeps its current meaning. See the MATLAB
user guide for the native replay helper ``loadBenchmarkOptions``.

.. _py_example3:

Example 3: useful option **load**
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

.. warning::

   Load only experiment archives from trusted sources. Python HDF5 archives
   may contain pickled objects; loading an untrusted file can execute code.
   The loader is not a safe reader for public uploads.

OptiProfiler provides a practically useful option named ``load``. This option allows you to load the results from a previous benchmarking run (without solving all the problems again) and use them to draw new profiles with different options. For example, if you have just run :ref:`Example 2 <py_example2>` and OptiProfiler has finished the job and successfully created the folder ``out`` in the current working directory, you can run:

.. code-block:: python

    from optiprofiler import benchmark

    scores = benchmark(
        load='latest',
        solvers_to_load=[0, 2],
        ptype='u',
        mindim=7,
        maxdim=9,
    )

This will directly draw the profiles for the **solver1** and **solver3** with the ``'noisy'`` feature and all the unconstrained problems with dimension between 7 and 9 selected from the previous run. The results will also be saved under the current directory with a new subfolder named ``noisy_<timestamp>`` with the new timestamp.

.. _py_example4:

Example 4: testing parametrized solvers
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

If you want to benchmark a solver with one variable parameter, you can define callables by looping over the parameter values. For example, if **solver** accepts the signature ``solver(fun, x0, para)``, and you want to benchmark it with the parameter ``para`` taking values from 1 to 3, you can run:

.. code-block:: python

    from optiprofiler import benchmark

    def make_solver(para):
        def solver_wrapper(fun, x0):
            return solver(fun, x0, para)
        return solver_wrapper

    solvers = [make_solver(i) for i in range(1, 4)]
    solver_names = [f'solver{i}' for i in range(1, 4)]
    scores = benchmark(solvers, solver_names=solver_names)

.. note::

    We use named functions (``def``) instead of lambda expressions here so
    that the benchmark can still run in parallel when ``n_jobs > 1``.
    See :ref:`py_callable_picklability` for the full list of affected
    callables and the rationale.

.. _py_example5:

Example 5: customizing the test suite
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

OptiProfiler allows you to customize the test suite by creating your own feature and loading your own problem library.
For example, if you want to create a new feature that adds noise to the objective function and perturbs the initial guess at the same time, you can try the following:

.. code-block:: python

    from optiprofiler import benchmark

    def mod_fun(x, rand_stream, problem):
        return problem.fun(x) + 1e-3 * rand_stream.standard_normal()

    def mod_x0(rand_stream, problem):
        return problem.x0 + 1e-3 * rand_stream.standard_normal(problem.n)

    scores = benchmark(
        [solver1, solver2],
        feature_name='custom',
        mod_fun=mod_fun,
        mod_x0=mod_x0,
    )

.. note::

    Again, ``mod_fun`` and ``mod_x0`` are defined with ``def`` rather than
    ``lambda`` so that the benchmark can run in parallel when
    ``n_jobs > 1``. See :ref:`py_callable_picklability` for details.

Problem libraries may be distributed as independent Python packages. See
:ref:`python_problem_libraries` for installation and lifecycle management and
:ref:`python_problem_library_providers` for the supported providers. Once a
compatible package is installed, its public library name can be passed directly
to ``plibs``; no filesystem path is needed. For example, a package named
``optiprofiler-myproblems`` may register the library ``'myproblems'``, after
which it can be used as follows:

.. code-block:: python

    scores = benchmark(
        [solver1, solver2],
        plibs=['s2mpj', 'myproblems'],
    )

OptiProfiler discovers these packages without importing them. Their optional
dependencies are imported only when the corresponding library is selected for
a benchmark. See :ref:`python_custom_problem_libraries` to create a local or
installable provider.

Library-specific experiment options
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Different problem libraries may expose entirely different options. Keep these
options separated by library name with ``plib_options`` instead of mixing them
with the common problem filters such as ``ptype``, ``mindim``, and ``maxdim``.
For example, using two hypothetical installable plugins:

.. code-block:: python

    scores = benchmark(
        [solver1, solver2],
        plibs=['myproblems', 'otherproblems'],
        plib_options={
            'myproblems': {
                'variant': 'large',
            },
            'otherproblems': {
                'data_split': 'validation',
            },
        },
    )

Each library validates only its own mapping. Options for a library not listed
in ``plibs`` and unknown option names are errors. ``plib_options`` is only for
runs that select problems from libraries; it is rejected with ``problem`` and
``load`` because those modes do not call a library adapter. The precedence is:

1. values in ``benchmark(..., plib_options=...)``;
2. process-level overrides from :func:`~optiprofiler.set_plib_config`;
3. environment or package configuration owned by the library;
4. the library's built-in defaults.

Use :func:`~optiprofiler.get_plib_config` to inspect the effective
process-level values. The per-experiment values do not mutate them. The raw
``plib_options`` supplied by the user are saved in ``options_user.pkl``. The
validated mapping with defaults filled in is saved in ``options_refined.pkl``
and, as one type-preserving serialized mapping, in the library's group in
``data_for_loading.h5``. These files are under the experiment's ``test_log``
directory.

Paths to external runtimes, licenses, compiled binaries, and caches are
installation concerns rather than benchmark semantics. A library should check
those in its availability hook instead of placing them in ``plib_options``.

For a local filesystem adapter or an independently distributed plugin, follow
:ref:`python_custom_problem_libraries`. The provider-authoring contract is kept
out of this quick-start tutorial so that ordinary benchmark usage remains
focused.

.. _py_example6:

Example 6: wrapping SciPy solvers with nonlinear constraints
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

(See also the file in the repository: ``python/examples/scipy_cobyqa_wrapper.py``)

For nonlinearly constrained problems, OptiProfiler calls each solver with the
signature

.. code-block:: python

    x = solver(fun, x0, xl, xu, aub, bub, aeq, beq, cub, ceq)

where ``cub(x) <= 0`` contains the nonlinear inequality constraints and
``ceq(x) = 0`` contains the nonlinear equality constraints. SciPy's
``minimize`` interface represents constraints with objects such as
``Bounds``, ``LinearConstraint``, and ``NonlinearConstraint``. The adapter in
the wrapper below is the conversion from OptiProfiler's callback signature to
SciPy's constraint objects: linear constraints become ``LinearConstraint``
objects, while ``cub`` and ``ceq`` are wrapped as ``NonlinearConstraint``
objects with bounds ``(-inf, 0)`` and ``(0, 0)``, respectively. The SciPy
documentation for `COBYQA <https://docs.scipy.org/doc/scipy/reference/optimize.minimize-cobyqa.html>`_
and the `optimization tutorial <https://docs.scipy.org/doc/scipy/tutorial/optimize.html>`_
show this object-based constraint interface; see also the API references for
`LinearConstraint <https://docs.scipy.org/doc/scipy/reference/generated/scipy.optimize.LinearConstraint.html>`_
and `NonlinearConstraint <https://docs.scipy.org/doc/scipy/reference/generated/scipy.optimize.NonlinearConstraint.html>`_.

.. code-block:: python

    import numpy as np
    from scipy.optimize import Bounds, LinearConstraint, NonlinearConstraint, minimize

    def scipy_cobyqa_wrapper(fun, x0, xl, xu, aub, bub, aeq, beq, cub, ceq,
                             maxfev=200):
        constraints = []

        if bub.size > 0:
            # OptiProfiler gives aub @ x <= bub; SciPy stores it as a
            # LinearConstraint with lower bound -inf and upper bound bub.
            constraints.append(LinearConstraint(aub, -np.inf, bub))
        if beq.size > 0:
            # Equality constraints use identical lower and upper bounds.
            constraints.append(LinearConstraint(aeq, beq, beq))

        c_ub_x0 = np.atleast_1d(cub(x0))
        if c_ub_x0.size > 0:
            # Convert cub(x) <= 0 to a SciPy NonlinearConstraint.
            constraints.append(NonlinearConstraint(cub, -np.inf, np.zeros_like(c_ub_x0)))
        c_eq_x0 = np.atleast_1d(ceq(x0))
        if c_eq_x0.size > 0:
            # Convert ceq(x) = 0 by using zero lower and upper bounds.
            constraints.append(NonlinearConstraint(ceq, np.zeros_like(c_eq_x0),
                                                   np.zeros_like(c_eq_x0)))

        result = minimize(
            fun,
            x0,
            method='COBYQA',
            bounds=Bounds(xl, xu),
            constraints=constraints,
            options={'maxfev': maxfev},
        )
        return result.x

Then pass the wrapper to ``benchmark`` as an ordinary solver. Since
``benchmark`` compares at least two solvers, this example compares two
COBYQA wrappers with different function-evaluation budgets:

.. code-block:: python

    def scipy_cobyqa_short(fun, x0, xl, xu, aub, bub, aeq, beq, cub, ceq):
        return scipy_cobyqa_wrapper(
            fun, x0, xl, xu, aub, bub, aeq, beq, cub, ceq, maxfev=100
        )

    def scipy_cobyqa_long(fun, x0, xl, xu, aub, bub, aeq, beq, cub, ceq):
        return scipy_cobyqa_wrapper(
            fun, x0, xl, xu, aub, bub, aeq, beq, cub, ceq, maxfev=200
        )

    scores = benchmark(
        [scipy_cobyqa_short, scipy_cobyqa_long],
        solver_names=['SciPy COBYQA short', 'SciPy COBYQA long'],
        ptype='n',
        problem_names=['HS10', 'HS11', 'HS12'],
        mindim=2,
        maxdim=5,
        max_eval_factor=500,
        plibs=['s2mpj'],
        draw_hist_plots='none',
        n_jobs=1,
    )

.. _py_cautions:

Cautions
--------

.. _py_callable_picklability:

Callable arguments must be picklable when running in parallel
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

When ``n_jobs > 1``, OptiProfiler dispatches problems to worker
processes via :mod:`multiprocessing`. The following callable
arguments are sent across process boundaries and must therefore be
picklable:

- the entries of ``solvers``;
- feature options: ``distribution``, ``noise_map``, ``mod_x0``,
  ``mod_affine``, ``mod_bounds``, ``mod_linear_ub``, ``mod_linear_eq``,
  ``mod_fun``, ``mod_cub``, ``mod_ceq``;
- profile options: ``merit_fun``, ``score_fun``, ``score_weight_fun``.

**Lambda expressions and locally-defined nested functions are not
picklable.** If any of the callables above is a lambda, OptiProfiler
detects the failure when serializing the worker arguments and silently
falls back to sequential mode (``n_jobs = 1``), which can be much
slower.

To enable parallel execution, define these callables as module-level
functions using ``def``. For parametrized solvers, use a closure
factory (see :ref:`py_example4`) or :func:`functools.partial` instead
of a lambda.
