.. _use:

Usage for MATLAB
================

OptiProfiler provides a :ref:`benchmark <matbenchmark>` function. This is the main entry point to the package. It benchmarks given solvers on the selected test suite.

We provide below simple examples on how to use OptiProfiler in MATLAB. For more details on the signature of the :ref:`benchmark <matbenchmark>` function, please refer to the :ref:`MATLAB API documentation <matlabapi>`.

Examples
--------

.. _example1:

Example 1: first example to try out
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
(See also the file in the repository: ``matlab/examples/example1.m``)

Let us first try to benchmark two callable optimization solvers **solver1** and **solver2** (e.g., **fminsearch** and **fminunc** in MATLAB Optimization Toolbox) on the default test suite.
(Note that each **solver** must accept signatures mentioned in the `Cautions` part of the :ref:`benchmark <matbenchmark>` function according to the type of problems you want to solve.)

To do this, run:

.. code-block:: matlab

    scores = benchmark({@solver1, @solver2})

This will benchmark the two solvers under the default test setting, which means ``'plain'`` feature (see :ref:`Feature <matfeature>`) and unconstrained problems from the default problem library whose dimension is smaller or equal to 2. It will also return the scores of the two solvers based on the profiles.

There will be a new folder named ``out`` in the current working directory, which contains a subfolder named ``plain_<timestamp>`` with all the detailed results.

.. figure:: images/mat_subfolder_structure.png
   :width: 80%
   :align: center
   :alt: Structure of the plain_<timestamp> subfolder

   Figure 1: Screenshot of the subfolder containing detailed results of the benchmarking run.

Additionally, a PDF file named ``summary.pdf`` is generated, summarizing all the performance profiles and data profiles.

.. figure:: images/mat_summary_pdf.jpg
   :width: 90%
   :align: center
   :alt: Summary PDF preview

   Figure 2: Screenshot of the summary PDF file summarizing all the performance profiles and data profiles.

The subfolder ``test_log`` contains diagnostic files for the experiment. In
particular, ``test_log/report.txt`` records selected problem names, timing
information, and special cases detected while building the profiles: problems
where ``merit_init = phi(x_0) = Inf`` (valid evaluations pass by convention,
but undefined objective or constraint values do not), solver runs that
terminated abnormally, and solver outputs that
were replaced by the initial point as an output-based penalty. The file
``test_log/log.txt`` contains the messages printed during the run.

If history plots are enabled, each problem library has a subfolder under
``history_plots``. The top-level ``PROBLEM.pdf`` files keep the combined view
with both raw and cumulative-minimum histories. The ``raw`` and ``cummin``
subfolders contain the corresponding single-view PDFs. The merged
``*_history_plots_summary.pdf`` file uses the top-level combined PDFs.

.. include:: history_display.rst

.. include:: matlab_runtime_output.rst

.. _example2:

Example 2: one step further by adding options
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
(See also the file in the repository: ``matlab/examples/example2.m``)

You can also add options to the benchmark function. For example, if you want to benchmark three solvers **solver1**, **solver2**, and **solver3** on the test suite with the ``'noisy'`` feature and all the unconstrained and bound-constrained problems with dimension between 6 and 10 from the default problem set, you can run:

.. code-block:: matlab

    options.ptype = 'ub';
    options.mindim = 6;
    options.maxdim = 10;
    options.feature_name = 'noisy';
    scores = benchmark({@solver1, @solver2, @solver3}, options)

This will create the corresponding folders ``out/noisy_<timestamp>`` and files as in :ref:`Example 1 <example1>`. More details on the options can be found in the :ref:`benchmark <matbenchmark>` function documentation.

For the deterministic noisy variant from Moré and Wild's benchmarking model,
set ``options.noise_mode = 'deterministic'``. If ``options.n_runs`` is not
provided, OptiProfiler uses one run for this deterministic feature unless
``options.solver_isrand`` marks at least one solver as randomized, in which
case OptiProfiler uses five runs as usual.

.. code-block:: matlab

    options.feature_name = 'noisy';
    options.noise_mode = 'deterministic';
    options.noise_map = 'chebyshev';
    scores = benchmark({@solver1, @solver2, @solver3}, options)

By default, **n_jobs** is set conservatively to about half of the available
workers instead of all workers. For the most reproducible timing experiments,
set ``options.n_jobs`` explicitly, for example ``options.n_jobs = 1`` for
sequential runs.


Example 3: useful option **load**
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
(See also the file in the repository: ``matlab/examples/example3.m``)

OptiProfiler provides a practically useful option named **load**. This option allows you to load the results from a previous benchmarking run (without solving all the problems again) and use them to draw new profiles with different options. For example, if you have just run :ref:`Example 2 <example2>` and OptiProfiler has finished the job and successfully created the folder ``out`` in the current working directory, you can run:

.. code-block:: matlab

    options.load = 'latest';
    options.solvers_to_load = [1, 3];
    options.ptype = 'u';
    options.mindim = 7;
    options.maxdim = 9;
    scores = benchmark(options)

This will directly draw the profiles for the **solver1** and **solver3** with the ``'noisy'`` feature and all the unconstrained problems with dimension between 7 and 9 selected from the previous run. The results will also be saved under the current directory with a new subfolder named ``noisy_<timestamp>`` with the new timestamp.


Reusable features, replay, and retained results
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

A ``Feature`` is a reusable specification, not an experiment or a random-stream
checkpoint. Its effective stages contain only feature-local options. For example:

.. code-block:: matlab

    stages = { ...
        struct('name', 'noisy', 'options', struct('noise_level', 1e-3)), ...
        struct('name', 'perturbed_x0', 'options', struct('perturbation_level', 1e-3))};
    options = struct('feature', Feature(stages), 'n_runs', 3, 'n_jobs', 1);
    scores = benchmark({@solver1, @solver2}, options);

Pass either ``feature`` or ``feature_name``, not both. ``n_runs`` belongs to the
benchmark options, not to a stage. Reusing the same Feature in another benchmark
does not reuse its previous runtime state. Identity stages such as ``plain`` are
normalized away. Repeated nonidentity stages retain their order and separate
options. With ``run_plain=true``, the primary experiment uses its own run count;
the plain-reference role retains the independent one-run policy. A saved run
axis is not a count of actual solver calls: deterministic runs can have copied
slots, which reports identify separately from actual executions.

There are two different ways to use an existing experiment:

* ``benchmark`` with ``load`` reanalyses saved numerical histories without
  executing solvers or loading the original provider problems. Current filters
  describe this selection; they are not a recovered original experiment plan.
  MATLAB keeps its existing sorted-unique ``solvers_to_load`` selection rule.
* ``loadBenchmarkOptions`` imports trusted native settings for a **new**
  benchmark. The caller supplies its solvers and output destination explicitly:

  .. code-block:: matlab

      [options, receipt] = loadBenchmarkOptions('test_log/options_refined.mat');
      options.benchmark_id = 'fresh-replay';
      options.savepath = pwd;
      scores = benchmark({@solver1, @solver2}, options);

New ``options_refined-v2`` files retain the canonical native Feature and a
separate ``n_runs``. The effective ``feature_specification`` is inspection-only
when the native Feature is present; inconsistent ordinary values, structure,
or callback positions/types are rejected. Native anonymous functions and
closures need not compare equal after a MAT round trip, so callback identity
and semantics are not compared. The returned receipt explicitly records this
limit, and only the canonical Feature supplies executable callbacks.
Canonical native state preserves resolved options and a separately retained
declaration, including an unknown declaration after a supported legacy import.
It does not reconstruct an original invocation route or resume a live runtime.
Saved output paths, load selectors and envelope fields are not forwarded by the
replay helper. Specification-only v2 structs are accepted as explicit new replay
inputs, not evidence of the original declaration.

Options saved by a ``load`` (re-plot) invocation describe that re-plot, not the
archived execution: ``benchmark`` writes no ``options_refined.mat`` for a load,
and ``loadBenchmarkOptions`` rejects an unversioned source with a nonempty
``load`` field (the ``options_user.mat`` of a load, or a flat
``options_refined.mat`` that older versions wrote for a load) with
``OptiProfiler:LoadInvocationNotReplayable`` instead of turning the load label
into a replay. A versioned ``options_refined-v2`` record carries its own
effective feature identity, so it is accepted as an explicit replay input and
its obsolete ``load``/``solvers_to_load`` selectors are simply not forwarded.
Replay the source experiment from its own ``test_log/options_refined.mat``.

Old flat native options sometimes omitted feature identity. Such a file requires
an explicit second argument, for example ``loadBenchmarkOptions(path, 'noisy')``;
neither folder names nor stamps are used to guess it. The receipt records this
as a present replay override, not a recovered historical request. Native MAT
files are trusted inputs and require their original callback/class dependencies;
they can execute code during loading. JSON callback descriptions cannot be used
to reconstruct executable functions.

New evaluation reports use ``optiprofiler.eval_report/2`` and archive provenance
uses ``feature_pipeline-v3``. Effective feature settings and per-role experiment
plans are separate. Historical payloads retain their original version and
unknown facts; a load report does not turn today's defaults into old execution
facts. JSON stage positions and same-kind occurrence numbers are zero-based;
MATLAB numerical solver/run indexes remain one-based as recorded by the report.
Complete observed runtime receipts stay in native archives, while compact
reports omit repeated per-stage runtime details. Long or private feature text
may be omitted from JSON with explicit UTF-8 byte counts and reasons; complete
native stamps remain available in native archives/settings.

Report feature declarations and effective stage lists longer than 256 entries
use a bounded object with ``values`` (the first 256 entries), ``total_items``,
and ``reason='metadata_item_limit'``. This is a report projection, not a limit
on feature execution or on the complete native Feature. Ordinary numerical
and problem-result arrays are not truncated by this feature-metadata rule.

Long generated display labels use a bounded prefix and CRC-32 suffix in MATLAB;
Python uses SHA-256 for that display suffix. These labels do not change seeds,
scientific feature identity, or output-directory uniqueness. Equal shortened
folder names across languages are not promised.

Feature option values are checked when a ``Feature`` is built. The numeric
options ``perturbation_level``, ``noise_level``, ``condition_factor``,
``mesh_size``, ``nan_rate`` and ``significant_digits`` must be finite real
scalars within their documented bounds; MATLAB stores them as double, and an
integer-class value must be exactly representable as double.
``perturbation_level`` is a nonnegative scalar: MATLAB has no coordinatewise
vector amplitudes. Logical options accept ``true``/``false`` or numeric 0/1 and
are stored as logical. Text choices such as ``mesh_type`` are lowercase and
case-sensitive. Option names are case-insensitive; two spellings of one name
(``noise_level`` and ``NOISE_LEVEL``) in one struct or name/value list are
rejected (``MATLAB:Feature:DuplicateOption``) because the configuration would
be ambiguous. Earlier versions accepted some values that could not describe a
valid experiment: NaN or infinite magnitudes, any ``perturbation_level``
(including negative, text or vector values), and integer or single classes that
then computed in integer or single arithmetic. Saved Features are revalidated
when executed again, including through ``loadBenchmarkOptions``. Invalid values
are rejected explicitly; supported integer and single scalars are instead
canonicalized to double, just like fresh input. ``benchmark`` with ``load`` still reanalyses
the saved histories, because it does not execute the saved Feature.

On the identity and single-stage strategies, ``FeaturedProblem`` now counts
nonlinear-constraint evaluations once per recorded query, as the composed
strategy and the Python package already did. Earlier versions reported the
number of constraint components while fewer queries had been recorded. They
could therefore serve stale constraint values once that number reached the
evaluation budget, and they passed repeated indices to stochastic constraint
modifiers. New executions record ``runtime_policy='matlab-legacy-single-v2'``.
Archives written under ``matlab-legacy-single-v1`` keep their recorded policy
when they are loaded. Numerical results are unchanged for constraint channels
with at most one component.

Composed features (two or more effective stages) record
``seed_policy='matlab-stage-horner32-v2'``. The per-stage, per-channel seeds
are derived exactly as under version 1 (an exact 32-bit Horner fold over the
run seed and the stage identity); version 2 seeds every per-query stream of a
composed view by folding the IEEE-754 words of the observed payload (values,
point, served index) with the same rule. Version 1 handed that payload to the
legacy product mixer, so a zero coordinate, a zero value or a zero counter
removed the dependence on the rest of the payload. Archives written under
version 1 keep their recorded policy string. The identity and single-stage
strategies are unchanged and keep their established legacy streams. The fold
is a finite 32-bit hash: it is not a statistical independence guarantee and
distinct payloads can still collide.

Migration note: ``is_stochastic`` is a read-only property of ``Feature`` (as
are ``name``, ``stages`` and ``declared``); the earlier method call syntax
``is_stochastic(F)`` is no longer available because MATLAB cannot expose one
name as both a property and a method. Use ``F.is_stochastic``.


Example 4: testing parametrized solvers
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
(See also the file in the repository: ``matlab/examples/example4.m``)

If you want to benchmark a solver with one variable parameter, you can define function handles by looping over the parameter values. For example, if **solver** accepts the signature ``@(fun, x0, para)``, and you want to benchmark it with the parameter ``para`` taking values from 1 to 3, you can run:

.. code-block:: matlab

    solvers = cell(1, 3);
    options.solver_names = cell(1, 3);
    for i = 1:3
        solvers{i} = @(fun, x0) solver(fun, x0, i);
        options.solver_names{i} = ['solver' num2str(i)];
    end
    scores = benchmark(solvers, options)


Example 5: customizing the test suite
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
(See also the file in the repository: ``matlab/examples/example5.m``)

OptiProfiler allows you to customize the test suite by creating your own feature and loading your own problem library.
For example, if you want to create a new feature that adds noise to the objective function and perturbs the initial guess at the same time, you can try the following:

.. code-block:: matlab

    options.feature_name = 'custom';
    options.mod_fun = @(x, rand_stream, problem) problem.fun(x) + 1e-3 * rand_stream.randn(1);
    options.mod_x0 = @(rand_stream, problem) problem.x0 + 1e-3 * rand_stream.randn(problem.n, 1);
    scores = benchmark({@solver1, @solver2}, options)

Installed providers are selected with ``options.plibs``.  See
:ref:`matlab_problem_libraries` for installing MatCUTEst and SOLAR, and
:ref:`matlab_custom_problem_libraries` to register a library stored in any
user-owned directory.  The detailed provider-authoring contract is kept out of
this solver tutorial so that the benchmark workflow remains focused.


Example 6: wrapping solvers with nonlinear constraints
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
(See also the file in the repository: ``matlab/examples/example6.m``)

For nonlinearly constrained problems, OptiProfiler calls each solver with the
signature

.. code-block:: matlab

    x = solver(fun, x0, xl, xu, aub, bub, aeq, beq, cub, ceq)

where ``cub(x) <= 0`` contains the nonlinear inequality constraints and
``ceq(x) = 0`` contains the nonlinear equality constraints. MATLAB solvers
such as ``fmincon`` instead expect one nonlinear constraint callback
``nonlcon`` returning both values. The small but important adapter is
``deal``: the expression ``@(x) deal(cub(x), ceq(x))`` evaluates
OptiProfiler's two callbacks and returns them as the two outputs expected by
``fmincon``, i.e., ``[c, ceq] = nonlcon(x)``. See the MathWorks
documentation for
`fmincon nonlinear constraints <https://www.mathworks.com/help/optim/ug/nonlinear-constraints.html>`_
and `deal <https://www.mathworks.com/help/matlab/ref/deal.html>`_.

.. code-block:: matlab

    function x = fmincon_short(fun, x0, xl, xu, aub, bub, aeq, beq, cub, ceq)
        x = fmincon_wrapper(fun, x0, xl, xu, aub, bub, aeq, beq, cub, ceq, 100);
    end

    function x = fmincon_long(fun, x0, xl, xu, aub, bub, aeq, beq, cub, ceq)
        x = fmincon_wrapper(fun, x0, xl, xu, aub, bub, aeq, beq, cub, ceq, 200);
    end

    function x = fmincon_wrapper(fun, x0, xl, xu, aub, bub, aeq, beq, cub, ceq, max_fun_evals)
        % Convert OptiProfiler's separate nonlinear callbacks to fmincon's
        % two-output callback: [c, ceq] = nonlcon(x).
        nonlcon = @(x) deal(cub(x), ceq(x));
        options = optimoptions('fmincon', 'MaxFunctionEvaluations', max_fun_evals);
        x = fmincon(fun, x0, aub, bub, aeq, beq, xl, xu, nonlcon, options);
    end

Then pass the wrappers to ``benchmark`` as ordinary solvers:

.. code-block:: matlab

    options.ptype = 'n';
    options.problem_names = {'HS10', 'HS11', 'HS12'};
    options.plibs = {'s2mpj'};
    options.mindim = 2;
    options.maxdim = 5;
    options.max_eval_factor = 500;
    options.draw_hist_plots = 'none';
    options.n_jobs = 1;
    options.solver_names = {'fmincon short', 'fmincon long'};
    scores = benchmark({@fmincon_short, @fmincon_long}, options)
