"""Numerical merit, reference targets, convergence work and profile scores.

No archive writes, rendering or report publication belongs to this layer.
Arrays retain the benchmark's problem/solver/run/evaluation axis conventions.
"""

import json

import numpy as np

from .utils import ProfileOption


def _default_score_weight_fun(x):
    return 1


def _default_score_fun(x):
    return np.mean(x[:, :, 0, 0], axis=1)


def compute_merit_values(merit_fun, fun_values, maxcv_values, maxcv_init):
    """
    Compute the merit function values.

    This is the Python counterpart of MATLAB's ``meritFunCompute.m`` and
    supports the same family of input shapes. The shape of ``fun_values``
    drives the broadcasting of ``maxcv_init``; the supported combinations
    are:

    +----------------+--------------------------------------+----------------------+
    | use case       | ``fun_values`` shape                 | ``maxcv_init`` shape |
    +================+======================================+======================+
    | merit history  | ``(n_problems, n_solvers, n_runs,    | ``(n_problems,       |
    |                | n_evals)``                           | n_runs)``            |
    +----------------+--------------------------------------+----------------------+
    | merit out      | ``(n_problems, n_solvers, n_runs)``  | ``(n_problems,       |
    |                |                                      | n_runs)``            |
    +----------------+--------------------------------------+----------------------+
    | merit inits    | ``(n_problems, n_runs)``             | ``(n_problems,       |
    |                |                                      | n_runs)``            |
    +----------------+--------------------------------------+----------------------+
    | single-problem | ``(n_solvers, n_runs, n_evals)``     | ``(n_runs,)``        |
    | history        |                                      |                      |
    +----------------+--------------------------------------+----------------------+
    | single-problem | ``(n_solvers, n_runs)``              | ``(n_runs,)``        |
    | out / inits    | or ``(n_runs,)``                     |                      |
    +----------------+--------------------------------------+----------------------+
    | legacy 1-D     | ``(n_problems, n_solvers, n_runs,    | ``(n_problems,)``    |
    | inits (loaded  | n_evals)`` or coarser                |                      |
    | from old .h5)  |                                      |                      |
    +----------------+--------------------------------------+----------------------+
    | scalar init    | any                                  | scalar               |
    +----------------+--------------------------------------+----------------------+

    Parameters
    ----------
    merit_fun : callable
        The merit function. It must have the signature
        ``merit_fun(fun_value, maxcv_value, maxcv_init)`` and accept
        either scalar or numpy array inputs (``np.vectorize`` is used
        internally to lift scalar implementations).
    fun_values : array-like
        The objective function values.
    maxcv_values : array-like
        The maximum constraint violation values; same shape as
        ``fun_values``.
    maxcv_init : float or array-like
        The initial maximum constraint violation value(s); see the
        table above for the supported shapes.

    Returns
    -------
    merit_values : np.ndarray
        The merit function values, with the same shape as ``fun_values``.
    """
    fun_values = np.asarray(fun_values)
    maxcv_values = np.asarray(maxcv_values)
    fun_shape = fun_values.shape

    if np.isscalar(maxcv_init):
        # Scalar broadcast: trivially equal at every entry.
        maxcv_init_b = np.full(fun_shape, maxcv_init)
    else:
        maxcv_init = np.asarray(maxcv_init)
        init_shape = maxcv_init.shape

        if init_shape == fun_shape:
            # Already the target shape (e.g. ``merit_inits`` recomputed
            # from per-run ``fun_inits`` / ``maxcv_inits``).
            maxcv_init_b = maxcv_init
        else:
            ndim_fun = len(fun_shape)
            # Mirror MATLAB ``meritFunCompute.m`` (type='multiple' /
            # 'single'). We pick the case by inspecting the shapes
            # explicitly rather than relying on numpy's right-aligned
            # broadcasting, which would silently mix up the
            # ``n_problems`` / ``n_runs`` axes.
            if ndim_fun == 4 and init_shape == (fun_shape[0], fun_shape[2]):
                # (P, S, R, E) <- (P, R)
                new_shape = (fun_shape[0], 1, fun_shape[2], 1)
            elif ndim_fun == 3 and init_shape == (fun_shape[0], fun_shape[2]):
                # (P, S, R) <- (P, R)
                new_shape = (fun_shape[0], 1, fun_shape[2])
            elif ndim_fun == 3 and init_shape == (fun_shape[1],):
                # (S, R, E) <- (R,)  (single-problem history)
                new_shape = (1, fun_shape[1], 1)
            elif ndim_fun == 2 and init_shape == (fun_shape[1],):
                # (S, R) <- (R,)  (single-problem out / inits)
                new_shape = (1, fun_shape[1])
            elif ndim_fun > 1 and init_shape == (fun_shape[0],):
                # Backward compatibility for legacy 1-D ``maxcv_inits``
                # loaded from old .h5 files: broadcast along the
                # leading axis only.
                new_shape = (fun_shape[0],) + (1,) * (ndim_fun - 1)
            elif ndim_fun == 1 and init_shape == fun_shape:
                new_shape = init_shape
            else:
                raise ValueError(
                    f"compute_merit_values: cannot align maxcv_init of shape "
                    f"{init_shape} with fun_values of shape {fun_shape}."
                )
            maxcv_init_b = np.broadcast_to(maxcv_init.reshape(new_shape), fun_shape)

    # A callback may return an integer first and a fraction or infinity later.
    # Do not let np.vectorize infer an integer/bool dtype from the first result.
    vectorized_merit_fun = np.vectorize(merit_fun, otypes=[float])
    merit_values = vectorized_merit_fun(fun_values, maxcv_values, maxcv_init_b)
    # A custom merit callback must not turn an invalid raw evaluation into a
    # finite reference value. Preserve the callback's policy for +/-Inf.
    invalid = np.isnan(fun_values) | np.isnan(maxcv_values)
    if np.any(invalid):
        merit_values = np.where(invalid, np.inf, merit_values)

    return merit_values


def _default_merit(fun_value, maxcv_value, maxcv_init):
    """
    The merit function varphi(x) is defined by the objective function f(x) and the maximum constraint violation v(x):
        varphi(x) = f(x)                        if v(x) <= v1
        varphi(x) = f(x) + 1e5 * (v(x) - v1)    if v1 < v(x) <= v2
        varphi(x) = np.inf                      if v(x) > v2
    where v1 = min(0.01, 1e-10 * max(1, v0)), v2 = max(0.1, 2 * v0), and v0 is the initial maximum constraint violation.

    Special case 1: _default_merit(np.nan, maxcv_value, maxcv_init) = np.inf
    Sepcial case 2: _default_merit(fun_value, np.nan, maxcv_init) = np.inf
    """
    if np.isnan(fun_value) or np.isnan(maxcv_value):
        return np.inf
    tolerance = 1e-10
    v1 = min(0.01, tolerance * max(1, maxcv_init))
    v2 = max(0.1, 2 * maxcv_init)
    if maxcv_value <= v1:
        return fun_value
    if v1 < maxcv_value <= v2:
        return fun_value + 1e5 * (maxcv_value - v1)
    else:
        return np.inf


def _mask_invalid_merits(results_plib):
    """Invalidate cached merits only where saved raw evaluations are NaN.

    Mutate and return one library record, including its nested plain results.
    Legacy helper records without raw values retain their existing merits.
    """
    if not results_plib:
        return results_plib
    for suffix in ('histories', 'outs', 'inits'):
        merit_key = 'merit_' + suffix
        if merit_key not in results_plib:
            continue
        merits = np.asarray(results_plib[merit_key])
        if suffix == 'inits' and merits.ndim == 1:
            for raw_key in ('fun_inits', 'maxcv_inits'):
                raw = np.asarray(results_plib.get(raw_key, []))
                if raw.ndim == 2:
                    merits = np.broadcast_to(merits[:, None], raw.shape).copy()
                    results_plib[merit_key] = merits
                    break
        invalid = np.zeros(merits.shape, dtype=bool)
        for prefix in ('fun_', 'maxcv_'):
            raw_key = prefix + suffix
            if raw_key not in results_plib:
                continue
            raw_invalid = np.isnan(results_plib[raw_key])
            if suffix == 'inits' and raw_invalid.ndim == 1 and merits.ndim == 2:
                raw_invalid = raw_invalid[:, None]
            invalid |= raw_invalid
        if np.any(invalid):
            results_plib[merit_key] = np.where(invalid, np.inf, merits)
    plain = results_plib.get('results_plib_plain')
    if plain:
        _mask_invalid_merits(plain)
    return results_plib


def reference_target(reference, merit_fun, maxcv_init, problem_name=None, run=None):
    """Convert one validated raw reference fact into the current merit space.

    The record itself was validated when it was constructed. A valid record can
    still be incompatible with the selected merit: the default merit, for
    example, is +Inf when the violation exceeds max(0.1, 2 * maxcv_init). Such a
    target is reported as an error instead of being replaced by the empirical
    target, which would silently change the scoring semantics.
    """
    value = merit_fun(reference.fun, reference.maxcv, maxcv_init)
    scalar = np.asarray(value)
    if scalar.ndim != 0 or scalar.dtype.kind not in 'iuf' or not np.isfinite(scalar):
        where = '' if problem_name is None else f" of problem '{problem_name}'"
        where += '' if run is None else f' (run {run})'
        raise ValueError(
            f'The supplied reference{where} is valid but incompatible with the selected merit function: '
            f'the merit function must produce a finite real scalar for a supplied reference, and it returned '
            f'{value!r} for kind={reference.kind!r}, fun={reference.fun!r}, maxcv={reference.maxcv!r} at '
            f'maxcv_init={float(maxcv_init)!r}. The default merit is infinite when the violation exceeds '
            'max(0.1, 2 * maxcv_init). Use a merit_fun that is finite at this reference, or exclude the problem.')
    return float(scalar)


def apply_reference_targets(result, minima, merit_fun):
    """Replace only supplied targets, converting raw facts with current merit.

    Store raw records, not converted merits: load may change merit_fun, and
    perturbed starts give different maxcv_init in each run. Missing fields in
    V1 archives mean the historical policy, never an assumed zero reference.
    """
    records = result.get('reference_records')
    if records is None:
        return minima
    from .opclasses import ProblemReference
    targets = np.array(minima, dtype=float, copy=True)
    if len(records) != targets.shape[0]:
        raise ValueError('Saved reference records do not match the problem rows.')
    inits = np.asarray(result['maxcv_inits'])
    names = result.get('problem_names')
    for row, serialized in enumerate(records):
        per_run = json.loads(serialized)
        if not isinstance(per_run, list) or len(per_run) != targets.shape[1]:
            raise ValueError('Saved reference records do not match the run count.')
        for run, record in enumerate(per_run):
            if record is None:
                continue
            reference = ProblemReference.from_record(record)
            initial_violation = inits[row] if inits.ndim == 1 else inits[row, run]
            name = None if names is None or row >= len(names) else names[row]
            targets[row, run] = reference_target(reference, merit_fun, initial_violation, name, run + 1)
    return targets


def carries_reference_records(results_plib):
    """Whether a library record contains at least one explicit reference fact."""
    records = results_plib.get('reference_records') if results_plib else None
    if records is None:
        return False
    for serialized in records:
        per_run = json.loads(serialized) if isinstance(serialized, (str, bytes)) else serialized
        if isinstance(per_run, list) and any(record is not None for record in per_run):
            return True
    return False


def compute_convergence_work(merit_histories, merit_outs, merit_inits, merit_mins,
                             n_evals, tolerance, invalid_initial, valid_histories, valid_outs):
    """Return first passing evaluation indices and work at passing solver outputs.

    Validity masks exclude raw NaN observations, padded history tails and
    output fallbacks (runs whose solver raised or returned nothing).
    The +Inf-initial convention and arithmetic order match profile scoring.
    An output reached without any evaluation is charged one unit of work, so
    that work stays positive in performance and data profile ratios.
    """
    n_problems, n_solvers, n_runs = merit_histories.shape[:3]
    initial_is_inf = np.isposinf(merit_inits)
    # Compute the number of function evaluations used by each
    # solver on each problem at each run to achieve convergence.
    work_hist = np.full((n_problems, n_solvers, n_runs), np.nan)
    work_out = np.full((n_problems, n_solvers, n_runs), np.nan)
    for i_problem in range(n_problems):
        for i_solver in range(n_solvers):
            for i_run in range(n_runs):
                if invalid_initial[i_problem, i_run] or np.isnan(merit_inits[i_problem, i_run]):
                    continue
                if initial_is_inf[i_problem, i_run]:
                    # Preserve the +Inf-initial convention, but only for
                    # actual evaluations without raw NaN values. A -Inf
                    # initial merit makes the minimum -Inf: exact hit below.
                    threshold = np.inf
                elif np.isfinite(merit_mins[i_problem, i_run]):
                    threshold = max(tolerance * merit_inits[i_problem, i_run] + (1.0 - tolerance) * merit_mins[i_problem, i_run], merit_mins[i_problem, i_run])
                else:
                    # A -Inf history minimum needs an exact -Inf hit;
                    # an undefined minimum cannot establish a target.
                    threshold = -np.inf
                passing = valid_histories[i_problem, i_solver, i_run] & (merit_histories[i_problem, i_solver, i_run] <= threshold)
                if np.any(passing):
                    work_hist[i_problem, i_solver, i_run] = np.flatnonzero(passing)[0] + 1
                if valid_outs[i_problem, i_solver, i_run] and merit_outs[i_problem, i_solver, i_run] <= threshold:
                    work_out[i_problem, i_solver, i_run] = max(n_evals[i_problem, i_solver, i_run], 1)

    return work_hist, work_out


def process_results(results_plibs, profile_options):
    """
    Processes results from results_plibs.
    Returns:
        merit_histories_merged, merit_outs_merged, merit_inits_merged, merit_mins_merged, n_evals_merged, problem_names_merged, problem_dims_merged
    """
    def unify_length(results_plibs, name):
        # Unify the length of the specified field in results_plibs (merit_histories or merit_histories_plain)
        max_length = max(r[name].shape[3] for r in results_plibs)
        for r in results_plibs:
            field = r[name]
            if field.shape[3] == max_length:
                r[name] = np.asarray(field, dtype=float)
                continue
            shape = list(field.shape)
            shape[3] = max_length
            field_unified = np.full(shape, np.nan)
            field_unified[:, :, :, :field.shape[3]] = field
            field_unified[:, :, :, field.shape[3]:] = field[:, :, :, -1:]
            r[name] = field_unified
        return results_plibs

    # Padding is local to this merge. Keep the caller's merit/raw shapes aligned
    # for saving or processing the same results again; arrays need no deep copy.
    results_plibs = [record.copy() for record in results_plibs]
    for record in results_plibs:
        if record.get('results_plib_plain'):
            record['results_plib_plain'] = record['results_plib_plain'].copy()

    merit_histories_merged = []
    merit_outs_merged = []
    merit_inits_merged = []
    n_evals_merged = []
    problem_names_merged = []
    problem_dims_merged = []

    # Old saved custom-merit arrays may predate raw-NaN validation.
    for results_plib in results_plibs:
        _mask_invalid_merits(results_plib)
    results_plibs = unify_length(results_plibs, 'merit_histories')
    for results_plib in results_plibs:
        merit_histories_merged.append(results_plib['merit_histories'])
        merit_outs_merged.append(results_plib['merit_outs'])
        merit_inits_merged.append(results_plib['merit_inits'])
        n_evals_merged.append(results_plib['n_evals'])
        problem_dims_merged.append(results_plib['problem_dims'])
        problem_names_merged.extend(results_plib['problem_names'])

    merit_histories_merged = np.concatenate(merit_histories_merged, axis=0)
    merit_outs_merged = np.concatenate(merit_outs_merged, axis=0)
    merit_inits_merged = np.concatenate(merit_inits_merged, axis=0)
    n_evals_merged = np.concatenate(n_evals_merged, axis=0)
    problem_dims_merged = np.concatenate(problem_dims_merged, axis=0)

    # Find the least merit value for each (problem, run). Mirrors MATLAB
    # ``processResults.m``:
    #
    #     merit_mins_merged = min(min(history, [], 4, omitnan), [], 2,
    #                              omitnan)
    #
    # i.e. collapse the evaluation axis (axis=3) and the solver axis
    # (axis=1), keeping the (problem, run) layout. ``merit_inits_merged``
    # is itself per-(problem, run) so the elementwise ``np.fmin`` is
    # well-defined and enforces the ``merit_min <= merit_init`` invariant
    # used by the Moré-Wild convergence test downstream.
    #
    # Backward compatibility: if ``merit_inits_merged`` was loaded from
    # an older .h5 file as a 1-D ``(n_problems,)`` vector, broadcast it
    # along the run axis before taking the elementwise min.
    # A custom merit can return NaN even for valid raw evaluations. Retain a
    # valid initial reference when every history merit is NaN, matching MATLAB
    # omitnan. fmin.reduce also leaves an entirely undefined slice as NaN without
    # an all-NaN warning. Raw-invalid initial points remain rejected by the
    # independent convergence-validity mask; this is not a feasibility repair.
    merit_mins_merged = np.fmin.reduce(merit_histories_merged, axis=(1, 3))
    if merit_inits_merged.ndim == 1:
        merit_inits_for_min = merit_inits_merged[:, None]
    else:
        merit_inits_for_min = merit_inits_merged
    merit_mins_merged = np.fmin(merit_mins_merged, merit_inits_for_min)

    # Match plain rows inside their original library container, before losing
    # library/variant identity in the merged arrays. Saved H5 data already keeps
    # this nesting, including when load filters solvers or problem rows.
    if profile_options.get(ProfileOption.RUN_PLAIN, False):
        first_problem = 0
        for results_plib in results_plibs:
            names = results_plib['problem_names']
            offset = first_problem
            first_problem += len(names)
            plain = results_plib.get('results_plib_plain')
            if not plain or not plain['problem_names']:
                continue
            plain_names = plain['problem_names']
            if len(set(names)) != len(names) or len(set(plain_names)) != len(plain_names):
                raise ValueError(
                    'The run_plain reference is ambiguous: duplicate problem names '
                    'within one problem library cannot be matched uniquely.'
                )
            # One reference per problem: reduce over all plain solvers, runs,
            # and evaluations. fmin ignores invalid plain values, as MATLAB's
            # min(..., omitnan) does, rather than erasing a valid featured value.
            plain_mins = np.fmin.reduce(plain['merit_histories'], axis=(1, 2, 3))
            plain_inits = plain['merit_inits']
            # Preserve the first-run convention and legacy one-dimensional inits.
            init_for_plain = plain_inits if plain_inits.ndim == 1 else plain_inits[:, 0]
            plain_mins = np.fmin(plain_mins, init_for_plain)
            plain_indices = {name: idx for idx, name in enumerate(plain_names)}
            for i_problem, name in enumerate(names):
                if name in plain_indices:
                    row = offset + i_problem
                    merit_mins_merged[row, :] = np.fmin(
                        merit_mins_merged[row, :], plain_mins[plain_indices[name]],
                    )

    # An explicit target takes precedence over both the cohort and run_plain.
    # Do not min() it with observed values: that would silently ignore a user
    # target reached or beaten by a solver. Actual run merits remain untouched.
    offset = 0
    for result in results_plibs:
        end = offset + len(result['problem_names'])
        if 'reference_records' in result:
            merit_mins_merged[offset:end] = apply_reference_targets(
                result, merit_mins_merged[offset:end], profile_options[ProfileOption.MERIT_FUN])
        offset = end
    return (merit_histories_merged, merit_outs_merged, merit_inits_merged, merit_mins_merged, n_evals_merged, problem_names_merged, problem_dims_merged)


def integrate(curve, profile_type, profile_options):
    """
    Compute the integral of the curve from different types of profiles.

    Parameters
    ----------
    curve : tuple, list, or numpy.ndarray
        The curve to integrate. Can be:
        - A tuple/list (x, y) where x and y are 1D arrays
        - A 2D array with shape (2, n_points) where curve[0] is x and curve[1] is y
    profile_type : str
        Type of profile: 'perf', 'data', or 'log_ratio'.
    profile_options : dict
        Profile options containing SCORE_WEIGHT_FUN.

    Returns
    -------
    float
        The integral value.
    """
    # Handle None or empty cases
    if curve is None:
        return 0.0
    if isinstance(curve, (tuple, list)):
        if len(curve) == 0:
            return 0.0
        # curve is (x, y) tuple/list
        x = np.asarray(curve[0])
        y = np.asarray(curve[1])
        if x.size == 0 or y.size == 0:
            return 0.0
    elif isinstance(curve, np.ndarray):
        if curve.size == 0:
            return 0.0
        # curve is 2D array with shape (2, n_points)
        x = curve[0, :]
        y = curve[1, :]
    else:
        return 0.0

    kernel = profile_options.get(ProfileOption.SCORE_WEIGHT_FUN, _default_score_weight_fun)
    integral = 0.0

    if profile_type in ('perf', 'data'):
        # The curve is a right-continuous step function.
        # Use rectangle rule (left endpoint): sum(dx * y[:-1] * kernel(x[:-1]))
        if len(x) < 2:
            return 0.0
        dx = np.diff(x)
        y_left = y[:-1]
        x_left = x[:-1]
        with np.errstate(invalid='ignore', divide='ignore', over='ignore'):
            contribution = dx * y_left * kernel(x_left)
        integral = np.sum(contribution[np.isfinite(contribution)])
    elif profile_type == 'log_ratio':
        # We do not modify the integral of log_ratio even if a score_weight_fun is provided.
        abs_y = np.abs(y)
        integral = np.sum(abs_y[np.isfinite(abs_y)])

    return integral


def compute_scores(curves, profile_options):
    """
    Compute the scores of the solvers for all the profiles.

    Parameters
    ----------
    curves : list
        List of curve dictionaries, one per tolerance level.
    profile_options : dict
        Profile options.

    Returns
    -------
    profile_scores : numpy.ndarray
        A 4D array containing scores for all profiles.
        Shape: (n_solvers, n_tols, 2, 2 or 3)
        - First dimension: solver index
        - Second dimension: tolerance index
        - Third dimension: 0 = history-based, 1 = output-based
        - Fourth dimension: 0 = perf, 1 = data, 2 = log_ratio (only if n_solvers == 2)
    """
    if not curves:
        return np.array([])

    n_tols = len(curves)
    n_solvers = len(curves[0]['hist']['perf'])

    if n_solvers == 2:
        profile_scores = np.ones((n_solvers, n_tols, 2, 3))
    else:
        profile_scores = np.ones((n_solvers, n_tols, 2, 2))

    for i_tol in range(n_tols):
        for i_solver in range(n_solvers):
            # Get curves for this solver and tolerance (use the last run, which is the mean)
            curve_hist_perf = curves[i_tol]['hist']['perf'][i_solver][-1]
            curve_hist_data = curves[i_tol]['hist']['data'][i_solver][-1]
            curve_out_perf = curves[i_tol]['out']['perf'][i_solver][-1]
            curve_out_data = curves[i_tol]['out']['data'][i_solver][-1]

            profile_scores[i_solver, i_tol, 0, 0] = integrate(curve_hist_perf, 'perf', profile_options)
            profile_scores[i_solver, i_tol, 0, 1] = integrate(curve_hist_data, 'data', profile_options)
            profile_scores[i_solver, i_tol, 1, 0] = integrate(curve_out_perf, 'perf', profile_options)
            profile_scores[i_solver, i_tol, 1, 1] = integrate(curve_out_data, 'data', profile_options)

            if n_solvers == 2:
                curve_hist_log_ratio = curves[i_tol]['hist']['log_ratio'][i_solver]
                curve_out_log_ratio = curves[i_tol]['out']['log_ratio'][i_solver]
                profile_scores[i_solver, i_tol, 0, 2] = integrate(curve_hist_log_ratio, 'log_ratio', profile_options)
                profile_scores[i_solver, i_tol, 1, 2] = integrate(curve_out_log_ratio, 'log_ratio', profile_options)

    # Normalize the profile_scores by dividing by the maximum with the same tolerance and types
    if profile_options.get(ProfileOption.NORMALIZED_SCORES, False):
        max_scores = np.max(profile_scores, axis=0, keepdims=True)
        max_scores[max_scores == 0] = 1
        profile_scores = profile_scores / max_scores

    return profile_scores
