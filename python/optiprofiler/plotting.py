import os
import re
import shutil
import inspect
import time
import warnings
from contextlib import redirect_stderr, redirect_stdout
from datetime import datetime
from inspect import signature
from multiprocessing import Pool
from pathlib import Path

import matplotlib
import numpy as np
from cycler import cycler
from matplotlib.figure import Figure
from matplotlib.colors import is_color_like
from matplotlib.lines import Line2D
from matplotlib.backends import backend_pdf
from matplotlib.ticker import MaxNLocator, FuncFormatter, FixedLocator, NullLocator

from .utils import ProfileOption


_COMPACT_LEGEND_SOLVER_THRESHOLD = 10
_COMPACT_LEGEND_MIN_ROWS_PER_COLUMN = 10
_PROFILE_LABELSIZE = 12
_PROFILE_TICK_LABELSIZE = 11
_HISTORY_LABELSIZE = 12
_HISTORY_TICK_LABELSIZE = 11
_HISTORY_LEGEND_FONTSIZE = 10
_HISTORY_YLABEL_PAD = 18
# A display limit, never an oracle/score/data limit. Squared deviations, run
# means and axis margins retain ample float64 headroom below overflow.
_HISTORY_DISPLAY_LIMIT = 1e100


def latex_escape_text(text):
    """Escape literal labels for Matplotlib when LaTeX rendering is enabled."""
    replacements = {
        '\\': r'\textbackslash{}',
        '_': r'\_',
        '#': r'\#',
        '%': r'\%',
        '&': r'\&',
        '$': r'\$',
        '{': r'\{',
        '}': r'\}',
        '^': r'\^{}',
        '~': r'\~{}',
    }
    return ''.join(replacements.get(char, char) for char in str(text))


def format_profile_text(text, profile_context):
    """Format literal plot text according to the active Matplotlib text backend."""
    if profile_context.get('text.usetex', False):
        return latex_escape_text(text)
    return str(text)


def draw_profiles(work, problem_dimensions, solver_names, tolerance_latex, i_tol, ax_summary_perf, ax_summary_data, ax_summary_log_ratio, is_summary, is_perf, is_data, is_log_ratio, profile_options, curves, *, _plot_sink=None):
    # Escaping is a LaTeX concern: without usetex the escapes would show up as
    # literal backslashes in the legends (e.g. 'scipy\_nelder\_mead').
    profile_context = set_profile_context(profile_options)
    solver_names = [format_profile_text(name, profile_context) for name in solver_names]
    n_solvers = work.shape[1]

    # Calculate the performance and data profiles.
    x_perf, y_perf, ratio_max_perf, x_data, y_data, ratio_max_data, curves = _get_extended_performances_data_profile_axes(work, problem_dimensions, profile_options, curves)
    if n_solvers == 2:
        x_log_ratio, y_log_ratio, ratio_max_log_ratio, n_solvers_fail, curves = _get_log_ratio_profile_axes(work, curves)

    if _plot_sink is not None:
        # Observe the exact arrays already used below, including full bands
        # and the zero/failure bars absent from the scoring-curve subsets.
        # Field names and x_transform labels are shared with MATLAB
        # drawProfiles.m/profilePresentation and pinned by plot_data.schema.json.
        semilogx = profile_options[ProfileOption.SEMILOGX]
        _plot_sink('performance', {
            'series': prepare_profile_plot_data(x_perf, y_perf, profile_options),
            'ratio_max': ratio_max_perf, 'failure_placeholder': 1.1 * ratio_max_perf,
            'n_runs': int(y_perf.shape[2]),
            'x_limits': [0.0 if semilogx else 1.0,
                         1.1 * (ratio_max_perf if np.isfinite(ratio_max_perf) else np.finfo(float).eps)],
            'y_limits': [0.0, 1.0],
            'x_transform': 'log2(work/best_work)' if semilogx else 'work/best_work'})
        _plot_sink('data', {
            'series': prepare_profile_plot_data(x_data, y_data, profile_options),
            'ratio_max': ratio_max_data, 'failure_placeholder': 1.1 * ratio_max_data,
            'n_runs': int(y_data.shape[2]),
            'x_limits': [0.0, 1.1 * (ratio_max_data if np.isfinite(ratio_max_data) else np.finfo(float).eps)],
            'y_limits': [0.0, 1.0],
            'x_transform': 'log2(1+work/(dimension+1))' if semilogx else 'work/(dimension+1)'})
        if n_solvers == 2:
            _plot_sink('log_ratio', prepare_log_ratio_plot_data(
                x_log_ratio, y_log_ratio, ratio_max_log_ratio, n_solvers_fail))
    
    if profile_options[ProfileOption.SCORE_ONLY]:
        return None, None, None, curves

    # File-only figures do not need pyplot's backend or GUI figure managers.
    fig_perf, fig_data = Figure(), Figure()
    ax_perf, ax_data = fig_perf.subplots(), fig_data.subplots()
    if n_solvers > 2:
        fig_log_ratio, ax_log_ratio = None, None
    else:
        fig_log_ratio = Figure()
        ax_log_ratio = fig_log_ratio.subplots()

    # Draw the performance profiles
    _draw_perf_detail(ax_perf, x_perf, y_perf, ratio_max_perf, solver_names, profile_options, tolerance_latex)
    _draw_data_detail(ax_data, x_data, y_data, ratio_max_data, solver_names, profile_options, tolerance_latex)
    if n_solvers == 2:
        _draw_log_ratio_detail(ax_log_ratio, x_log_ratio, y_log_ratio, ratio_max_log_ratio, n_solvers_fail, solver_names, profile_options, tolerance_latex)
    

    # Create the figures in the summary.

    if is_summary:
        if is_perf and ax_summary_perf is not None:
            _draw_perf_detail(ax_summary_perf[i_tol], x_perf, y_perf, ratio_max_perf, solver_names, profile_options, tolerance_latex)
        if is_data and ax_summary_data is not None:
            _draw_data_detail(ax_summary_data[i_tol], x_data, y_data, ratio_max_data, solver_names, profile_options, tolerance_latex)
        if is_log_ratio and ax_summary_log_ratio is not None:
            _draw_log_ratio_detail(ax_summary_log_ratio[i_tol], x_log_ratio, y_log_ratio, ratio_max_log_ratio, n_solvers_fail, solver_names, profile_options, tolerance_latex)

    return fig_perf, fig_data, fig_log_ratio, curves


def _draw_perf_detail(ax_perf, x_perf, y_perf, ratio_max_perf, solver_names, profile_options, tolerance_latex):
    handles = _draw_performance_data_profiles(ax_perf, x_perf, y_perf, solver_names, profile_options)
    ax_perf.tick_params(axis='both', which='major', labelsize=_PROFILE_TICK_LABELSIZE)
    # MATLAB tolerates Inf in set(ax, 'XLim', ...) but matplotlib does not.
    if not np.isfinite(ratio_max_perf):
        ratio_max_perf = np.finfo(float).eps
    # Set x-axis limits.
    if profile_options[ProfileOption.SEMILOGX]:
        ax_perf.set_xlim(0.0, 1.1 * ratio_max_perf)
    else:
        ax_perf.set_xlim(1.0, 1.1 * ratio_max_perf)
    # Modify x-axis ticks labels of the performance profiles (matching MATLAB).
    ticks, tick_labels = perf_ticks(1.1 * ratio_max_perf, profile_options[ProfileOption.SEMILOGX])
    ax_perf.set_xticks(ticks)
    ax_perf.set_xticklabels(tick_labels)
    # Set x-axis labels.
    if profile_options[ProfileOption.XLABEL_PERFORMANCE_PROFILE]:
        ax_perf.set_xlabel(profile_options[ProfileOption.XLABEL_PERFORMANCE_PROFILE], fontsize=_PROFILE_LABELSIZE)
    # Set y-axis labels.
    if profile_options[ProfileOption.YLABEL_PERFORMANCE_PROFILE]:
        ylabel_str = profile_options[ProfileOption.YLABEL_PERFORMANCE_PROFILE] % tolerance_latex if '%s' in profile_options[ProfileOption.YLABEL_PERFORMANCE_PROFILE] else profile_options[ProfileOption.YLABEL_PERFORMANCE_PROFILE]
        ax_perf.set_ylabel(ylabel_str, fontsize=_PROFILE_LABELSIZE, labelpad=10)
    _place_solver_legend(ax_perf, x_perf.shape[1], default_loc='lower right', handles=handles)




def _draw_data_detail(ax_data, x_data, y_data, ratio_max_data, solver_names, profile_options, tolerance_latex):
    handles = _draw_performance_data_profiles(ax_data, x_data, y_data, solver_names, profile_options)
    ax_data.tick_params(axis='both', which='major', labelsize=_PROFILE_TICK_LABELSIZE)
    # MATLAB tolerates Inf in set(ax, 'XLim', ...) but matplotlib does not.
    if not np.isfinite(ratio_max_data):
        ratio_max_data = np.finfo(float).eps
    # Set x-axis limits.
    ax_data.set_xlim(0.0, 1.1 * ratio_max_data)
    # Modify x-axis ticks labels of the data profiles (matching MATLAB).
    ticks, tick_labels = data_ticks(1.1 * ratio_max_data, profile_options[ProfileOption.SEMILOGX])
    ax_data.set_xticks(ticks)
    ax_data.set_xticklabels(tick_labels)
    # Set x-axis labels.
    if profile_options[ProfileOption.XLABEL_DATA_PROFILE]:
        ax_data.set_xlabel(profile_options[ProfileOption.XLABEL_DATA_PROFILE], fontsize=_PROFILE_LABELSIZE)
    # Set y-axis labels.
    if profile_options[ProfileOption.YLABEL_DATA_PROFILE]:
        ylabel_str = profile_options[ProfileOption.YLABEL_DATA_PROFILE] % tolerance_latex if '%s' in profile_options[ProfileOption.YLABEL_DATA_PROFILE] else profile_options[ProfileOption.YLABEL_DATA_PROFILE]
        ax_data.set_ylabel(ylabel_str, fontsize=_PROFILE_LABELSIZE, labelpad=10)
    _place_solver_legend(ax_data, x_data.shape[1], default_loc='lower right', handles=handles)



def _draw_log_ratio_detail(ax_log_ratio, x_log_ratio, y_log_ratio, ratio_max_log_ratio, n_solvers_fail, solver_names, profile_options, tolerance_latex):
    _draw_log_ratio_profiles(ax_log_ratio, x_log_ratio, y_log_ratio, ratio_max_log_ratio, n_solvers_fail, solver_names, profile_options)
    ax_log_ratio.tick_params(axis='both', which='major', labelsize=_PROFILE_TICK_LABELSIZE)
    # Set x-axis labels.
    if profile_options[ProfileOption.XLABEL_LOG_RATIO_PROFILE]:
        ax_log_ratio.set_xlabel(profile_options[ProfileOption.XLABEL_LOG_RATIO_PROFILE], fontsize=_PROFILE_LABELSIZE)
    # Set y-axis labels.
    if profile_options[ProfileOption.YLABEL_LOG_RATIO_PROFILE]:
        ylabel_str = profile_options[ProfileOption.YLABEL_LOG_RATIO_PROFILE] % tolerance_latex if '%s' in profile_options[ProfileOption.YLABEL_LOG_RATIO_PROFILE] else profile_options[ProfileOption.YLABEL_LOG_RATIO_PROFILE]
        ax_log_ratio.set_ylabel(ylabel_str, fontsize=_PROFILE_LABELSIZE, labelpad=14)







def prepare_profile_plot_data(x, y, profile_options):
    """Full step vertices and across-run bands shared by report and renderer."""
    n_solvers = x.shape[1]
    y_mean = np.mean(y, 2)

    if profile_options[ProfileOption.ERRORBAR_TYPE] == 'minmax':
        y_lower = np.min(y, 2)
        y_upper = np.max(y, 2)
    elif profile_options[ProfileOption.ERRORBAR_TYPE] == 'meanstd':
        y_std = np.std(y, 2)
        y_lower = np.maximum(0.0, y_mean - y_std)
        y_upper = np.minimum(1.0, y_mean + y_std)
    else:
        raise ValueError("Unknown {ProfileOption.ERRORBAR_TYPE.value}: {profile_options[ProfileOption.ERRORBAR_TYPE]}")

    return [{'solver_index': i + 1,
             'x': np.repeat(x[:, i], 2)[1:],
             'mean': np.repeat(y_mean[:, i], 2)[:-1],
             'lower': np.repeat(y_lower[:, i], 2)[:-1],
             'upper': np.repeat(y_upper[:, i], 2)[:-1],
             'band_visible': y.shape[2] > 1,
             'geometry': 'step'} for i in range(n_solvers)]


def _draw_performance_data_profiles(ax, x, y, solver_names, profile_options):
    handles = []
    profile_context = set_profile_context(profile_options)
    line_colors = profile_options[ProfileOption.LINE_COLORS]
    line_styles = profile_options[ProfileOption.LINE_STYLES]
    line_widths = profile_options[ProfileOption.LINE_WIDTHS]
    series = prepare_profile_plot_data(x, y, profile_options)

    with matplotlib.rc_context(profile_context):
        for i_solver, item in enumerate(series):
            x_stairs, y_mean_stairs = item['x'], item['mean']
            y_lower_stairs, y_upper_stairs = item['lower'], item['upper']

            color = line_colors[i_solver % len(line_colors)]
            line_style = line_styles[i_solver % len(line_styles)]
            line_width = line_widths[i_solver % len(line_widths)]

            handles.extend(ax.plot(x_stairs, y_mean_stairs, label=solver_names[i_solver], color=color, linestyle=line_style, linewidth=line_width))
            if item['band_visible']:
                ax.fill_between(x_stairs, y_lower_stairs, y_upper_stairs, color=color, alpha=0.2)
        
        # Set Y-axis ticks to match MATLAB style: major ticks at 0, 0.2, 0.4, 0.6, 0.8, 1.0
        ax.set_ylim(0.0, 1.0)
        ax.set_yticks([0.0, 0.2, 0.4, 0.6, 0.8, 1.0])
        # Set minor ticks at 0.1, 0.3, 0.5, 0.7, 0.9
        ax.set_yticks([0.1, 0.3, 0.5, 0.7, 0.9], minor=True)
        # Show ticks on both sides (left and right)
        ax.yaxis.set_ticks_position('both')
        ax.tick_params(axis='y', which='both', direction='in', right=True)
        ax.tick_params(axis='x', which='both', direction='in')
        # Remove tick labels on the right side
        ax.tick_params(axis='y', which='both', labelleft=True, labelright=False)
    return handles


def prepare_log_ratio_plot_data(x, y, ratio_max, n_solvers_equal):
    """Describe all sorted bar positions, including invisible ties and sentinels.

    The historic renderer has no problem identity after sorting. Do not assign
    anonymous bars to problems; target_work retains the unsorted identities.
    """
    opacity = np.ones(len(x))
    if n_solvers_equal:
        opacity[:n_solvers_equal] = 0.5
        opacity[-n_solvers_equal:] = 0.5
    n_below = int(np.sum(y < 0) - n_solvers_equal)
    n_above = int(np.sum(y > 0) - n_solvers_equal)
    groups = []
    if n_solvers_equal > 0:
        groups.extend([{'start_index': 1, 'end_index': int(n_solvers_equal), 'solver_index': 1, 'opacity': 0.5},
                       {'start_index': len(x) - int(n_solvers_equal) + 1, 'end_index': len(x), 'solver_index': 2, 'opacity': 0.5}])
    if n_below > 0:
        groups.append({'start_index': int(n_solvers_equal) + 1, 'end_index': int(n_solvers_equal) + n_below,
                       'solver_index': 1, 'opacity': 1.0})
    if n_above > 0:
        groups.append({'start_index': len(x) - int(n_solvers_equal) - n_above + 1,
                       'end_index': len(x) - int(n_solvers_equal), 'solver_index': 2, 'opacity': 1.0})
    # Shared log-ratio record vocabulary (MATLAB drawProfiles.m emits the same
    # keys). MATLAB retains bar identity (problem_mapping='bar_sources');
    # this renderer sorts anonymously, so problem_mapping is null with a
    # reason rather than an identity guessed from sorted values.
    return {'series': [{'geometry': 'bar', 'x': x, 'y': y,
                        'visible': y != 0, 'opacity': opacity,
                        'solver_index_by_sign': {'negative': 1, 'positive': 2}}],
            'bar_groups': groups,
            'ratio_max': ratio_max, 'x_transform': 'sorted_bar_position',
            'x_limits': [0.5, len(x) + 0.5],
            'y_limits': [-1.1 * ratio_max, 1.1 * ratio_max],
            'y_transform': 'log2(work_solver_1/work_solver_2)',
            'problem_mapping': None, 'problem_mapping_reason': 'legacy_bar_sort_does_not_retain_identity',
            'both_failed_pairs': int(n_solvers_equal),
            'tie_pairs': int(np.count_nonzero(y == 0)),
            'failure_placeholder': 1.1 * ratio_max}

def _draw_log_ratio_profiles(ax, x, y, ratio_max, n_solvers_equal, solver_names, profile_options):
    profile_context = set_profile_context(profile_options)
    n_problems = x.shape[0]
    if not np.isfinite(ratio_max) or ratio_max <= 0:
        ratio_max = np.finfo(float).eps
    presentation = prepare_log_ratio_plot_data(x, y, ratio_max, n_solvers_equal)

    with matplotlib.rc_context(profile_context):
        bar_colors = profile_options[ProfileOption.BAR_COLORS][:2]
        for group in presentation['bar_groups']:
            selection = slice(group['start_index'] - 1, group['end_index'])
            opacity = {'alpha': 0.5} if group['opacity'] == 0.5 else {}
            ax.bar(x[selection], y[selection], color=bar_colors[group['solver_index'] - 1], **opacity)

        # Remove x-axis ticks (matching MATLAB's xticks([]))
        ax.set_xticks([])
        with warnings.catch_warnings():
            warnings.filterwarnings('ignore', category=UserWarning)
            ax.set_xlim(0.5, n_problems + 0.5)
        ax.set_ylim(-1.1 * ratio_max, 1.1 * ratio_max)
        ax.text((n_problems + 1) / 2, -ratio_max, solver_names[0], horizontalalignment='center', verticalalignment='bottom', fontsize=24)
        ax.text((n_problems + 1) / 2, ratio_max, solver_names[1], horizontalalignment='center', verticalalignment='top', fontsize=24)
        # Show ticks on both sides (left and right) and set direction to 'in'
        ax.yaxis.set_ticks_position('both')
        ax.tick_params(axis='y', which='both', direction='in', right=True, labelleft=True, labelright=False)
        ax.tick_params(axis='x', which='both', direction='in')


def _get_extended_performances_data_profile_axes(work, problem_dimensions, profile_options, curves):
    n_problems, n_solvers, n_runs = work.shape

    def denominator_perf(i_problem, i_run):
        return np.nanmin(work[i_problem, :, i_run], initial=np.inf)
    x_perf, y_perf, ratio_max_perf = _get_performance_data_profile_axes(work, denominator_perf)
    if profile_options[ProfileOption.SEMILOGX]:
        # We output the log2(x_perf) and log2(ratio_max_perf). This is because we want the x-axis is in log2 scale in the performance profile.
        x_perf[np.isfinite(x_perf)] = np.log2(x_perf[np.isfinite(x_perf)])
        if ratio_max_perf > np.finfo(float).eps:
            ratio_max_perf = max(np.log2(ratio_max_perf), np.finfo(float).eps)
    x_perf[np.isinf(x_perf)] = 1.1 * ratio_max_perf
    if profile_options[ProfileOption.SEMILOGX]:
        x_perf = np.vstack([np.zeros((1, n_solvers)), x_perf])
    else:
        x_perf = np.vstack([np.ones((1, n_solvers)), x_perf])
    y_perf = np.vstack([np.zeros((1, n_solvers, n_runs)), y_perf])
    if n_problems > 0:
        x_perf = np.vstack([x_perf, np.full((1, n_solvers), ratio_max_perf * 1.1)])
        y_perf = np.vstack([y_perf, y_perf[-1, np.newaxis, :, :]])

    def denominator_data(i_problem, i_run):
        return problem_dimensions[i_problem] + 1
    x_data, y_data, ratio_max_data = _get_performance_data_profile_axes(work, denominator_data)
    if profile_options[ProfileOption.SEMILOGX]:
        # We output the log2(1 + x_data) and log2(1 + ratio_max_data). This is because we want the x-axis is in log2 scale in the data profile.
        x_data[np.isfinite(x_data)] = np.log2(1.0 + x_data[np.isfinite(x_data)])
        if ratio_max_data > np.finfo(float).eps:
            ratio_max_data = max(np.log2(1.0 + ratio_max_data), np.finfo(float).eps)
    x_data[np.isinf(x_data)] = 1.1 * ratio_max_data
    x_data = np.vstack([np.zeros((1, n_solvers)), x_data])
    y_data = np.vstack([np.zeros((1, n_solvers, n_runs)), y_data])
    if n_problems > 0:
        x_data = np.vstack([x_data, np.full((1, n_solvers), ratio_max_data * 1.1)])
        y_data = np.vstack([y_data, y_data[-1, np.newaxis, :, :]])

    # Store the curves in the `curves` dictionary.
    for i_solver in range(n_solvers):
        for i_run in range(n_runs):
            curves['perf'][i_solver][i_run] = (x_perf[:, i_solver], y_perf[:, i_solver, i_run])
            curves['data'][i_solver][i_run] = (x_data[:, i_solver], y_data[:, i_solver, i_run])
        y_mean_perf = np.mean(y_perf[:, i_solver, :], axis=1)
        y_mean_data = np.mean(y_data[:, i_solver, :], axis=1)
        curves['perf'][i_solver][n_runs] = (x_perf[:, i_solver], y_mean_perf)
        curves['data'][i_solver][n_runs] = (x_data[:, i_solver], y_mean_data)

    return x_perf, y_perf, ratio_max_perf, x_data, y_data, ratio_max_data, curves


def _get_performance_data_profile_axes(work, denominator):
    """
    Calculate the axes of the performance and data profiles.
    """
    n_problems, n_solvers, n_runs = work.shape

    # Calculate the x-axis values.
    x = np.full((n_runs, n_problems, n_solvers), np.nan)
    for i_run in range(n_runs):
        for i_problem in range(n_problems):
            x[i_run, i_problem, :] = work[i_problem, :, i_run] / denominator(i_problem, i_run)
    ratio_max = np.nanmax(x, initial=np.finfo(float).eps)
    x[np.isnan(x)] = np.inf
    x = np.sort(x, 1)
    x = np.reshape(x, (n_problems * n_runs, n_solvers))
    sort_x = np.argsort(x, 0, 'stable')
    x = np.take_along_axis(x, sort_x, 0)

    # Find the index of the element in x[:, i_solver] that is the last element that is less than or equal to ratio_max.
    index_ratio_max = np.full(n_solvers, np.nan)
    for i_solver in range(n_solvers):
        mask = x[:, i_solver] <= ratio_max
        if np.any(mask):
            index_ratio_max[i_solver] = np.max(np.where(mask)[0])

    # Calculate the y-axis values.
    y = np.full((n_problems * n_runs, n_solvers, n_runs), np.nan)
    for i_solver in range(n_solvers):
        for i_run in range(n_runs):
            if n_problems > 0:
                y[i_run * n_problems:(i_run + 1) * n_problems, i_solver, i_run] = np.linspace(1 / n_problems, 1.0, n_problems)
                y[:, i_solver, i_run] = np.take_along_axis(y[:, i_solver, i_run], sort_x[:, i_solver], 0)
            for i_problem in range(n_problems * n_runs):
                if np.isnan(y[i_problem, i_solver, i_run]):
                    y[i_problem, i_solver, i_run] = y[i_problem - 1, i_solver, i_run] if i_problem > 0 else 0.0

    # Calculate the max_ratio_y, which is the maximum y-value before the corresponding x-value exceeds ratio_max.
    ratio_max_y = np.zeros((n_solvers, n_runs))
    for i_solver in range(n_solvers):
        for i_run in range(n_runs):
            if not np.isnan(index_ratio_max[i_solver]):
                ratio_max_y[i_solver, i_run] = y[int(index_ratio_max[i_solver]), i_solver, i_run]

    # Correct the y-values using the ratio_max_y.
    for i_solver in range(n_solvers):
        for i_run in range(n_runs):
            for i_problem in range(n_problems * n_runs):
                y[i_problem, i_solver, i_run] = min(y[i_problem, i_solver, i_run], ratio_max_y[i_solver, i_run])

    return x, y, ratio_max


def _get_log_ratio_profile_axes(work, curves):
    """
    Calculate the axes of the log-ratio profiles.
    """
    n_problems, n_solvers, n_runs = work.shape
    work_flat = np.reshape(np.swapaxes(work, 1, 2), (n_problems * n_runs, n_solvers))
    y_log_ratio = np.full(n_problems * n_runs, np.nan)
    log_ratio_finite = (
        np.isfinite(work_flat[:, 0])
        & np.isfinite(work_flat[:, 1])
        & (work_flat[:, 0] > 0)
        & (work_flat[:, 1] > 0)
    )
    y_log_ratio[log_ratio_finite] = np.log2(work_flat[log_ratio_finite, 0]) - np.log2(work_flat[log_ratio_finite, 1])

    ratio_max = np.nanmax(np.abs(y_log_ratio[log_ratio_finite]), initial=np.finfo(float).eps)
    if ratio_max == 0 or not np.isfinite(ratio_max):
        ratio_max = np.finfo(float).eps
    
    fail1 = ~np.isfinite(work_flat[:, 0]) | (work_flat[:, 0] <= 0)
    fail2 = ~np.isfinite(work_flat[:, 1]) | (work_flat[:, 1] <= 0)

    mask1 = fail1 & ~fail2
    y_log_ratio[mask1] = 1.1 * ratio_max
    mask2 = ~fail1 & fail2
    y_log_ratio[mask2] = -1.1 * ratio_max

    mask_both_fail = fail1 & fail2
    n_solvers_fail = np.sum(mask_both_fail)
    y_log_ratio[mask_both_fail] = 1.1 * ratio_max
    if n_solvers_fail > 0:
        y_log_ratio = np.concatenate([y_log_ratio, -1.1 * ratio_max * np.ones(n_solvers_fail)])

    y_log_ratio_sorted = np.sort(y_log_ratio)
    x_log_ratio = np.arange(1, len(y_log_ratio_sorted) + 1)

    # Store the curves in the `curves` dictionary.
    curves['log_ratio'][0] = (x_log_ratio[y_log_ratio_sorted < 0], y_log_ratio_sorted[y_log_ratio_sorted < 0])
    curves['log_ratio'][1] = (x_log_ratio[y_log_ratio_sorted > 0], y_log_ratio_sorted[y_log_ratio_sorted > 0])

    return x_log_ratio, y_log_ratio_sorted, ratio_max, n_solvers_fail, curves


def _perf_formatter(x, _):
    if x.is_integer():
        return str(int(2 ** x))
    else:
        return f'$2^{{{f"{x:.8f}".rstrip("0").rstrip(".")}}}$'


def _data_formatter(x, _):
    if x.is_integer():
        return str(int(2 ** x - 1))
    else:
        return f'$2^{{{f"{x:.8f}".rstrip("0").rstrip(".")}}}-1$'


def _format_tick_value(value):
    if np.isclose(value, round(value)):
        return str(int(round(value)))
    return f'{value:.2f}'.rstrip('0').rstrip('.')


def perf_ticks(ratio_cut_perf, is_semilogx):
    """
    Generate ticks and tick labels for performance profiles.
    Matches MATLAB's perfTicks function.
    """
    if is_semilogx:
        if ratio_cut_perf >= 5:
            max_power = int(np.floor(ratio_cut_perf))
            ticks = np.linspace(0, max_power, 6)
            ticks[1:-1] = np.round(ticks[1:-1])
            ticks = np.unique(ticks)
        elif ratio_cut_perf >= 1:
            max_power = int(np.floor(ratio_cut_perf))
            ticks = np.arange(0, max_power + 1)
        elif ratio_cut_perf >= 1e-3:
            ticks = np.array([0, ratio_cut_perf])
        else:
            ticks = np.array([0])
        tick_labels = [_format_tick_value(2 ** t) for t in ticks]
    else:
        if ratio_cut_perf >= 5:
            max_power = int(np.floor(ratio_cut_perf))
            ticks = np.linspace(1, max_power, 5)
            ticks[1:-1] = np.round(ticks[1:-1])
            ticks = np.unique(ticks)
        elif ratio_cut_perf >= 2:
            max_power = int(np.floor(ratio_cut_perf))
            ticks = np.arange(1, max_power + 1)
        elif ratio_cut_perf >= 1 + 1e-3:
            ticks = np.array([1, ratio_cut_perf])
        else:
            ticks = np.array([1])
        tick_labels = [str(int(t)) if t == int(t) else f'{t:.2f}' for t in ticks]
    return ticks, tick_labels


def data_ticks(ratio_cut_data, is_semilogx):
    """
    Generate ticks and tick labels for data profiles.
    Matches MATLAB's dataTicks function.
    """
    if is_semilogx:
        # Data profiles are plotted at z = log2(1 + alpha), where
        # alpha = n_eval / (n + 1). Choose ticks in the original alpha
        # units and then map them to z so that labels remain literal
        # simplex-gradient budgets (0, 1, 2, 4, ...).
        ratio_cut_alpha = max(0.0, 2 ** ratio_cut_data - 1.0)
        if ratio_cut_alpha >= 1:
            max_power = int(np.floor(np.log2(ratio_cut_alpha)))
            if max_power + 1 <= 5:
                powers = np.arange(0, max_power + 1)
            else:
                powers = np.unique(np.concatenate([[0], np.floor(np.linspace(1, max_power, 4)).astype(int)]))
            alpha_ticks = np.concatenate([[0], 2.0 ** powers])
            ticks = np.log2(1.0 + alpha_ticks)
            tick_labels = [str(int(tick)) for tick in alpha_ticks]
        elif ratio_cut_alpha >= 1e-1:
            ticks = np.array([0, ratio_cut_data])
            tick_labels = ['0', f'{ratio_cut_alpha:.2f}'.rstrip('0').rstrip('.')]
        else:
            ticks = np.array([0])
            tick_labels = ['0']
    else:
        if ratio_cut_data >= 5:
            max_power = int(np.floor(ratio_cut_data))
            ticks = np.linspace(0, max_power, 5)
            ticks[1:-1] = np.round(ticks[1:-1])
            ticks = np.unique(ticks)
        elif ratio_cut_data >= 1:
            max_power = int(np.floor(ratio_cut_data))
            ticks = np.arange(0, max_power + 1)
        elif ratio_cut_data >= 1e-1:
            ticks = np.array([0, ratio_cut_data])
        else:
            ticks = np.array([0])
        tick_labels = [str(int(t)) if t == int(t) else f'{t:.2f}' for t in ticks]
    return ticks, tick_labels


def set_profile_context(profile_options):
    profile_context = {
        'font.family': 'serif',
        'font.size': 16,
        'axes.labelsize': 16,
        'axes.titlesize': 18,
        'legend.fontsize': 14,
        'xtick.labelsize': 14,
        'ytick.labelsize': 14,
        'text.usetex': True if shutil.which('latex') else False,
    }
    return profile_context


def format_float_scientific_latex(x, digits=3):
    """
    Format a floating-point number as scientific notation in LaTeX.
    
    Parameters:
    -----------
    x : float
        The number to format
    digits : int, optional
        Number of significant digits. Default is 3.
    
    Returns:
    --------
    raw : str
        Raw scientific notation
    formatted : str
        LaTeX formatted scientific notation
    """
    # Handle zero
    if x == 0:
        return "0", "0"
    
    # Format with numpy
    raw = np.format_float_scientific(x, precision=digits-1, trim='-', exp_digits=0)
    
    # More robust regex that handles various scientific notation formats:
    # - Captures optional negative sign
    # - Handles coefficient with or without decimal points
    # - Accepts both positive and negative exponents with optional sign
    match = re.compile(r'^(?P<sign>-)?(?P<coefficient>[0-9]+\.?[0-9]*)e(?P<exponent>[-+]?[0-9]+)$').match(raw)
    
    if not match:
        # Fallback method if regex matching fails
        sign = '-' if x < 0 else ''
        abs_x = abs(x)
        exponent = int(np.floor(np.log10(abs_x)))
        coefficient = abs_x / (10 ** exponent)
        
        # Limit significant digits
        coefficient = round(coefficient, digits-1)
        
        # Check if rounding caused coefficient to be >= 10
        if coefficient >= 10:
            coefficient /= 10
            exponent += 1
            
        # Format coefficient to desired precision
        coefficient_str = f"{coefficient:.{digits-1}f}".rstrip('0').rstrip('.')
        raw = f"{sign}{coefficient_str}e{exponent}"
    else:
        sign = match.group('sign') or ''
        coefficient = match.group('coefficient')
        exponent = int(match.group('exponent'))
        
        # Check if coefficient is approximately 10
        if float(coefficient) >= 9.99999:
            coefficient = '1'
            exponent += 1
            raw = f"{sign}1e{exponent}"
    
    # Handle special cases
    if exponent == 0:
        # No need for exponent part when it's 0
        return f"{sign}{coefficient}", f"{sign}{coefficient}"
    elif coefficient == '1':
        # Simplify to just 10^{exponent} when coefficient is 1
        return raw, f"{sign}10^{{{exponent}}}"
    
    # Standard scientific notation in LaTeX
    return raw, f"{sign}{coefficient} \\times 10^{{{exponent}}}"


def draw_hist(fun_histories, maxcv_histories, merit_histories, fun_init, maxcv_init, merit_init, solver_names, list_axs_summary, is_cum, ptype, problem_n, n_eval, profile_options, default_height, show_xlabel=True):
    """
    Draws the history plots of the function values, the maximum constraint violation, 
    and the merit function values.
    
    Parameters:
    -----------
    fun_histories : numpy.ndarray
        Histories of function values, shape ``(n_solvers, n_runs, n_evals)``.
    maxcv_histories : numpy.ndarray
        Histories of maximum constraint violations, same shape as
        ``fun_histories``.
    merit_histories : numpy.ndarray
        Histories of merit function values, same shape as ``fun_histories``.
    fun_init, maxcv_init, merit_init : array-like or float
        Initial values. May be either a scalar (legacy) or a per-run
        vector of length ``n_runs`` so that ``process_hist_y_axes`` can
        substitute the matching run's initial value for non-finite
        history entries (mirroring MATLAB's ``drawHist.m`` /
        ``processHistYaxes.m``).
    solver_names : list
        Names of the solvers
    list_axs_summary : list
        List of axes for plotting
    is_cum : bool
        Whether to use cumulative minimum
    ptype : str
        Problem type ('u' for unconstrained)
    problem_n : int
        Problem dimension
    n_eval : numpy.ndarray
        Number of evaluations
    profile_options : dict
        Options for plotting
    default_height : float
        Default height for sizing
    """
    
    fun_histories, fun_note = process_hist_y_axes(fun_histories, fun_init, return_note=True)
    maxcv_histories, maxcv_note = process_hist_y_axes(maxcv_histories, maxcv_init, return_note=True)
    merit_histories, merit_note = process_hist_y_axes(merit_histories, merit_init, return_note=True)
    
    # Convert default_height from inches to pixels for fontsize calculation.
    # MATLAB uses pixels directly, while matplotlib uses inches.
    dpi = matplotlib.rcParams.get('figure.dpi', 100)
    default_height_px = default_height * dpi
    profile_context = set_profile_context(profile_options)
    
    def set_hist_ylabel(ax, base_label, y_shift, formatted_shift):
        """Helper function to set ylabel with dynamic fontsize."""
        if y_shift > 0:
            label = f"{base_label} shifted above by ${formatted_shift}$"
        else:
            label = base_label
        fontsize = min(10, 1.2 * default_height_px / len(f"{base_label} shifted above by ${formatted_shift}$"))
        with matplotlib.rc_context(profile_context):
            ax.set_ylabel(label, fontsize=min(fontsize, _HISTORY_LABELSIZE), labelpad=_HISTORY_YLABEL_PAD)
    
    # Draw and label function value histories.
    y_shift_fun = compute_y_shift(fun_histories, profile_options)
    draw_fun_maxcv_merit_hist(list_axs_summary[0], fun_histories, solver_names, is_cum, problem_n, y_shift_fun, n_eval, profile_options, show_xlabel)
    _, formatted_fun_shift = format_float_scientific_latex(y_shift_fun)
    base_label_fun = "Cummin of function values" if is_cum else "Function values"
    set_hist_ylabel(list_axs_summary[0], base_label_fun, y_shift_fun, formatted_fun_shift)
    _annotate_history_display(list_axs_summary[0], fun_note)
    
    # Return early for unconstrained problems.
    if ptype == 'u':
        return
    
    # Draw and label maxcv and merit histories.
    y_shift_maxcv = compute_y_shift(maxcv_histories, profile_options)
    y_shift_merit = compute_y_shift(merit_histories, profile_options)
    
    draw_fun_maxcv_merit_hist(list_axs_summary[1], maxcv_histories, solver_names, is_cum, problem_n, y_shift_maxcv, n_eval, profile_options, show_xlabel)
    _, formatted_maxcv_shift = format_float_scientific_latex(y_shift_maxcv)
    base_label_maxcv = "Cummin of maximum constraint violations" if is_cum else "Maximum constraint violations"
    set_hist_ylabel(list_axs_summary[1], base_label_maxcv, y_shift_maxcv, formatted_maxcv_shift)
    _annotate_history_display(list_axs_summary[1], maxcv_note)
    
    draw_fun_maxcv_merit_hist(list_axs_summary[2], merit_histories, solver_names, is_cum, problem_n, y_shift_merit, n_eval, profile_options, show_xlabel)
    _, formatted_merit_shift = format_float_scientific_latex(y_shift_merit)
    base_label_merit = "Cummin of merit function values" if is_cum else "Merit function values"
    set_hist_ylabel(list_axs_summary[2], base_label_merit, y_shift_merit, formatted_merit_shift)
    _annotate_history_display(list_axs_summary[2], merit_note)


def _annotate_history_display(ax, note):
    if note:
        ax.text(0.01, 0.01, note, transform=ax.transAxes, fontsize=7,
                ha='left', va='bottom', usetex=False,
                bbox={'facecolor': 'white', 'alpha': 0.85, 'edgecolor': 'none'})


def compute_y_shift(history, profile_options):
    """
    Compute the shift of the y-axis.
    
    Parameters:
    -----------
    history : numpy.ndarray
        History data
    profile_options : dict
        Options for plotting
        
    Returns:
    --------
    y_shift : float
        The shift value for the y-axis
    """
    y_shift = 0
    if profile_options[ProfileOption.ERRORBAR_TYPE] == 'meanstd':
        y_mean = np.mean(history, axis=1)
        # Preserve Python's population standard deviation (N denominator).
        # MATLAB histories use N-1 for N>1; changing this would alter old bands.
        y_std = np.std(history, axis=1, ddof=0)
        y_lower = y_mean - y_std
        y_min = np.min(y_lower)
    else:
        y_min = np.min(history)
    
    # Strictly positive tiny values are valid on a log axis. Adding eps to
    # them would erase their relative variation, so shift only nonpositive data.
    if np.any(np.diff(history.flatten())) and y_min <= 0:
        y_shift = max(np.finfo(float).eps - y_min, np.spacing(-y_min) - y_min)
    
    return y_shift


def process_hist_y_axes(value_histories, value_inits, return_note=False):
    """
    Process the value_histories of the y-axis data.

    ``value_histories`` has shape ``(n_solvers, n_runs, n_evals)`` and
    ``value_inits`` is either a scalar (legacy) or a per-run vector of
    length ``n_runs``. The legacy 2-D input accepts only a scalar initial value.
    Only the returned display copy is clipped to +/-1e100
    before statistics/axis arithmetic. Non-finite entries are shown above the
    matching run's finite range; absent finite data and initial value, use 1.
    The optional note must be displayed whenever this transformation occurs.
    Raw histories, solver evaluations and scores are never changed here.
    """
    original = np.asarray(value_histories, dtype=float)
    processed = np.clip(original, -_HISTORY_DISPLAY_LIMIT, _HISTORY_DISPLAY_LIMIT)
    notes = []
    n_clipped = np.count_nonzero(np.isfinite(original) & (np.abs(original) > _HISTORY_DISPLAY_LIMIT))
    n_nonfinite = np.count_nonzero(~np.isfinite(original))
    if n_clipped:
        notes.append(f'Display clipped at +/-1e100: {n_clipped} entries')
    if n_nonfinite:
        notes.append(f'Nonfinite placeholders: {n_nonfinite} entries')
    # Preserve the legacy scalar/2-D helper interface as one range; actual
    # history plotting always uses a separate range for each run.
    by_run = original.ndim == 3
    n_runs = original.shape[1] if by_run else 1
    inits = np.broadcast_to(np.asarray(value_inits, dtype=float).reshape(-1), (n_runs,))
    n_empty_runs = 0
    for i_run in range(n_runs):
        source_run = original[:, i_run, :] if by_run else original
        slice_run = processed[:, i_run, :] if by_run else processed
        mask_run = ~np.isfinite(source_run)
        if not np.any(mask_run):
            continue
        finite = slice_run[~mask_run]
        if np.isfinite(inits[i_run]):
            finite = np.append(finite, np.clip(inits[i_run], -_HISTORY_DISPLAY_LIMIT, _HISTORY_DISPLAY_LIMIT))
        if finite.size:
            low, high = np.min(finite), np.max(finite)
            # Clip first: even high-low and the placeholder cannot overflow.
            # A constant finite run still needs a visibly distinct placeholder.
            gap = max(0.5 * (high - low), abs(high) * 0.05, np.finfo(float).eps)
            replacement = high + gap
        else:
            replacement = 1.0
            n_empty_runs += 1
        slice_run[mask_run] = replacement
    if n_empty_runs:
        notes.append(f'No finite reference in {n_empty_runs} run(s): shown at 1')
    note = '\n'.join(notes)
    return (processed, note) if return_note else processed


def prepare_history_plot_data(y, is_cum, y_shift, n_eval, profile_options):
    """Prepare indices, mean, lower and upper for each solver without axes.

    ``y`` is the processed display copy shaped (solver, run, evaluation).
    Returned indices are 1-based actual evaluation numbers, not array offsets.
    All arithmetic below preserves the historic renderer's operation order.
    """
    
    n_solvers = y.shape[0]
    n_runs = y.shape[1]
    # Keep (solver, evaluation) even for one solver or one evaluation.
    y_mean = np.mean(y, axis=1)

    if profile_options[ProfileOption.ERRORBAR_TYPE] == 'minmax':
        y_lower = np.min(y, axis=1)
        y_upper = np.max(y, axis=1)
    else:
        # Match compute_y_shift and the established Python population bands;
        # MATLAB retains its sample normalization (see history_display docs).
        y_std = np.std(y, axis=1, ddof=0)
        y_lower = y_mean - y_std
        y_upper = y_mean + y_std

    # Compute the band before translating, exactly as compute_y_shift does.
    # Recomputing variance after a large translation can reintroduce negative
    # lower bands through cancellation even when the selected shift was safe.
    y_mean = y_mean + y_shift
    y_lower = y_lower + y_shift
    y_upper = y_upper + y_shift

    if is_cum:
        y_mean = np.minimum.accumulate(y_mean, axis=1)
        y_lower = np.minimum.accumulate(y_lower, axis=1)
        y_upper = np.minimum.accumulate(y_upper, axis=1)

    """
    Block Aggregation for Large Evaluation Histories
    
    When the number of evaluation points is very large, directly plotting all
    y-values becomes computationally expensive (both time and memory). To
    avoid this, we apply block aggregation, which means that we group
    neighboring points into small segments ("blocks") and keep only one
    representative value from each block. This greatly reduces the total
    number of points we need to draw, while still keeping the main shape of
    the curve.
    
    There are three ways (modes) to get the value for each block:
    
    1. 'min'
        - Find the smallest y_mean value inside the block.
        - Its position (x index) is used as the x coordinate (abs_idx).
        - The corresponding y_lower and y_upper at that same position are used.
        - This highlights the best (lowest) performance in that part of the curve.

    2. 'mean'
        - Use the middle x position of the block (the median of start and end indices).
        - Compute the average of all y_mean, y_lower, and y_upper values in the block.
        - This produces a smooth averaged curve that represents the general trend.

    3. 'max'
        - Find the largest y_mean value inside the block.
        - Its position (x index) is used as the x coordinate (abs_idx).
        - The corresponding y_lower and y_upper at that same position are used.
        - This shows the worst (highest) value in that part of the curve.
    
    The first and last points are always kept. Only the points between them
    are affected by the block aggregation.
    
    After doing this, we only keep around `n_blocks` (1000) points instead of
    all evaluations. This keeps the figure clear and makes plotting faster.
    """
    max_eval = y.shape[2]
    n_blocks = 1000

    # We initialize the cell array to store the x indices and y values to be plotted.
    x_indices = [[] for _ in range(n_solvers)]
    y_values_m = [[] for _ in range(n_solvers)]
    y_values_l = [[] for _ in range(n_solvers)]
    y_values_u = [[] for _ in range(n_solvers)]

    # We will build blocks only for 2:(max_eval-1) so the first and last points
    # are excluded from block aggregation and will be appended unconditionally
    # at the end.
    if max_eval > 2:
        inner_len = max_eval - 2
        n_blocks_eff = min(n_blocks, inner_len)
        q = inner_len // n_blocks_eff
        r = inner_len % n_blocks_eff
        blocks = np.full(n_blocks_eff, q, dtype=int)
        if r > 0:
            # Distribute the remaining r elements evenly across the blocks.
            # Choose approximately r equally spaced block indices using linspace,
            # and add one extra element to each of these blocks so that
            # the total number of elements sums back to inner_len.
            idxs = np.round(np.linspace(0, n_blocks_eff - 1, r)).astype(int)
            blocks[idxs] += 1

        # We aggregate over inner blocks (indices 1, ..., max_eval-2).
        for i_block in range(n_blocks_eff):
            idx_start = int(sum(blocks[:i_block])) + 1  # Add 1 to account skip the first point.
            idx_end = int(idx_start + blocks[i_block])
            idx = np.arange(idx_start, idx_end)
            for i_solver in range(n_solvers):
                i_eval = int(np.max(n_eval[i_solver, :]))
                if profile_options[ProfileOption.HIST_AGGREGATION] == 'min':
                    y_value_m = np.nanmin(y_mean[i_solver, idx])
                    rel_idx = int(np.nanargmin(y_mean[i_solver, idx]))
                    # NumPy array offsets are zero based, but plotted
                    # evaluation numbers and n_eval are one based.
                    abs_idx = idx[rel_idx] + 1
                    y_value_l = y_lower[i_solver, abs_idx - 1]
                    y_value_u = y_upper[i_solver, abs_idx - 1]
                elif profile_options[ProfileOption.HIST_AGGREGATION] == 'mean':
                    abs_idx = (idx_start + idx_end + 1) // 2
                    y_value_m = np.nanmean(y_mean[i_solver, idx])
                    y_value_l = np.nanmean(y_lower[i_solver, idx])
                    y_value_u = np.nanmean(y_upper[i_solver, idx])
                elif profile_options[ProfileOption.HIST_AGGREGATION] == 'max':
                    y_value_m = np.nanmax(y_mean[i_solver, idx])
                    rel_idx = int(np.nanargmax(y_mean[i_solver, idx]))
                    abs_idx = idx[rel_idx] + 1
                    y_value_l = y_lower[i_solver, abs_idx - 1]
                    y_value_u = y_upper[i_solver, abs_idx - 1]
                if abs_idx > i_eval:
                    continue
                x_indices[i_solver].append(abs_idx)
                y_values_m[i_solver].append(y_value_m)
                y_values_l[i_solver].append(y_value_l)
                y_values_u[i_solver].append(y_value_u)

    # We add the first and last indices unconditionally.
    for i_solver in range(n_solvers):
        i_eval = int(np.max(n_eval[i_solver, :]))
        if i_eval > 0:
            # Add the first point.
            x_indices[i_solver] = [1] + x_indices[i_solver]
            y_values_m[i_solver] = [y_mean[i_solver, 0]] + y_values_m[i_solver]
            y_values_l[i_solver] = [y_lower[i_solver, 0]] + y_values_l[i_solver]
            y_values_u[i_solver] = [y_upper[i_solver, 0]] + y_values_u[i_solver]

            # Add the last point.                
            if i_eval != 1:
                x_indices[i_solver].append(i_eval)
                y_values_m[i_solver].append(y_mean[i_solver, i_eval - 1])
                y_values_l[i_solver].append(y_lower[i_solver, i_eval - 1])
                y_values_u[i_solver].append(y_upper[i_solver, i_eval - 1])

    # This is the single numerical preparation used by both the renderer and
    # EvalReport. Keep aggregation/translation/cummin/block order unchanged;
    # these arrays are a display representation, never new scoring inputs.
    return x_indices, y_values_m, y_values_l, y_values_u


def prepare_history_panels(histories, initials, ptype, problem_n, n_eval, profile_options):
    """Numeric history panels for the same observed inputs given to draw_hist.

    No solver/merit callback runs here. The renderer can supply its own already
    computed merit values, which matters when a user callback is stateful and
    the earlier scoring-time merit observations are not identical.
    """
    panels = []
    for channel in ('objective', 'constraint', 'merit'):
        if ptype == 'u' and channel != 'objective':
            continue
        values, initial = histories.get(channel), initials.get(channel)
        if not isinstance(values, np.ndarray) or initial is None or not values.size:
            continue
        shown, note = process_hist_y_axes(values, initial, return_note=True)
        shift = compute_y_shift(shown, profile_options)
        for mode in ('raw', 'cummin'):
            indices, means, lower, upper = prepare_history_plot_data(
                shown, mode == 'cummin', shift, n_eval, profile_options)
            series = [{'solver_index': solver + 1,
                       'evaluation_indices': indices[solver],
                       'x': np.asarray(indices[solver]) / (problem_n + 1),
                       'mean': means[solver], 'lower': lower[solver], 'upper': upper[solver],
                       'geometry': 'point' if len(indices[solver]) == 1 else 'line',
                       'band_visible': shown.shape[1] > 1 and len(indices[solver]) > 1}
                      for solver in range(shown.shape[0])]
            # Same panel vocabulary as MATLAB prepareEvalReportHistory.m.
            # Block aggregation is keyed on the PADDED history length
            # (max_eval = ceil(max_eval_factor*dimension)): above 1002 the
            # renderer keeps about 1000 interior points, so an array position
            # in these series is not an evaluation number; readers must use
            # evaluation_indices.
            panels.append({'kind': 'history', 'channel': channel, 'mode': mode,
                           'series': series, 'x_transform': 'evaluation_index/(dimension+1)',
                           'y_scale': 'log' if any(np.any(m) and np.any(np.diff(m)) for m in means) else 'linear',
                           'y_shift': shift, 'display_limit': _HISTORY_DISPLAY_LIMIT,
                           'display_note': note or None,
                           'nonfinite_policy': 'per_run_above_finite_range_with_initial_fallback',
                           'aggregation': profile_options[ProfileOption.HIST_AGGREGATION],
                           'aggregation_trigger': 'padded_history_length_above_1002',
                           'padded_length': int(shown.shape[2]),
                           'errorbar_type': profile_options[ProfileOption.ERRORBAR_TYPE],
                           'n_runs': int(shown.shape[1]), 'std_ddof': 0})
    return panels


def draw_fun_maxcv_merit_hist(ax, y, solver_names, is_cum, problem_n, y_shift, n_eval, profile_options, show_xlabel=True):
    """Render the shared, pure history preparation without changing its data."""
    profile_context = set_profile_context(profile_options)
    line_colors = profile_options[ProfileOption.LINE_COLORS]
    line_styles = profile_options[ProfileOption.LINE_STYLES]
    line_widths = profile_options[ProfileOption.LINE_WIDTHS]
    n_solvers, n_runs = y.shape[:2]
    x_indices, y_values_m, y_values_l, y_values_u = prepare_history_plot_data(
        y, is_cum, y_shift, n_eval, profile_options)
    xl_lim = 1 / (problem_n + 1)
    xr_lim = 1 / (problem_n + 1)
    is_log_scale = False
    handles = []

    with matplotlib.rc_context(profile_context):
        for i_solver in range(n_solvers):
            # Truncate the histories according to the function evaluations of each solver.
            i_x = x_indices[i_solver]
            i_y_mean = y_values_m[i_solver]
            i_y_lower = y_values_l[i_solver]
            i_y_upper = y_values_u[i_solver]
            i_eval = len(i_x)
            if i_eval == 0:
                continue
            x = np.array(i_x) / (problem_n + 1)
            xr_lim = max(xr_lim, x[-1])

            color = line_colors[i_solver % len(line_colors)]
            line_style = line_styles[i_solver % len(line_styles)]
            line_width = line_widths[i_solver % len(line_widths)]

            if i_eval == 1:
                handles.extend(ax.plot(x, i_y_mean, color=color, marker='o', label=solver_names[i_solver]))
            elif i_eval > 1:
                handles.extend(ax.plot(x, i_y_mean, color=color, linestyle=line_style, linewidth=line_width, label=solver_names[i_solver]))
                if n_runs > 1:
                    ax.fill_between(x, i_y_lower, i_y_upper, color=color, alpha=0.2)
            if np.any(i_y_mean) and np.any(np.diff(i_y_mean)):
                is_log_scale = True

    # When the function values are not all zero and there is at least some change in the function values, use log scale for the y-axis.
    if is_log_scale:
        ax.set_yscale('log')
        shown = np.concatenate([np.asarray(values) for group in
                                (y_values_m, y_values_l, y_values_u)
                                for values in group if len(values)])
        positive = shown[shown > 0]
        if positive.size:
            log_low, log_high = np.log10(np.min(positive)), np.log10(np.max(positive))
            margin = max(0.05 * (log_high - log_low), 0.05)
            # Explicit finite limits prevent autoscale's silent [1, 10]
            # fallback. Do not floor valid tiny values in the displayed data.
            smallest = np.nextafter(0.0, 1.0)
            low = max(smallest, 10.0 ** max(np.log10(smallest), log_low - margin))
            high = 10.0 ** min(102.0, log_high + margin)
            ax.set_ylim(low, high)
            # Automatic LogLocator may request decades outside float64 even
            # when all data/limits are finite. Restrict only extreme axes;
            # ordinary histories keep Matplotlib's usual tick placement.
            if log_low < -250 or log_high - log_low > 200:
                first, last = int(np.ceil(np.log10(low))), int(np.floor(np.log10(high)))
                step = max(1, int(np.ceil((last - first) / 8)))
                ticks = np.power(10.0, np.arange(first, last + 1, step, dtype=float))
                ax.yaxis.set_major_locator(FixedLocator(ticks))
                ax.yaxis.set_minor_locator(NullLocator())
    
    # Show ticks on both sides (left and right) and set direction to 'in' (matching MATLAB).
    ax.yaxis.set_ticks_position('both')
    ax.tick_params(axis='y', which='both', direction='in', right=True, labelleft=True, labelright=False, pad=4)
    ax.tick_params(axis='x', which='both', direction='in')
    # Add box around the plot (matching MATLAB's box(ax, 'on')).
    for spine in ax.spines.values():
        spine.set_visible(True)
    if abs(xl_lim - xr_lim) < np.finfo(float).eps:
        xr_lim = xl_lim + np.finfo(float).eps
    
    ax.tick_params(axis='both', which='major', labelsize=_HISTORY_TICK_LABELSIZE)
    ax.set_xlim([xl_lim, xr_lim])
    if show_xlabel:
        ax.set_xlabel(profile_options[ProfileOption.XLABEL_DATA_PROFILE], fontsize=_HISTORY_LABELSIZE, labelpad=12)
    else:
        ax.set_xlabel('')
    legend = _place_solver_legend(ax, n_solvers, default_loc='upper right', handles=handles)
    if legend is not None:
        for text in legend.get_texts():
            text.set_fontsize(_HISTORY_LEGEND_FONTSIZE)


def _compact_legend_fontsize(n_solvers):
    if n_solvers <= _COMPACT_LEGEND_SOLVER_THRESHOLD:
        return None
    if n_solvers <= 20:
        return 10
    if n_solvers <= 30:
        return 9
    return 8


def _legend_columns(n_solvers):
    if n_solvers <= _COMPACT_LEGEND_SOLVER_THRESHOLD:
        return 1
    return max(1, n_solvers // _COMPACT_LEGEND_MIN_ROWS_PER_COLUMN)


def summary_legend_extra_width(n_solvers, default_width, solver_names=None):
    if n_solvers <= _COMPACT_LEGEND_SOLVER_THRESHOLD:
        return 0.0
    if solver_names:
        max_label_length = max(len(str(name).replace(r'\_', '_')) for name in solver_names)
    else:
        max_label_length = 12
    per_column_fraction = min(0.45, max(0.22, 0.08 + 0.015 * max_label_length))
    return default_width * (0.05 + per_column_fraction * _legend_columns(n_solvers))


def _normalize_axis_values(values, axis_limits, axis_scale):
    values = np.asarray(values, dtype=float)
    limits = np.asarray(axis_limits, dtype=float)
    finite_mask = np.isfinite(values)
    if axis_scale == 'log':
        finite_mask &= values > 0.0
        finite_mask &= np.all(limits > 0.0)
        if not np.any(finite_mask):
            return np.array([])
        values = np.log10(values[finite_mask])
        limits = np.log10(limits)
    else:
        values = values[finite_mask]
    span = limits[1] - limits[0]
    if not np.isfinite(span) or abs(span) < np.finfo(float).eps:
        return np.array([])
    normalized = (values - limits[0]) / span
    normalized = normalized[np.isfinite(normalized)]
    return np.clip(normalized, 0.0, 1.0)


def _choose_legend_location(ax, default_loc, handles=None):
    if handles is None:
        handles, _ = ax.get_legend_handles_labels()
    if not handles:
        return default_loc

    xlim = ax.get_xlim()
    ylim = ax.get_ylim()
    x_values = []
    y_values = []
    for handle in handles:
        if not isinstance(handle, Line2D):
            continue
        x_norm = _normalize_axis_values(handle.get_xdata(), xlim, ax.get_xscale())
        y_norm = _normalize_axis_values(handle.get_ydata(), ylim, ax.get_yscale())
        n = min(x_norm.size, y_norm.size)
        if n == 0:
            continue
        x_values.append(x_norm[:n])
        y_values.append(y_norm[:n])
    if not x_values:
        return default_loc

    x_all = np.concatenate(x_values)
    y_all = np.concatenate(y_values)
    if x_all.size == 0:
        return default_loc

    targets = {
        'upper right': (0.85, 0.85),
        'upper left': (0.15, 0.85),
        'lower right': (0.85, 0.15),
        'lower left': (0.15, 0.15),
    }
    scores = {}
    for loc, (target_x, target_y) in targets.items():
        distance2 = (x_all - target_x) ** 2 + (y_all - target_y) ** 2
        scores[loc] = float(np.sum(np.exp(-distance2 / 0.06)))
    return min(targets, key=lambda loc: (scores[loc], loc != default_loc))


def _place_solver_legend(ax, n_solvers, default_loc='lower right', handles=None):
    if handles is None:
        handles, _ = ax.get_legend_handles_labels()
    if not handles:
        return None
    labels = [handle.get_label() for handle in handles]
    # Solver identifiers beginning with '_' are literal names, not requests
    # to suppress their curves. Older Matplotlib also filters explicit labels:
    # construct with neutral labels, then restore the exact visible names.
    legend_args = (handles, [str(index) for index in range(len(handles))])
    if n_solvers > _COMPACT_LEGEND_SOLVER_THRESHOLD:
        legend = ax.legend(
            *legend_args,
            loc='center left',
            bbox_to_anchor=(1.02, 0.5),
            ncol=_legend_columns(n_solvers),
            borderaxespad=0.0,
            fontsize=_compact_legend_fontsize(n_solvers),
            handlelength=2.0,
            columnspacing=0.9,
            labelspacing=0.35,
        )
    else:
        legend = ax.legend(*legend_args, loc=_choose_legend_location(ax, default_loc, handles))
    if legend is not None:
        for text, label in zip(legend.get_texts(), labels):
            text.set_text(label)
        legend.set_in_layout(True)
    return legend
