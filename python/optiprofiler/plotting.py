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

import numpy as np
from cycler import cycler
from matplotlib import pyplot as plt
from matplotlib.colors import is_color_like
from matplotlib.lines import Line2D
from matplotlib.backends import backend_pdf
from matplotlib.ticker import MaxNLocator, FuncFormatter

from .utils import ProfileOption



def draw_profiles(work, problem_dimensions, solver_names, tolerance_latex, i_tol, ax_summary_perf, ax_summary_data, ax_summary_log_ratio, is_summary, is_perf, is_data, is_log_ratio, profile_options, curves):
    solver_names = [name.replace('_', r'\_') for name in solver_names]
    n_solvers = work.shape[1]

    # Create the individual figures.
    fig_perf, ax_perf = plt.subplots()
    fig_data, ax_data = plt.subplots()
    if n_solvers > 2:
        fig_log_ratio, ax_log_ratio = None, None
    else:
        fig_log_ratio, ax_log_ratio = plt.subplots()

    # Calculate the performance and data profiles.
    x_perf, y_perf, ratio_max_perf, x_data, y_data, ratio_max_data, curves = _get_extended_performances_data_profile_axes(work, problem_dimensions, profile_options, curves)
    if n_solvers == 2:
        x_log_ratio, y_log_ratio, ratio_max_log_ratio, n_solvers_fail, curves = _get_log_ratio_profile_axes(work, curves)
    
    if profile_options[ProfileOption.SCORE_ONLY]:
        return fig_perf, fig_data, fig_log_ratio, curves

    

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
        







    # if ax_summary_perf_hist is not None:
    #     if i_tolerance == 0:
    #         _draw_performance_data_profiles(ax_summary_perf_hist[i_tolerance], x_perf_hist, y_perf_hist, labels)
    #     else:
    #         _draw_performance_data_profiles(ax_summary_perf_hist[i_tolerance], x_perf_hist, y_perf_hist)
    # if ax_summary_perf_out is not None:
    #     _draw_performance_data_profiles(ax_summary_perf_out[i_tolerance], x_perf_out, y_perf_out)
    # perf_formatter = FuncFormatter(_perf_formatter)
    # ax_perf_hist.xaxis.set_major_formatter(perf_formatter)
    # ax_perf_out.xaxis.set_major_formatter(perf_formatter)
    # if ax_summary_perf_hist is not None:
    #     ax_summary_perf_hist[i_tolerance].xaxis.set_major_formatter(perf_formatter)
    # if ax_summary_perf_out is not None:
    #     ax_summary_perf_out[i_tolerance].xaxis.set_major_formatter(perf_formatter)
    # with warnings.catch_warnings():
    #     warnings.filterwarnings('ignore', category=UserWarning)
    #     ax_perf_hist.set_xlim(0.0, 1.05 * ratio_max_perf_hist)
    #     ax_perf_out.set_xlim(0.0, 1.05 * ratio_max_perf_out)
    #     if ax_summary_perf_hist is not None:
    #         ax_summary_perf_hist[i_tolerance].set_xlim(0.0, 1.05 * ratio_max_perf_hist)
    #     if ax_summary_perf_out is not None:
    #         ax_summary_perf_out[i_tolerance].set_xlim(0.0, 1.05 * ratio_max_perf_out)
    # ax_perf_hist.set_xlabel('Performance ratio')
    # ax_perf_out.set_xlabel('Performance ratio')
    # if ax_summary_perf_hist is not None:
    #     ax_summary_perf_hist[i_tolerance].set_xlabel('Performance ratio')
    # if ax_summary_perf_out is not None:
    #     ax_summary_perf_out[i_tolerance].set_xlabel('Performance ratio')
    # ax_perf_hist.set_ylabel(f'Performance profiles {tolerance_label}')
    # ax_perf_out.set_ylabel(f'Performance profiles {tolerance_label}')
    # if ax_summary_perf_hist is not None:
    #     ax_summary_perf_hist[i_tolerance].set_ylabel(f'Performance profiles {tolerance_label}')
    # if ax_summary_perf_out is not None:
    #     ax_summary_perf_out[i_tolerance].set_ylabel(f'Performance profiles {tolerance_label}')

    # # Draw the data profiles.
    # _draw_performance_data_profiles(ax_data_hist, x_data_hist, y_data_hist, labels)
    # _draw_performance_data_profiles(ax_data_out, x_data_out, y_data_out, labels)
    # if ax_summary_data_hist is not None:
    #     _draw_performance_data_profiles(ax_summary_data_hist[i_tolerance], x_data_hist, y_data_hist)
    # if ax_summary_data_out is not None:
    #     _draw_performance_data_profiles(ax_summary_data_out[i_tolerance], x_data_out, y_data_out)
    # data_formatter = FuncFormatter(_data_formatter)
    # ax_data_hist.xaxis.set_major_formatter(data_formatter)
    # ax_data_out.xaxis.set_major_formatter(data_formatter)
    # if ax_summary_data_hist is not None:
    #     ax_summary_data_hist[i_tolerance].xaxis.set_major_formatter(data_formatter)
    # if ax_summary_data_out is not None:
    #     ax_summary_data_out[i_tolerance].xaxis.set_major_formatter(data_formatter)
    # ax_data_hist.set_xlim(0.0, 1.05 * ratio_max_data_hist)
    # ax_data_out.set_xlim(0.0, 1.05 * ratio_max_data_out)
    # if ax_summary_data_hist is not None:
    #     ax_summary_data_hist[i_tolerance].set_xlim(0.0, 1.05 * ratio_max_data_hist)
    # if ax_summary_data_out is not None:
    #     ax_summary_data_out[i_tolerance].set_xlim(0.0, 1.05 * ratio_max_data_out)
    # ax_data_hist.set_xlabel('Number of simplex gradients')
    # ax_data_out.set_xlabel('Number of simplex gradients')
    # if ax_summary_data_hist is not None:
    #     ax_summary_data_hist[i_tolerance].set_xlabel('Number of simplex gradients')
    # if ax_summary_data_out is not None:
    #     ax_summary_data_out[i_tolerance].set_xlabel('Number of simplex gradients')
    # ax_data_hist.set_ylabel(f'Data profiles {tolerance_label}')
    # ax_data_out.set_ylabel(f'Data profiles {tolerance_label}')
    # if ax_summary_data_hist is not None:
    #     ax_summary_data_hist[i_tolerance].set_ylabel(f'Data profiles {tolerance_label}')
    # if ax_summary_data_out is not None:
    #     ax_summary_data_out[i_tolerance].set_ylabel(f'Data profiles {tolerance_label}')

    # # Draw the log-ratio profiles.
    # if n_solvers <= 2:
    #     _draw_log_ratio_profiles(ax_log_ratio_hist, np.copy(work_hist), labels)
    #     _draw_log_ratio_profiles(ax_log_ratio_out, np.copy(work_out), labels)
    #     if ax_summary_log_ratio_hist is not None:
    #         _draw_log_ratio_profiles(ax_summary_log_ratio_hist[i_tolerance], np.copy(work_hist), labels)
    #     if ax_summary_log_ratio_out is not None:
    #         _draw_log_ratio_profiles(ax_summary_log_ratio_out[i_tolerance], np.copy(work_out), labels)
    #     ax_log_ratio_hist.set_xlabel('Problem')
    #     ax_log_ratio_out.set_xlabel('Problem')
    #     if ax_summary_log_ratio_hist is not None:
    #         ax_summary_log_ratio_hist[i_tolerance].set_xlabel('Problem')
    #     if ax_summary_log_ratio_out is not None:
    #         ax_summary_log_ratio_out[i_tolerance].set_xlabel('Problem')
    #     ax_log_ratio_hist.set_ylabel(f'Log-ratio profiles {tolerance_label}')
    #     ax_log_ratio_out.set_ylabel(f'Log-ratio profiles {tolerance_label}')
    #     if ax_summary_log_ratio_hist is not None:
    #         ax_summary_log_ratio_hist[i_tolerance].set_ylabel(f'Log-ratio profiles {tolerance_label}')
    #     if ax_summary_log_ratio_out is not None:
    #         ax_summary_log_ratio_out[i_tolerance].set_ylabel(f'Log-ratio profiles {tolerance_label}')
    return fig_perf, fig_data, fig_log_ratio, curves


def _draw_perf_detail(ax_perf, x_perf, y_perf, ratio_max_perf, solver_names, profile_options, tolerance_latex):
    _draw_performance_data_profiles(ax_perf, x_perf, y_perf, solver_names, profile_options)
    # Set x-axis limits.
    if profile_options[ProfileOption.SEMILOGX]:
        ax_perf.set_xlim(0.0, 1.1 * ratio_max_perf)
    else:
        ax_perf.set_xlim(1.0, 1.1 * ratio_max_perf)
    # Modify x-axis ticks labels of the performance profiles.




def _draw_data_detail(ax_data, x_data, y_data, ratio_max_data, solver_names, profile_options, tolerance_latex):
    _draw_performance_data_profiles(ax_data, x_data, y_data, solver_names, profile_options)
    # Set x-axis limits.
    ax_data.set_xlim(0.0, 1.1 * ratio_max_data)
    # Modify x-axis ticks labels of the data profiles.



def _draw_log_ratio_detail(ax_log_ratio, x_log_ratio, y_log_ratio, ratio_max_log_ratio, n_solvers_fail, solver_names, profile_options, tolerance_latex):
    _draw_log_ratio_profiles(ax_log_ratio, x_log_ratio, y_log_ratio, ratio_max_log_ratio, n_solvers_fail, solver_names, profile_options)
    # Set x-axis limits.







def _draw_performance_data_profiles(ax, x, y, solver_names, profile_options):
    profile_context = set_profile_context(profile_options)
    line_colors = profile_options[ProfileOption.LINE_COLORS]
    line_styles = profile_options[ProfileOption.LINE_STYLES]
    line_widths = profile_options[ProfileOption.LINE_WIDTHS]
    n_solvers = x.shape[1]
    n_runs = y.shape[2]
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

    with plt.rc_context(profile_context):
        for i_solver in range(n_solvers):
            x_stairs = np.repeat(x[:, i_solver], 2)[1:]
            y_mean_stairs = np.repeat(y_mean[:, i_solver], 2)[:-1]
            y_lower_stairs = np.repeat(y_lower[:, i_solver], 2)[:-1]
            y_upper_stairs = np.repeat(y_upper[:, i_solver], 2)[:-1]

            color = line_colors[i_solver % len(line_colors)]
            line_style = line_styles[i_solver % len(line_styles)]
            line_width = line_widths[i_solver % len(line_widths)]

            ax.plot(x_stairs, y_mean_stairs, label=solver_names[i_solver], color=color, linestyle=line_style, linewidth=line_width)
            if n_runs > 1:
                ax.fill_between(x_stairs, y_lower_stairs, y_upper_stairs, color=color, alpha=0.2)
        ax.xaxis.set_major_locator(MaxNLocator(5, integer=True))
        ax.yaxis.set_ticks_position('both')
        ax.yaxis.set_major_locator(MaxNLocator(5, prune='lower'))
        ax.yaxis.set_minor_locator(MaxNLocator(10))
        ax.set_ylim(0.0, 1.0)
        ax.tick_params(which='both', direction='in')
        ax.legend(loc='lower right')


def _draw_log_ratio_profiles(ax, x, y, ratio_max, n_solvers_equal, solver_names, profile_options):
    profile_context = set_profile_context(profile_options)
    n_problems = x.shape[0]
    n_below = np.sum(y < 0) - n_solvers_equal
    n_above = np.sum(y > 0) - n_solvers_equal

    with plt.rc_context(profile_context):
        bar_colors = profile_options[ProfileOption.BAR_COLORS][:2]
        if n_solvers_equal > 0:
            ax.bar(x[:n_solvers_equal], y[:n_solvers_equal], color=bar_colors[0], alpha=0.5)
            ax.bar(x[-n_solvers_equal:], y[-n_solvers_equal:], color=bar_colors[1], alpha=0.5)
        if n_below > 0:
            ax.bar(x[n_solvers_equal:n_solvers_equal + n_below], y[n_solvers_equal:n_solvers_equal + n_below], color=bar_colors[0])
        if n_above > 0:
            ax.bar(x[-(n_solvers_equal + n_above):-n_solvers_equal], y[-(n_solvers_equal + n_above):-n_solvers_equal], color=bar_colors[1])

        ax.text((n_problems + 1) / 2, -ratio_max, solver_names[0], horizontalalignment='center', verticalalignment='bottom')
        ax.text((n_problems + 1) / 2, ratio_max, solver_names[1], horizontalalignment='center', verticalalignment='top')
        with warnings.catch_warnings():
            warnings.filterwarnings('ignore', category=UserWarning)
            ax.set_xlim(0.5, n_problems + 0.5)
        ax.set_ylim(-1.1 * ratio_max, 1.1 * ratio_max)
        ax.tick_params(which='both', direction='in')



def _get_extended_performances_data_profile_axes(work, problem_dimensions, profile_options, curves):
    n_problems, n_solvers, n_runs = work.shape

    denominator_perf = lambda i_problem, i_run: np.nanmin(work[i_problem, :, i_run], initial=np.inf)
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

    denominator_data = lambda i_problem, i_run: problem_dimensions[i_problem] + 1
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
    log_ratio_finite = np.isfinite(work_flat[:, 0]) & np.isfinite(work_flat[:, 1])
    y_log_ratio[log_ratio_finite] = np.log2(work_flat[log_ratio_finite, 0]) - np.log2(work_flat[log_ratio_finite, 1])

    ratio_max = np.nanmax(np.abs(y_log_ratio[log_ratio_finite]), initial=np.finfo(float).eps)
    if ratio_max == 0 or np.isnan(ratio_max):
        ratio_max = np.finfo(float).eps
    
    mask1 = np.isnan(work_flat[:, 0]) & np.isfinite(work_flat[:, 1])
    y_log_ratio[mask1] = 1.1 * ratio_max
    mask2 = np.isfinite(work_flat[:, 0]) & np.isnan(work_flat[:, 1])
    y_log_ratio[mask2] = -1.1 * ratio_max

    mask_both_fail = np.isnan(work_flat[:, 0]) & np.isnan(work_flat[:, 1])
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


def draw_hist(fun_histories, maxcv_histories, merit_histories, fun_init, maxcv_init, merit_init, solver_names, list_axs_summary, is_cum, ptype, problem_n, n_eval, profile_options, default_height):
    """
    Draws the history plots of the function values, the maximum constraint violation, 
    and the merit function values.
    
    Parameters:
    -----------
    fun_histories : numpy.ndarray
        Histories of function values
    maxcv_histories : numpy.ndarray
        Histories of maximum constraint violations
    merit_histories : numpy.ndarray
        Histories of merit function values
    fun_init, maxcv_init, merit_init : float
        Initial values
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
    
    fun_histories = process_hist_y_axes(fun_histories, fun_init)
    maxcv_histories = process_hist_y_axes(maxcv_histories, maxcv_init)
    merit_histories = process_hist_y_axes(merit_histories, merit_init)
    
    # Define the shift of the y-axis. Shift the y-axis if there is value that is too close to zero.
    y_shift_fun = compute_y_shift(fun_histories, profile_options)
    
    # First, draw the histories of function values.
    draw_fun_maxcv_merit_hist(list_axs_summary[0], fun_histories, solver_names, is_cum, problem_n, y_shift_fun, n_eval, profile_options)
    _, formatted_fun_shift = format_float_scientific_latex(y_shift_fun)
    
    if is_cum:
        if y_shift_fun > 0:
            y_label = f"Cummin of function values shifted above by ${formatted_fun_shift}$"
            list_axs_summary[0].set_ylabel(y_label)
        else:
            list_axs_summary[0].set_ylabel("Cummin of function values")
    else:
        if y_shift_fun > 0:
            y_label = f"Function values shifted above by ${formatted_fun_shift}$"
            list_axs_summary[0].set_ylabel(y_label)
        else:
            list_axs_summary[0].set_ylabel("Function values")
    
    # If the problem is unconstrained, do not draw the histories of maximum constraint violations and merit function values.
    if ptype == 'u':
        return
    
    # Do the same for the maximum constraint violations and the merit function values.
    y_shift_maxcv = compute_y_shift(maxcv_histories, profile_options)
    y_shift_merit = compute_y_shift(merit_histories, profile_options)
    
    # Second, draw the histories of maximum constraint violations and merit function values.
    draw_fun_maxcv_merit_hist(list_axs_summary[1], maxcv_histories, solver_names, is_cum, problem_n, y_shift_maxcv, n_eval, profile_options)
    _, formatted_maxcv_shift = format_float_scientific_latex(y_shift_maxcv)
    
    draw_fun_maxcv_merit_hist(list_axs_summary[2], merit_histories, solver_names, is_cum, problem_n, y_shift_merit, n_eval, profile_options)
    _, formatted_merit_shift = format_float_scientific_latex(y_shift_merit)
    
    if is_cum:
        if y_shift_maxcv > 0:
            y_label = f"Cummin of maximum constraint violations shifted above by ${formatted_maxcv_shift}$"
            list_axs_summary[1].set_ylabel(y_label)
        else:
            list_axs_summary[1].set_ylabel("Cummin of maximum constraint violations")
        
        if y_shift_merit > 0:
            y_label = f"Cummin of merit function values shifted above by ${formatted_merit_shift}$"
            list_axs_summary[2].set_ylabel(y_label)
        else:
            list_axs_summary[2].set_ylabel("Cummin of merit function values")
    else:
        if y_shift_maxcv > 0:
            y_label = f"Maximum constraint violations shifted above by ${formatted_maxcv_shift}$"
            list_axs_summary[1].set_ylabel(y_label)
        else:
            list_axs_summary[1].set_ylabel("Maximum constraint violations")
        
        if y_shift_merit > 0:
            y_label = f"Merit function values shifted above by ${formatted_merit_shift}$"
            list_axs_summary[2].set_ylabel(y_label)
        else:
            list_axs_summary[2].set_ylabel("Merit function values")


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
        y_std = np.std(history, axis=1, ddof=0)  # ddof=0 for MATLAB compatibility
        y_lower = y_mean - y_std
        y_min = np.min(y_lower)
    else:
        y_min = np.min(history)
    
    # Shift the y-axis if there is value that is smaller than eps and the values are not all the same.
    if np.any(np.diff(history.flatten())) and y_min < np.finfo(float).eps:
        y_shift = max(np.finfo(float).eps - y_min, np.spacing(-y_min) - y_min)
    
    return y_shift


def process_hist_y_axes(value_histories, value_init):
    """
    Process the value_histories of the y-axis data.
    Replace non-finite values with value_init, then set those positions to value_min + 1.5 * (value_max - value_min).
    """
    value_histories = np.array(value_histories, copy=True)
    mask_hist_nan_inf = ~np.isfinite(value_histories)
    value_histories[mask_hist_nan_inf] = value_init
    value_max = np.max(value_histories)
    value_min = np.min(value_histories)
    value_histories_processed = value_histories.copy()
    value_histories_processed[mask_hist_nan_inf] = value_min + 1.5 * (value_max - value_min)
    return value_histories_processed


def draw_fun_maxcv_merit_hist(ax, y, solver_names, is_cum, problem_n, y_shift, n_eval, profile_options):
    """
    Draws figures of histories of function values, maximum constraint violation, or merit function values.
    
    Parameters:
    -----------
    ax : matplotlib.axes.Axes
        The axes to draw on
    y : numpy.ndarray
        The history data with shape (n_solvers, n_runs, n_evals)
    solver_names : list
        Names of solvers
    is_cum : bool
        Whether to use cumulative minimum
    problem_n : int
        Problem dimension
    y_shift : float
        Shift for y-axis values
    n_eval : numpy.ndarray
        Number of evaluations for each solver and run
    profile_options : dict
        Options for plotting
    """
    
    profile_context = set_profile_context(profile_options)
    line_colors = profile_options[ProfileOption.LINE_COLORS]
    line_styles = profile_options[ProfileOption.LINE_STYLES]
    line_widths = profile_options[ProfileOption.LINE_WIDTHS]
    
    y = y + y_shift
    n_solvers = y.shape[0]
    n_runs = y.shape[1]
    y_mean = np.nanmean(y, axis=1).squeeze()

    if profile_options[ProfileOption.ERRORBAR_TYPE] == 'minmax':
        y_lower = np.nanmin(y, axis=1).squeeze()
        y_upper = np.nanmax(y, axis=1).squeeze()
    else:
        y_std = np.std(y, axis=1, ddof=0).squeeze()
        y_lower = y_mean - y_std
        y_upper = y_mean + y_std

    if is_cum:
        y_mean = np.minimum.accumulate(y_mean, axis=1)
        y_lower = np.minimum.accumulate(y_lower, axis=1)
        y_upper = np.minimum.accumulate(y_upper, axis=1)

    # We use block minimization aggregation in case there are too many points, which will cost a lot of time and memory to plot.
    max_eval = y.shape[2]
    n_blocks = 1000
    q = max_eval // n_blocks
    r = max_eval % n_blocks
    blocks = np.full(n_blocks, q)
    blocks[:r] += 1
    # We initialize the cell array to store the x indices and y values to be plotted.
    x_indices = [[] for _ in range(n_solvers)]
    y_values_m = [[] for _ in range(n_solvers)]
    y_values_l = [[] for _ in range(n_solvers)]
    y_values_u = [[] for _ in range(n_solvers)]

    # We only do block minimization when the number of evaluations exceeds n_blocks.
    if max_eval > n_blocks:
        for i_block in range(n_blocks):
            idx_start = int(sum(blocks[:i_block]))
            idx_end = int(idx_start + blocks[i_block])
            idx = np.arange(idx_start, idx_end)
            for i_solver in range(n_solvers):
                i_eval = int(np.max(n_eval[i_solver, :]))
                if profile_options[ProfileOption.HIST_AGGREGATION] == 'min':
                    y_value_m = np.nanmin(y_mean[i_solver, idx])
                    rel_idx = int(np.nanargmin(y_mean[i_solver, idx]))
                    abs_idx = idx[rel_idx]
                    y_value_l = y_lower[i_solver, abs_idx]
                    y_value_u = y_upper[i_solver, abs_idx]
                elif profile_options[ProfileOption.HIST_AGGREGATION] == 'mean':
                    abs_idx = (idx_start + idx_end) // 2
                    y_value_m = np.nanmean(y_mean[i_solver, idx])
                    y_value_l = np.nanmean(y_lower[i_solver, idx])
                    y_value_u = np.nanmean(y_upper[i_solver, idx])
                elif profile_options[ProfileOption.HIST_AGGREGATION] == 'max':
                    y_value_m = np.nanmax(y_mean[i_solver, idx])
                    rel_idx = int(np.nanargmax(y_mean[i_solver, idx]))
                    abs_idx = idx[rel_idx]
                    y_value_l = y_lower[i_solver, abs_idx]
                    y_value_u = y_upper[i_solver, abs_idx]
                if abs_idx > i_eval:
                    continue
                x_indices[i_solver].append(abs_idx)
                y_values_m[i_solver].append(y_value_m)
                y_values_l[i_solver].append(y_value_l)
                y_values_u[i_solver].append(y_value_u)
    else:
        for i_solver in range(n_solvers):
            i_eval = int(np.max(n_eval[i_solver, :]))
            if i_eval > 0:
                x_indices[i_solver] = list(range(min(i_eval, max_eval)))
                y_values_m[i_solver] = y_mean[i_solver, :i_eval].tolist()
                y_values_l[i_solver] = y_lower[i_solver, :i_eval].tolist()
                y_values_u[i_solver] = y_upper[i_solver, :i_eval].tolist()

    # We add the first and last indices if they are not included.
    for i_solver in range(n_solvers):
        i_eval = int(np.max(n_eval[i_solver, :]))
        if i_eval > 0:
            if len(x_indices[i_solver]) == 0:
                x_indices[i_solver] = [0]
                y_values_m[i_solver] = [y_mean[i_solver, 0]]
                y_values_l[i_solver] = [y_lower[i_solver, 0]]
                y_values_u[i_solver] = [y_upper[i_solver, 0]]
            if x_indices[i_solver][0] != 0:
                x_indices[i_solver] = [0] + x_indices[i_solver]
                y_values_m[i_solver] = [y_mean[i_solver, 0]] + y_values_m[i_solver]
                y_values_l[i_solver] = [y_lower[i_solver, 0]] + y_values_l[i_solver]
                y_values_u[i_solver] = [y_upper[i_solver, 0]] + y_values_u[i_solver]
            if x_indices[i_solver][-1] != i_eval - 1 and i_eval != 1:
                x_indices[i_solver].append(i_eval - 1)
                y_values_m[i_solver].append(y_mean[i_solver, i_eval - 1])
                y_values_l[i_solver].append(y_lower[i_solver, i_eval - 1])
                y_values_u[i_solver].append(y_upper[i_solver, i_eval - 1])

    xl_lim = 1 / (problem_n + 1)
    xr_lim = 1 / (problem_n + 1)
    is_log_scale = False

    with plt.rc_context(profile_context):
        for i_solver in range(n_solvers):
            # Truncate the histories according to the function evaluations of each solver.
            i_x = x_indices[i_solver]
            i_y_mean = y_values_m[i_solver]
            i_y_lower = y_values_l[i_solver]
            i_y_upper = y_values_u[i_solver]
            i_eval = len(i_x)
            x = (np.array(i_x) + 1) / (problem_n + 1)
            xr_lim = max(xr_lim, x[-1])

            color = line_colors[i_solver % len(line_colors)]
            line_style = line_styles[i_solver % len(line_styles)]
            line_width = line_widths[i_solver % len(line_widths)]

            if i_eval == 1:
                ax.plot(x, i_y_mean, color=color, marker='o', label=solver_names[i_solver])
            elif i_eval > 1:
                ax.plot(x, i_y_mean, color=color, linestyle=line_style, linewidth=line_width, label=solver_names[i_solver])
                if n_runs > 1:
                    ax.fill_between(x, i_y_lower, i_y_upper, color=color, alpha=0.2)
            if np.any(i_y_mean) and np.any(np.diff(i_y_mean)):
                is_log_scale = True

    # When the function values are not all zero and there is at least some change in the function values, use log scale for the y-axis.
    if is_log_scale:
        ax.set_yscale('log')
    
    ax.yaxis.set_ticks_position('both')
    ax.tick_params(which='both', direction='in')
    ax.legend(loc='upper right')

    if abs(xl_lim - xr_lim) < np.finfo(float).eps:
        xr_lim = xl_lim + np.finfo(float).eps
    
    ax.set_xlim([xl_lim, xr_lim])
    ax.set_xlabel(profile_options[ProfileOption.XLABEL_DATA_PROFILE])
