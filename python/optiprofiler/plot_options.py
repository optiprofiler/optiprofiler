"""Native Matplotlib styling, kept outside numerical/scoring configuration."""
import copy
from collections.abc import Mapping

import matplotlib
from matplotlib.figure import Figure
from matplotlib.legend import Legend
from matplotlib.text import Text


_GROUPS = {'figure', 'axes', 'xaxis', 'yaxis', 'xlabel', 'ylabel', 'title',
           'text', 'line', 'legend', 'tick_params', 'rc'}


def normalize_plot_options(value):
    """Own and validate native property maps before any solver executes."""
    if not isinstance(value, Mapping):
        raise TypeError('plot_options must be a mapping of object groups to native property mappings.')
    unknown = set(value) - _GROUPS
    if unknown:
        raise ValueError(f'Unknown plot_options group(s): {sorted(unknown)!r}. Allowed groups: {sorted(_GROUPS)!r}.')
    if any(not isinstance(properties, Mapping) for properties in value.values()):
        raise TypeError('Each plot_options group must be a native property mapping.')
    result = copy.deepcopy({group: dict(properties) for group, properties in value.items()})
    if not result:
        return result
    # Backend selection is runtime policy, not styling. In particular a GUI
    # backend must not be initialized in file-only plotting workers, and
    # rc_context does not restore the backend on exit.
    if set(result.get('rc', {})) & {'backend', 'backend_fallback', 'interactive'}:
        raise ValueError('plot_options.rc cannot change the managed rendering backend or interactive mode.')
    if set(result.get('figure', {})) & {'canvas', 'number'}:
        raise ValueError('plot_options.figure cannot replace the managed figure canvas or number.')
    # The benchmark owns solver identities and proxy handles. Letting a style
    # map supply these constructor arguments would fail only during export or
    # silently label another solver's curve as the current solver.
    if set(result.get('legend', {})) & {'handles', 'labels'}:
        raise ValueError('plot_options.legend cannot replace solver handles or labels; use solver_names instead.')
    figure = Figure()
    try:
        with matplotlib.rc_context(result.get('rc', {})):
            axes = figure.subplots()
            axes.plot([0, 1], [0, 1], label='validation')
            axes.legend(**result.get('legend', {}))
            apply_plot_options(figure, result)
    except (AttributeError, TypeError, ValueError) as cause:
        raise ValueError(f'Invalid native property in plot_options: {cause}') from cause
    finally:
        figure.clear()
    return result


def apply_plot_options(figure, options):
    """Apply user properties after package defaults, without touching scores.

    Groups use the installed Matplotlib version's native setters. Legend
    constructor options are applied by the legend factory, not by recreating
    it here (which could lose proxy handles or literal underscore labels).
    """
    if not options:
        return
    figure.set(**options.get('figure', {}))
    # Legends are constructed before generic text styling, but their explicit
    # font/color choices must still win (as do xlabel/ylabel/title groups).
    # Retain the resolved native values instead of rebuilding solver handles.
    legend_options = options.get('legend', {})
    legend_styles = []
    for legend in figure.findobj(match=Legend):
        for text in legend.get_texts():
            style = {}
            if 'prop' in legend_options:
                style['fontproperties'] = text.get_fontproperties().copy()
            elif 'fontsize' in legend_options:
                style['fontsize'] = text.get_fontsize()
            if 'labelcolor' in legend_options:
                style['color'] = text.get_color()
            legend_styles.append((text, style))
        title = legend.get_title()
        if 'title_fontproperties' in legend_options:
            legend_styles.append((title, {'fontproperties': title.get_fontproperties().copy()}))
        elif 'title_fontsize' in legend_options:
            legend_styles.append((title, {'fontsize': title.get_fontsize()}))
    for text in figure.findobj(match=Text):
        text.set(**options.get('text', {}))
        if text.get_gid() == 'optiprofiler-title':
            text.set(**options.get('title', {}))
        elif text.get_gid() == 'optiprofiler-ylabel':
            text.set(**options.get('ylabel', {}))
    for axes in figure.axes:
        axes.set(**options.get('axes', {}))
        axes.xaxis.set(**options.get('xaxis', {}))
        axes.yaxis.set(**options.get('yaxis', {}))
        axes.xaxis.label.set(**options.get('xlabel', {}))
        axes.yaxis.label.set(**options.get('ylabel', {}))
        axes.title.set(**options.get('title', {}))
        axes.tick_params(**options.get('tick_params', {}))
        for line in axes.lines:
            line.set(**options.get('line', {}))
    for text, style in legend_styles:
        text.set(**style)
