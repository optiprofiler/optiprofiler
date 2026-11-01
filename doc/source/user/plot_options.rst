Native plot styling
===================

``plot_options`` changes presentation, not the feature, seed, numerical
histories, convergence tests, scores, or EvalReport plot coordinates.
Its default is empty, so existing figures retain their default style.
Use the same option when loading saved experiments to redraw them without
running solvers or modifying the source archive.

Python uses native Matplotlib property names:

.. code-block:: python

    benchmark(solvers, plibs=['s2mpj'], plot_options={
        'xlabel': {'fontsize': 18},
        'ylabel': {'fontsize': 18},
        'axes': {'facecolor': 'white'},
        'legend': {'fontsize': 12, 'loc': 'lower right'},
        'tick_params': {'labelsize': 14},
        'rc': {'text.usetex': False},
    })

The groups ``figure``, ``axes``, ``xaxis``, ``yaxis``, ``xlabel``, ``ylabel``,
``title``, ``text``, and ``line`` use the corresponding object's setters.
Use properties accepted by the installed object's ``set`` method, rather
than arbitrary Python attributes. ``legend`` uses Matplotlib legend
constructor properties; ``tick_params`` uses ``Axes.tick_params`` arguments.
``rc`` overrides the package's rendering context locally; it does not change
global defaults. Figure canvas/number, rendering backend/interactive mode,
and legend handles/labels remain managed by OptiProfiler. Use ``solver_names``
to name solver curves. Parallel history rendering requires serializable values.
Specific label/title groups and explicit legend font/color settings take
precedence over the general ``text`` group.

MATLAB uses native graphics property names instead:

.. code-block:: matlab

    options.plot_options = struct( ...
        'axes', struct('FontSize',14), ...
        'xlabel', struct('FontSize',18), ...
        'ylabel', struct('FontSize',18), ...
        'legend', struct('FontSize',12,'Location','southeast'));
    benchmark(solvers, options);

MATLAB supports the same object groups except ``tick_params`` and ``rc``;
use Axes or ruler properties for ticks. Names are checked before solver
execution without creating figures. Values of every group except ``figure``
are also checked before solver execution, by MATLAB itself on unparented
objects of the same classes (no figure or window is created). Creating a
figure would open a window, so of the ``figure`` values only ``Position``,
``OuterPosition``, ``InnerPosition``, ``PaperPosition``, ``PaperSize``,
``Color`` and ``Name`` are checked before solver execution, against the rules
MATLAB applies to them; the other ``figure`` values and object-valued entries
such as ``Parent`` are checked when they are applied. If such a value is
refused while a profile is exported, ``benchmark`` raises the error after
closing every figure it opened, the report's ``rendering`` stage fails, and the
saved numerical data remain loadable. A history plot uses the SVG fallback
instead of raising, so a direct Problem finishes with a ``partial`` report.
Styling is applied after the package defaults to
single profiles, combined summaries, and history figures. Properties apply
to every matching object, so a common explicit axis limit or position also
applies to every panel. Portable SVG/HTML output, used when native graphics
are unavailable, retains its portable style and reports that native property
maps were not applied. Property availability follows the installed plotting
library/MATLAB version; the two languages do not translate property names.

Axis label text
---------------

The text of ``xlabel_data_profile``, ``xlabel_log_ratio_profile`` and
``xlabel_performance_profile`` is passed to the renderer verbatim, so write
``\%`` for a percent sign under LaTeX. MATLAB uses its LaTeX interpreter
unless the ``xlabel`` or ``text`` group sets another one: a bare ``%`` is
invalid syntax there, and the label is drawn as raw text. Python uses LaTeX
when it is installed, unless the ``rc`` group sets ``text.usetex`` to false:
under LaTeX a bare ``%`` starts a comment and cuts the label, and without LaTeX
``\%`` keeps its backslash. ``$50\%$`` works in both cases.

In ``ylabel_data_profile``, ``ylabel_log_ratio_profile`` and
``ylabel_performance_profile``, ``%s`` is replaced by the current tolerance.
Keep it inside ``$...$``, as the defaults do, because the tolerance is LaTeX
math. Outside math mode MATLAB draws the label as raw text, and Python shows
literal text such as ``10^{-1}`` without LaTeX and raises a LaTeX error under
LaTeX.
