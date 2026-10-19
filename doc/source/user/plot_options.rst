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
execution without creating figures; native property values are checked by
MATLAB when applying them. Styling is applied after the package defaults to
single profiles, combined summaries, and history figures. Properties apply
to every matching object, so a common explicit axis limit or position also
applies to every panel. Portable SVG/HTML output, used when native graphics
are unavailable, retains its portable style and reports that native property
maps were not applied. Property availability follows the installed plotting
library/MATLAB version; the two languages do not translate property names.
