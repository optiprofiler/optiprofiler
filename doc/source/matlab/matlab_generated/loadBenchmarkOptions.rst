.. _matloadbenchmarkoptions:

loadBenchmarkOptions
====================

``[options, receipt] = loadBenchmarkOptions(source, legacy_feature_name)``
    Imports trusted native configuration for a new benchmark. The second
    argument is optional and only supplies identity missing from a legacy file.

``source`` is a path to a MAT file containing ``options_refined`` or a loaded
scalar native options structure. The returned ``options`` has one canonical
``feature`` and a separate experiment ``n_runs`` when retained. It excludes
saved load selectors, output targets and envelope-only fields. The caller must
provide solvers separately and choose the new output destination.

For new ``options_refined-v2`` archives, the retained native Feature is the
canonical authority. ``feature_specification`` is its inspection duplicate;
disagreement is an error rather than an implicit choice. Canonical effective
options and known/unknown declaration provenance survive native transport.
A specification-only v2 input describes a new replay request, not a recovered
historical declaration.

Supported unversioned native configurations use the narrow legacy Feature
import boundary. A retained root ``n_runs`` overrides an old child residue.
Old flat settings without feature identity require ``legacy_feature_name``;
the loader never infers identity from a path, stamp, or collection of options.
The ``receipt`` records retained facts and any explicit present override,
without claiming that an old resolved count was the original user request.

This is **not** the ``benchmark`` ``load`` option: it prepares a fresh execution,
not history reanalysis or a resumed random-stream/runtime checkpoint. Native
MAT loading is for trusted inputs and may execute code. Original callback/class
dependencies must be available; JSON callback descriptions are not executable
reconstruction recipes. Unsupported versions and inconsistent native states
raise explicit errors.

See :ref:`use` for examples of fresh replay and saved-result reanalysis.
