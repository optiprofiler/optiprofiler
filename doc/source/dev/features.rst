.. py:module:: optiprofiler

Problem features
================

.. currentmodule:: optiprofiler

.. autosummary::
    :toctree: generated/

    Feature

Validation and historical replay
--------------------------------

Magnitude options ``noise_level`` and ``condition_factor`` must be finite,
nonnegative real scalars. ``perturbation_level`` also retains Python's
coordinatewise vector amplitudes: a finite, nonnegative vector of length one
or the problem dimension. Its dimension is checked when the problem is known.
Finite lists/tuples and NumPy vectors follow the same arithmetic. ``mesh_size``
must be finite and strictly positive; ``nan_rate`` must be finite and between
zero and one. Python accepts NumPy real scalars, but not booleans or arrays for
scalar settings. These checks reject an invalid configuration, not the
intentional nonfinite observations produced by a feature such as ``random_nan``.

Execution or replay of an old invalid configuration now fails explicitly. This
includes boolean magnitudes (formerly accepted as integers), negative
``perturbation_level`` and nonfinite magnitudes; no replacement value is guessed.
Loading retained numerical results for reporting does not recreate the feature
and remains separate from replaying its executable configuration.

Python accepts case-insensitive ``mesh_type`` on fresh input and stores its
normalized spelling in the effective specification. Historical native/refined
configurations need separate treatment: older Python code accepted mixed-case
spellings but executed the absolute grid for all of them. The compatibility
reader preserves that actual grid as ``absolute`` instead of silently changing
an archived experiment. It does not rewrite the historical declaration or the
source archive. Fresh ``Feature('quantized', mesh_type='RELATIVE')`` correctly
means the relative grid; old serialized effective ``'RELATIVE'`` retains the
absolute grid that was originally evaluated.

New ``options_user.pkl`` dictionaries carry ``schema='options_user-v2'`` in
addition to the original input. This distinguishes raw current input from old
unmarked flat archives when using ``legacy_compat.replay_arguments``; both the
``feature_name`` shorthand and structured ``feature`` route are supported.
Use ``options_refined.pkl`` to replay the resolved defaults and run count.
