.. _matsolarmatlabselect:

solar_matlab_select
===================

**solar_matlab_select(**\ *options*\ **)**
    **solar_matlab_select** selects SOLAR problems satisfying given criteria.

--------------------------------------------------------------------------

**problem_names** = **solar_matlab_select**\(**options**) returns the names of
selected SOLAR problems as a cell array. More details about the upstream SOLAR
simulator can be found at <https://github.com/bbopt/solar>.

**options** is a struct with the following fields:

    - **ptype**: problem type, using any combination of ``'u'``, ``'b'``,
      ``'l'``, and ``'n'``. Default is ``'ubln'``.
    - **mindim** and **maxdim**: minimum and maximum dimension.
    - **minb** and **maxb**: minimum and maximum number of bound constraints.
    - **minlcon** and **maxlcon**: minimum and maximum number of linear
      constraints.
    - **minnlcon** and **maxnlcon**: minimum and maximum number of nonlinear
      constraints.
    - **mincon** and **maxcon**: minimum and maximum number of linear and
      nonlinear constraints.
    - **excludelist**: list of problem names to exclude.

SOLAR 8 and 9 are multiobjective and are not returned by the first scalar
OptiProfiler selector. SOLAR 11 is disabled for now because the upstream
runtime returns an empty output at the documented initial point in the current
recorded revision.

The SOLAR MATLAB adapter vendors a slim SOLAR runtime under LGPL-2.1. The
runtime keeps its upstream license, README, and manifest under
``runtime/solar/``. Some SOLAR problems call a comparatively expensive external
C++ simulator; use small evaluation budgets when testing.
