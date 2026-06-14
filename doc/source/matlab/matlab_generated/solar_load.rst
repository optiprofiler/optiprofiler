.. _matsolarload:

solar_load
=================

**solar_load(**\ *problem_name*\ **)**
    **solar_load** converts a SOLAR problem name to a Problem class instance.

--------------------------------------------------------------------------

**problem** = **solar_load**\(**problem_name**) returns a **Problem**
class instance **problem** corresponding to a selected scalar SOLAR problem.
More details about the upstream SOLAR simulator can be found at
<https://github.com/bbopt/solar>.

The SOLAR MATLAB adapter vendors a slim SOLAR runtime under LGPL-2.1. The
runtime keeps its upstream license, README, and manifest under
``runtime/solar/``. Some SOLAR problems call a comparatively expensive external
C++ simulator; use small evaluation budgets when testing.

You may use **solar_select** to get the problem names you want.
