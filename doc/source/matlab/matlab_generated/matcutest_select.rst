.. _matmatcutestselect:

matcutest_select
================

**matcutest_select(**\ *options*\ **)**
    **matcutest_select** selects the problems in MatCUTEst that satisfy given criteria.

--------------------------------------------------------------------------

**problem_names** = **matcutest_select**\(**options**) returns the names of selected problems from MatCUTEst that satisfy the criteria in **options** as a cell array **problem_names**. More details about MatCUTEst can be found in the official website: <https://github.com/matcutest>.

**options** is a struct with the following fields:

    - **ptype**: the type of the problems to be selected. It should be a string or char consisting of any combination of ``'u'`` (unconstrained), ``'b'`` (bound constrained), ``'l'`` (linearly constrained), and ``'n'`` (nonlinearly constrained), such as ``'b'``, ``'ul'``, ``'ubn'``. Default is ``'ubln'``.

    - **mindim**: the minimum dimension of the problems to be selected. Default is ``1``.

    - **maxdim**: the maximum dimension of the problems to be selected. Default is ``Inf``.

    - **minb**: the minimum number of bound constraints of the problems to be selected. Default is ``0``.

    - **maxb**: the maximum number of bound constraints of the problems to be selected. Default is ``Inf``.

    - **mincon**: the minimum number of linear and nonlinear constraints of the problems to be selected. Default is ``0``.

    - **maxcon**: the maximum number of linear and nonlinear constraints of the problems to be selected. Default is ``Inf``.

    - **excludelist**: the list of problems to be excluded. Default is not to exclude any problem.

Note that MatCUTEst is only available in Linux.