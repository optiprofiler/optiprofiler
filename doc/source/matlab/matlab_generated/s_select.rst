.. _matsselect:

s_select
========

**s_select(**\ *options*\ **)**
    **s_select** selects the problems in S2MPJ that satisfy given criteria.

--------------------------------------------------------------------------

**problem_names** = **s_select**\(**options**) returns the names of selected problems from S2MPJ that satisfy the criteria in **options** as a cell array **problem_names**. More details about S2MPJ can be found in the website https://github.com/GrattonToint/S2MPJ.

**options** is a struct with the following fields:

    - **ptype**: the type of the problems to be selected. It should be a string or char consisting of any combination of ``'u'`` (unconstrained), ``'b'`` (bound constrained), ``'l'`` (linearly constrained), and ``'n'`` (nonlinearly constrained), such as ``'b'``, ``'ul'``, ``'ubn'``. Default is ``'ubln'``.

    - **mindim**: the minimum dimension of the problems to be selected. Default is ``1``.

    - **maxdim**: the maximum dimension of the problems to be selected. Default is ``Inf``.

    - **minb**: the minimum number of bound constraints of the problems to be selected. Default is ``0``.

    - **maxb**: the maximum number of bound constraints of the problems to be selected. Default is ``Inf``.

    - **mincon**: the minimum number of linear and nonlinear constraints of the problems to be selected. Default is ``0``.

    - **maxcon**: the maximum number of linear and nonlinear constraints of the problems to be selected. Default is ``Inf``.

    - **oracle**: the oracle provided by the problem. If it is 0, then the problem should provide zeroth-order information. If it is 1, then the problem should provide both zeroth-order and first-order information. If it is 2, then the problem should provide zeroth-order, first-order, and second-order information. Default is ``0``.

    - **excludelist**: the list of problems to be excluded. Default is not to exclude any problem.

Two things to note:

1. All the information about the problems can be found in a csv file named 'probinfo.csv' in the same directory as this function.

2. The problem name may appear in the form of 'problem_name_dim_m_con' where 'problem_name' is the name of the problem, 'dim' is the dimension of the problem, and 'm_con' is the number of linear and nonlinear constraints of the problem. This case only happens when this problem can accept extra arguments to change the dimension or the number of constraints. This information is stored in the 'probinfo.csv' file as the last few columns.