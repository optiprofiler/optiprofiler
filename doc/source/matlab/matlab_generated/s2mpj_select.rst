.. _mats2mpjselect:

s2mpj_select
============

**s2mpj_select(**\ *options*\ **)**
    **s2mpj_select** selects the problems in S2MPJ that satisfy given criteria.

--------------------------------------------------------------------------

**problem_names** = **s2mpj_select**\(**options**) returns the names of selected problems from S2MPJ that satisfy the criteria in **options** as a cell array **problem_names**. More details about S2MPJ can be found in the official website <https://github.com/GrattonToint/S2MPJ>.

**options** is a struct with the following fields:

    - **ptype**: the type of the problems to be selected. It should be a string or char consisting of any combination of ``'u'`` (unconstrained), ``'b'`` (bound constrained), ``'l'`` (linearly constrained), and ``'n'`` (nonlinearly constrained), such as ``'b'``, ``'ul'``, ``'ubn'``. Default is ``'ubln'``.

    - **mindim**: the minimum dimension of the problems to be selected. Default is ``1``.

    - **maxdim**: the maximum dimension of the problems to be selected. Default is ``Inf``.

    - **minb**: the minimum number of bound constraints of the problems to be selected. Default is ``0``.

    - **maxb**: the maximum number of bound constraints of the problems to be selected. Default is ``Inf``.

    - **minlcon**: the minimum number of linear constraints of the problems to be selected. Default is ``0``.

    - **maxlcon**: the maximum number of linear constraints of the problems to be selected. Default is ``Inf``.

    - **minnlcon**: the minimum number of nonlinear constraints of the problems to be selected. Default is ``0``.

    - **maxnlcon**: the maximum number of nonlinear constraints of the problems to be selected. Default is ``Inf``.

    - **mincon**: the minimum number of linear and nonlinear constraints of the problems to be selected. Default is ``0``.

    - **maxcon**: the maximum number of linear and nonlinear constraints of the problems to be selected. Default is ``Inf``.

    - **oracle**: the oracle provided by the problem. If it is 0, then the problem should provide zeroth-order information. If it is 1, then the problem should provide both zeroth-order and first-order information. If it is 2, then the problem should provide zeroth-order, first-order, and second-order information. Default is ``0``.

    - **excludelist**: the list of problems to be excluded. Default is not to exclude any problem.

Three things to note:

1. All the information about the problems can be found in a csv file named ``probinfo.csv`` in the same directory as this function.

2. The problem name may appear in the form of 'problem_name_dim_m_con' where 'problem_name' is the name of the problem, 'dim' is the dimension of the problem, and 'm_con' is the number of linear and nonlinear constraints of the problem. This case only happens when this problem can accept extra arguments to change the dimension or the number of constraints. This information is stored in the ``probinfo.csv`` file as the last few columns.

3. There is a file ``variable_size.txt`` in the same directory as this function. This file can be used to set the ``variable_size`` option to ``'default'``, ``'min'``, ``'max'``, or ``'all'`` (without quotes in the file). If this file does not exist or is empty, the ``variable_size`` option will be set to ``'default'``. ``variable_size`` is used to determine how to select the problems with variable dimension and/or number of constraints. The options are:

    - 'default': Only consider the default dimension and number of constraints for each problem.
    - 'min':     For problems with variable dimension and/or constraints, select the one with the smallest dimension and, among those, the smallest number of constraints that satisfies the options (priority: smaller dimension, then fewer constraints).
    - 'max':     For problems with variable dimension and/or constraints, select the one with the largest dimension and, among those, the largest number of constraints that satisfies the options (priority: larger dimension, then more constraints).
    - 'all':     For problems with variable dimension and/or constraints, include all configurations that satisfy the options.