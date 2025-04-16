.. _matsload:

s_load
=======

**s_load(**\ *problem_name*\ **)**
    **s_load** coverts a problem in S2MPJ to a Problem class instance.

--------------------------------------------------------------------------

Users only need to use the following signature to call this function:

**problem** = **s_load**\(**problem_name**) returns a **Problem** class instance **problem** that corresponds to the problem named **problem_name** in S2MPJ. More details about S2MPJ can be found in the website https://github.com/GrattonToint/S2MPJ.

There are two ways to get **problem_name** you want.

1. Use the function **s_select** to get the problem names you want.

2. Look for a csv file named 'probinfo.csv' in the same directory as this function. The csv file contains the information of all the problems in S2MPJ.

Note that problem name may appear in the form of 'problem_name_dim_m_con' where 'problem_name' is the name of the problem, 'dim' is the dimension of the problem, and 'm_con' is the number of linear and nonlinear constraints of the problem. This case only happens when this problem can accept extra arguments to change the dimension or the number of constraints. This information is stored in the 'probinfo.csv' file as the last few columns.