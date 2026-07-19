.. _matlab_custom_problem_libraries:

Custom MATLAB problem libraries
===============================

A custom provider can live in any user-owned writable directory.  It must
define a selector and loader with stable MATLAB function names:

* ``myproblems_select(options)`` returns the matching problem names;
* ``myproblems_load(name)`` returns a :ref:`Problem <matproblem>`.

Register the root and canonical functions once:

.. code-block:: matlab

    registration = struct( ...
        'name', 'myproblems', ...
        'root', '/path/to/myproblems', ...
        'select_function', 'myproblems_select', ...
        'load_function', 'myproblems_load');
    registerProblemLibrary(registration);

Then select it like a built-in provider:

.. code-block:: matlab

    options.plibs = {'s2mpj', 'myproblems'};
    scores = benchmark({@solver1, @solver2}, options);

An identical repeated registration is harmless.  The same public name cannot
silently be redirected to a different root or different functions; unregister
the old provider before registering a replacement.

.. code-block:: matlab

    unregisterProblemLibrary('myproblems');

Unregistering removes the registry and MATLAB path entries but preserves the
provider directory and its data.  Generic tools may call
``resolveProblemLibrary('myproblems')`` to obtain validated function handles.
Ordinary benchmark users need only ``options.plibs``.
