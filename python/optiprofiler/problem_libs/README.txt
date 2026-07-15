# Using Custom Problem Libraries in OptiProfiler

This guide explains the two supported ways to integrate an optimization
problem library with OptiProfiler: an independently installable plugin package,
or an unpackaged local library loaded from an explicit path. The bundled
`custom` library remains a reference implementation for problem construction.

## Recommended: an installable plugin package

An independent package registers its public library name in the
`optiprofiler.problem_libraries` entry-point group. For example:

```toml
[project.entry-points."optiprofiler.problem_libraries"]
my_problems = "optiprofiler_my_problems.plugin:get_problem_library"
```

The target must be a zero-argument factory returning the versioned protocol
object:

```python
from optiprofiler import ProblemLibraryPlugin

from .problems import (
    get_default_options,
    load_problem,
    select_problems,
    validate_options,
)

def get_problem_library():
    return ProblemLibraryPlugin(
        name='my_problems',
        api_version=1,
        select=select_problems,
        load=load_problem,
        check_available=None,
        get_default_options=get_default_options,
        validate_options=validate_options,
    )
```

This API-v1 protocol is currently a development feature and has not yet been
published in an OptiProfiler release. For the future release, plugin metadata
should declare `optiprofiler>=1.4.0` (or a later minimum version if it uses a
newer protocol feature), so an ordinary `pip` installation can resolve a
compatible released core. Until then, plugin CI should test against a checkout
of the development core rather than presenting `1.4.0` as available on a
package index.

The required `select(problem_options, library_options)` callable returns an
ordered sequence of problem-name strings matching the common selection options
and the validated library-specific mapping, not one problem-name string. The
required `load(problem_name, library_options)` callable returns an
OptiProfiler `Problem`. Both callbacks receive the same effective library
options.

A configurable plugin provides both `get_default_options()` and
`validate_options(options)`; providing only one is invalid. The first returns
raw package, environment, or file defaults. OptiProfiler overlays process and
per-run values before calling the second function once to reject unknown or
invalid values and return a normalized mapping. Thus an explicit per-run value
can replace an invalid lower-priority value. Libraries with no configurable
options omit both callbacks. Raw defaults, overrides, and effective options
must use string keys and be pickleable for isolation, experiment storage, and
multiprocessing.

The optional zero-argument `check_available` callable
may raise an informative exception when an external runtime is unavailable;
OptiProfiler calls it once before selecting problems. The optional
`collect_info` callable is reserved for library-maintenance operations and is
not currently called by `benchmark`.

Discovery reads only entry-point metadata. Once selected, the plugin factory is
loaded in the main process. Before each problem is loaded, OptiProfiler
reconstructs the plugin in the current process, which may be a spawned worker.
Keep the factory lightweight and safe to import under multiprocessing `spawn`.
The callbacks must not rely on mutable state shared through one plugin
instance. OptiProfiler resolves library options once in the main process and
passes an isolated copy of the same effective mapping explicitly to the
selector and each worker loader. Treat that mapping as read-only.

The entry-point name, `ProblemLibraryPlugin.name`, and the name in `plibs` must
match exactly. A public name must be an ASCII Python identifier: it must start
with a letter or underscore and contain only letters, digits, and underscores.
Protocol version `1` is currently supported. OptiProfiler
rejects duplicate installed providers, as well as a collision between an
installed provider and a bundled library. This makes resolution independent of
package installation order. An explicit local path is the only override and
takes precedence for that benchmark run.

After installing the package, no path option is needed:

```python
benchmark(solvers, plibs=['my_problems'])
```

## Library-specific experiment options

Each problem library owns its option names and validation. Group options by
library so that unrelated adapters never interpret each other's fields. For
example, using two hypothetical installable plugins:

```python
benchmark(
    solvers,
    plibs=['my_problems', 'other_problems'],
    plib_options={
        'my_problems': {'variant': 'large'},
        'other_problems': {'data_split': 'validation'},
    },
)
```

The precedence is per-run `plib_options`, process-level `set_plib_config`
overrides, library-owned environment or package configuration, then library
defaults. The raw and effective mappings are stored in `options_user.pkl` and
`options_refined.pkl`; the effective mapping is also stored as one
type-preserving serialized value in `data_for_loading.h5`. `plib_options` is
rejected with `problem` and `load`, because those modes do not call a library
adapter.

Runtime installation details such as external executable paths, licenses, and
caches belong in `check_available`, not in experiment options.

## Local development without packaging

If you installed OptiProfiler via pip or conda, the package is installed in a
read-only location such as `site-packages`. Do not add a local library to that
directory. During development, use `custom_problem_libs_path` to specify an
external directory containing your custom problem libraries.

If an external library has the same name as a built-in library, the external
library takes precedence for that benchmark run. This is useful for testing a
new adapter version without modifying the installed OptiProfiler package.

Example:
```python
from optiprofiler import benchmark

benchmark(
    solvers,
    plibs=['my_problems'],  # Your custom library name
    custom_problem_libs_path='/path/to/your/problem_libs',  # External path
    mindim=1,
    maxdim=10
)
```

Where `/path/to/your/problem_libs/` should have this structure:
```
/path/to/your/problem_libs/
├── my_problems/                    # Your custom library
│   ├── my_problems_tools.py       # Must contain my_problems_load() and my_problems_select()
│   └── ...
└── another_lib/                    # Another custom library (optional)
    └── another_lib_tools.py
```

## Overview

OptiProfiler allows benchmarking solvers using custom problem libraries. To use your own problem library, you need to implement two key functions in a tools file:

1. A loading function that retrieves problems by name (`*_load`)
2. A selection function that filters problems based on criteria (`*_select`)

The `custom` library demonstrates one possible implementation approach.

## Detailed Steps

### 1. Create Problem Library Folder

First, create a new subfolder (e.g., `your_problem_lib/`) in any working
directory outside the installed OptiProfiler package:

    /path/to/problem_libs/
    └── your_problem_lib/          <-- Your custom problem library folder
        ├── your_problem_lib_tools.py
        └── python_problems/       <-- Optional organization of problem files

### 2. Implement Core Functions

#### 2.1 Tools File (refer to `custom_tools.py`)

Create a python file named `your_problem_lib_tools.py` inside your library folder. This file must contain two functions:

The tools file and functions must use the library name exactly. OptiProfiler
does not infer the library name from other `*_tools.py` files. For example, a
library folder named `solar` must contain `solar_tools.py` with
`solar_load` and `solar_select`.

1.  `your_problem_lib_load(problem_name)`:
    -   Accepts a problem name as input.
    -   Returns a valid `Problem` class instance (from `optiprofiler.opclasses`).

2.  `your_problem_lib_select(options)`:
    -   Accepts a dictionary `options` containing filtering criteria (e.g., `mindim`, `maxdim`, `ptype`, etc.).
    -   Returns a list of problem names that satisfy the criteria.

These original one-argument functions remain compatible when the adapter has
no explicit `plib_options`. A configurable local adapter additionally defines
`your_problem_lib_get_default_options()` and
`your_problem_lib_validate_options(options)`. Its selector then accepts
`(problem_options, library_options)`, and its loader accepts
`library_options` as a keyword argument. A legacy local adapter without those
callbacks cannot accept `plib_options` or provider-scoped `set_plib_config`
overrides; configure it externally or adopt the explicit callbacks.

The internal implementation is flexible. The example in `custom_tools.py` demonstrates an approach where problem information is stored in a CSV file (`probinfo_python.csv`) for efficient filtering.

### 3. Create Problem Definitions

You have complete flexibility in how you organize and define your optimization problems. Some options include:

- Individual Python files (one per problem), as seen in `python_problems/`.
- A single file containing multiple problem definitions.
- Problem definitions stored in other formats.

The only requirement is that your `your_problem_lib_load` function can retrieve these problems by name and return them as `Problem` objects.

### 4. Using Your Custom Problem Library

In your benchmarking script, use your custom problem library by setting the `plibs` option in `problem_options`:

```python
from optiprofiler import benchmark

# Set options
problem_options = {
    'plibs': ['your_problem_lib'],  # Use your custom library
    'mindim': 1,
    'maxdim': 100
}

# Run benchmark
benchmark(solvers, **problem_options)
```

When the library is outside the installed package, also pass its parent folder
or its own folder:

```python
benchmark(
    solvers,
    plibs=['your_problem_lib'],
    custom_problem_libs_path='/path/to/problem_libs',
)
```

## Problem Class Properties

When creating Problem objects, you have to set at least two properties:

- `fun`: The objective function callable
- `x0`: The initial guess for the optimization variable

You can create the Problem object using the constructor:

```python
from optiprofiler.opclasses import Problem
problem = Problem(fun=your_objective_function, x0=initial_guess)
```

## Example Implementation

The `custom` folder provides a reference implementation that you can study and adapt for your own problem library. It demonstrates:

1. One way to organize problem files
2. How to implement the required functions

You are encouraged to examine the files in the `custom` folder to understand the implementation details and adapt them to your specific needs.

You may also want to view our website for more information on how to use OptiProfiler: www.optprof.com

## Built-in adapter note: SOLAR

The `solar` problem library wraps the upstream SOLAR black-box
optimization simulator. It vendors a slim runtime under `runtime/solar/`,
including the upstream LGPL-2.1 license, README, and manifest recording the
exact upstream commit. The adapter builds a local `solar` executable on first
use when the binary is missing. Some SOLAR problems are much slower than
ordinary algebraic test problems, so use small `max_eval_factor` values when
trying it for the first time.
