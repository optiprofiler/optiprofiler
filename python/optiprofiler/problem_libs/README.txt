# Using Custom Problem Libraries in OptiProfiler

This guide explains how to create and integrate your own optimization problem library into OptiProfiler using the `custom` as a reference.

## Quick Start: The Simplest Way to Create and Use Custom Problems

If you want to quickly create and use your own optimization problems, follow these steps:

1. Create Your Problem Files:
   - Use the existing `.py` files (`custom1.py`, `custom2.py`, `custom3.py`, `custom4.py`) in the `optiprofiler/problem_libs/custom/python_problems/` folder as examples to define your own problems. Each file should define a function that returns a dictionary containing the problem attributes (fun, x0, etc.).
   - Save your new `.py` files in the `optiprofiler/problem_libs/custom/python_problems/` folder.
   - You may delete the existing example files if you don't need them.

2. Run `custom_get_info`:
   - You can run the `custom_get_info` function from `optiprofiler.problem_libs.custom.custom_tools` to update the problem information.
   - This function will automatically generate/update the `probinfo_python.csv` file in the `optiprofiler/problem_libs/custom/` folder.
   - Example code to run it:
     ```python
     from optiprofiler.problem_libs.custom.custom_tools import custom_get_info
     custom_get_info()
     ```

That's it! You can now use your custom problems in OptiProfiler by specifying `problem_options={'plibs': ['custom']}` in your benchmarking script.

## Overview

OptiProfiler allows benchmarking solvers using custom problem libraries. To use your own problem library, you need to implement two key functions in a tools file:

1. A loading function that retrieves problems by name (`*_load`)
2. A selection function that filters problems based on criteria (`*_select`)

The `custom` library demonstrates one possible implementation approach.

## Detailed Steps

### 1. Create Problem Library Folder

First, create a new subfolder (e.g., `your_problem_lib/`) within the `optiprofiler/problem_libs/` folder:

    optiprofiler/
    ├── problem_libs/
    │   ├── custom/               <-- Example implementation
    │   │   ├── custom_tools.py   <-- Contains load and select functions
    │   │   ├── probinfo_python.csv
    │   │   └── python_problems/  <-- Python implementation of example problems
    │   ├── your_problem_lib/     <-- Your custom problem library folder
    │   │   ├── your_problem_lib_tools.py
    │   │   └── python_problems/  <-- Your problem files
    │   ├── s2mpj/                <-- Built-in problem library
    │   └── ...
    └── ...

### 2. Implement Core Functions

#### 2.1 Tools File (refer to `custom_tools.py`)

Create a python file named `your_problem_lib_tools.py` inside your library folder. This file must contain two functions:

1.  `your_problem_lib_load(problem_name)`:
    -   Accepts a problem name as input.
    -   Returns a valid `Problem` class instance (from `optiprofiler.opclasses`).

2.  `your_problem_lib_select(options)`:
    -   Accepts a dictionary `options` containing filtering criteria (e.g., `mindim`, `maxdim`, `ptype`, etc.).
    -   Returns a list of problem names that satisfy the criteria.

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