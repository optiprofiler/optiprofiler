# OptiProfiler Repository Refactoring Plan

This document outlines the steps to restructure the OptiProfiler repository to better support both MATLAB and Python environments, ensuring proper packaging and dependency management.

## Phase 1: Preparation & Directory Structure

1.  **Create Branch**
    *   Ensure you are on a new branch (e.g., `reorganize`) to avoid affecting the main codebase during the process.

2.  **Create Target Directories**
    *   Create `matlab/optiprofiler/problem_libs` to house MATLAB-specific problem sets.
    *   Create `python/optiprofiler/problem_libs` to house Python-specific problem sets and wrappers.

## Phase 2: Git Submodule Migration

*Note: Submodules cannot be simply moved; they must be removed and re-added.*

3.  **Migrate `matcutest` (MATLAB/Linux only)**
    *   Remove existing submodule: `git rm problems/matcutest`
    *   **Note**: We will NOT add it as a submodule. Instead, `setup.m` will dynamically clone it if needed.

4.  **Migrate `s2mpj` (Shared Library)**
    *   Remove existing submodule: `git rm problems/s2mpj`
    *   Add to new location (inside `src` subdirectory):
        ```bash
        git submodule add https://github.com/GrattonToint/S2MPJ.git python/optiprofiler/problem_libs/s2mpj/src
        ```

## Phase 3: Code Migration

5.  **Migrate `custom` Problems**
    *   Move MATLAB custom problems from `problems/custom/...` to `matlab/optiprofiler/problem_libs/custom`.
    *   Move Python custom problems from `problems/custom/...` to `python/optiprofiler/problem_libs/custom`.

6.  **Migrate `pycutest` Wrapper**
    *   Move `problems/pycutest` content to `python/optiprofiler/problem_libs/pycutest` (or as a module).

7.  **Cleanup**
    *   Remove the now empty `problems` directory in the root.

## Phase 4: Update Configuration & Scripts

8.  **Update MATLAB `setup.m`**
    *   **Path Updates**:
        *   Update `matcutest` path to `matlab/optiprofiler/problem_libs/matcutest/src`.
    *   **S2MPJ Handling (Copy Strategy)**:
        *   Source: `../../python/optiprofiler/problem_libs/s2mpj/src`.
        *   Destination: `matlab/optiprofiler/problem_libs/s2mpj`.`.
    *   **S2MPJ Handling (Copy Strategy)**:
        *   Source: `../../python/optiprofiler/problem_libs/s2mpj/src`.
        *   Destination: `matlab/optiprofiler/problem_libs/s2mpj`.
        *   **Logic**: In `setup.m`, copy the contents from Source to Destination. This ensures MATLAB has its own local copy of the library.
        *   Add Destination to MATLAB path.
    *   **Installation Logic**:
        *   Implement check: `if isunix && ~ismac`.
        *   **Dynamic Clone**: Check if `matcutest` directory exists. If not, use `system('git clone ...')` to download it to `matlab/optiprofiler/problem_libs/matcutest the `options` structure.

9.  **Update Python Configuration**
    *   Ensure `python/optiprofiler/problem_libs/__init__.py` exists.
    *   Verify imports work correctly with the new structure.

## Phase 5: Verification

10. **Verify MATLAB**
    *   Run `setup.m` and check if paths (especially the cross-referenced `s2mpj`) are added correctly.
    *   Run `testOptiProfiler`.

11. **Verify Python**
    *   Run `pip install -e .` in the `python` directory.
    *   Verify `import optiprofiler.problem_libs.s2mpj` works.
