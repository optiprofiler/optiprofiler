#!/usr/bin/env bash

set -e
set -x

# Run the extra tests on Linux only (pycutest not available on macOS ARM64 in CI)
if [[ "$RUNNER_OS" == "Linux" ]]; then
  export PYCUTEST_CACHE=$GITHUB_WORKSPACE
  python -m pip install --progress-bar=off ".[extra,tests]"
  # Measure only the optiprofiler package. Using --cov=. from the repo root
  # includes thousands of unrelated .py files (submodules, embedded problems),
  # which makes coverage XML generation extremely slow on Windows CI.
  python -m pytest --cov=optiprofiler --cov-report=xml --run-extra
else
  python -m pip install --progress-bar=off ".[tests]"
  python -m pytest --cov=optiprofiler --cov-report=xml
fi
