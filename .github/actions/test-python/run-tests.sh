#!/usr/bin/env bash

set -e
set -x

# Run the extra tests on Linux and macOS
if [[ "$RUNNER_OS" == "Linux" || "$RUNNER_OS" == "macOS" ]]; then
  export PYCUTEST_CACHE=$GITHUB_WORKSPACE
  python -m pip install --progress-bar=off .[extra,tests]
  python -m pytest --cov=. --cov-report=xml --run-extra
else
  python -m pip install --progress-bar=off .[tests]
  python -m pytest --cov=. --cov-report=xml
fi
