#!/usr/bin/env bash

set -e
set -x

python -m pip install --progress-bar=off ".[tests]"
python -m pip install --progress-bar=off "pytest-timeout>=2.2,<3"
# Measure only the core package. External problem libraries own their adapter
# and runtime tests and are exercised separately by integration CI.
# A result timeout inside a worker test does not bound pool startup/teardown.
# Print each case and dump all thread stacks on a hang, including on Windows,
# rather than consume the runner's six-hour default without useful evidence.
python -m pytest -v --timeout=180 --timeout-method=thread --cov=optiprofiler --cov-report=xml
