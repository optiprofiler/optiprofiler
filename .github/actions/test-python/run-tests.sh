#!/usr/bin/env bash

set -e
set -x

python -m pip install --progress-bar=off ".[tests]"
# Measure only the core package. External problem libraries own their adapter
# and runtime tests and are exercised separately by integration CI.
python -m pytest --cov=optiprofiler --cov-report=xml
