#!/usr/bin/env bash

set -e
set -x

python -m pip install --progress-bar=off .[tests]
if [[ "$OSTYPE" == "darwin"* ]]; then
  python -m pytest --cov=. --cov-report=xml --run-extra
else
  python -m pytest --cov=. --cov-report=xml
fi
