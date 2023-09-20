#!/usr/bin/env bash

set -e
set -x


if [[ "$OSTYPE" == "darwin"* ]]; then
  python -m pip install --progress-bar=off .[extra,tests]
  python -m pytest --cov=. --cov-report=xml --run-extra
else
  python -m pip install --progress-bar=off .[tests]
  python -m pytest --cov=. --cov-report=xml
fi
