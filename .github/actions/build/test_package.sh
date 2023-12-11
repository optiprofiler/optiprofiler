#!/usr/bin/env bash

set -e
set -x

for f in dist/*.tar.gz dist/*.whl; do
    python -m venv venv
    source venv/bin/activate
    python -m pip install --progress-bar=off "$f[tests]"
    python -c "import OptiProfiler; OptiProfiler.show_versions()"
    python -m pytest --pyargs OptiProfiler
    deactivate
    rm -r venv
done
