#!/usr/bin/env bash

set -e
set -x

# Runners no longer provide conda (ubuntu-latest is Ubuntu 26.04), so install
# the client from PyPI. 1.13 needs the removed pkg_resources, and 1.16 turns on
# telemetry that reports the first arguments of each command to Anaconda; it
# is disabled here in any case.
python -m venv venv
source venv/bin/activate
python -m pip install --progress-bar=off 'anaconda-client==1.15.0'
export ANACONDA_TELEMETRY_ENABLED=false
# The client reads the token from ANACONDA_API_TOKEN rather than from -t, which
# xtrace printed (and telemetry would report); export it without xtrace.
set +x
export ANACONDA_API_TOKEN="${ANACONDA_TOKEN:?is not set}"
set -x
# The build also retains a release-manifest.json acceptance receipt. Upload
# distributions only, not that auxiliary metadata file. A distribution that
# already exists is an error, as on PyPI (skip-existing: false in build.yml):
# --force would replace it, so that a re-run could leave the two channels with
# different files for one version.
anaconda upload -u optiprofiler dist/*.whl dist/*.tar.gz
