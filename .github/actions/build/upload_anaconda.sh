#!/usr/bin/env bash

set -e
set -x

export PATH=$CONDA/bin:$PATH
conda create -n venv -y python="$(python -c 'import platform; print(platform.python_version())')"
source activate venv
conda install -y anaconda-client
# The build also retains a release-manifest.json acceptance receipt. Upload
# distributions only, not that auxiliary metadata file.
anaconda -t "$ANACONDA_TOKEN" upload --force -u optiprofiler dist/*.whl dist/*.tar.gz
