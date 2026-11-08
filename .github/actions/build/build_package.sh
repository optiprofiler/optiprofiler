#!/usr/bin/env bash

set -e
set -x

python -m venv venv
source venv/bin/activate
python -m pip install --progress-bar=off build
python -m build
deactivate
rm -r venv
# Both artifacts must carry the bundled S2MPJ revision of this checkout's gitlink.
python tools/check_release.py --dist-dir dist
