import sys
from pathlib import Path

# PEP0440 compatible formatted version, see:
# https://www.python.org/dev/peps/pep-0440/
#
# Final release markers:
#   X.Y.0   # For first release after an increment in Y
#   X.Y.Z   # For bugfix releases
#
# Admissible pre-release markers:
#   X.YaN   # Alpha release
#   X.YbN   # Beta release
#   X.YrcN  # Release Candidate
#
# Dev branch marker is: 'X.Y.dev' or 'X.Y.devN' where N is an integer.
# 'X.Y.dev0' is the canonical version of 'X.Y.dev'.
__version__ = '1.0.dev0'

# Add the project root to sys.path to allow import of 'problems'.
def _add_project_root_to_path():
    """
    Add the project root directory to sys.path.
    
    This allows importing 'problems' (from project root) when the package is
    installed in editable mode.
    """
    project_root = Path(__file__).resolve().parent.parent.parent
    project_root_str = str(project_root)
    
    if project_root_str not in sys.path:
        sys.path.insert(0, project_root_str)

# Allow to import 'problems' if 'optiprofiler' is imported.
_add_project_root_to_path()

# Public API of the optiprofiler package
from .modules import Feature, Problem, FeaturedProblem
from .profiles import benchmark
from .utils import show_versions

__all__ = [
    'Problem',
    'Feature',
    'FeaturedProblem',
    'benchmark',
    'show_versions',
]
