import os
import sys


def add_optiprofiler():
    """
    Add the 'python' directory to the system path.
    """
    # Get the current directory of this file
    current_dir = os.path.dirname(os.path.abspath(__file__))

    # Get the root directory of the project
    root_dir = os.path.dirname(current_dir)

    # Construct the path to the 'python' directory
    python_dir = os.path.join(root_dir, 'python')

    # Check if the 'python' directory exists
    if python_dir not in sys.path:
        sys.path.insert(0, python_dir)
    return True