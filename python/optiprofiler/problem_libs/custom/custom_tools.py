import os
import sys
import importlib
import pandas as pd
import numpy as np

# Import Problem class from optiprofiler
from optiprofiler.opclasses import Problem


def custom_load(problem_name):
    """
    This is a toy example to show how to write a custom problem loader.
    """
    # Add the path 'python_problems' to sys.path
    current_dir = os.path.dirname(os.path.abspath(__file__))
    problems_dir = os.path.join(current_dir, 'python_problems')
    if problems_dir not in sys.path:
        sys.path.append(problems_dir)
    problem_module = importlib.import_module(f'{problem_name}')
    problem_handle = getattr(problem_module, problem_name)
    p_dict = problem_handle()
    return Problem(**p_dict)


def custom_select(options):
    """
    This is a toy example to show how to write a custom problem selector. We will use the csv file `probinfo_python.csv` created by `custom_tools.custom_get_info` to select the problems.
    """
    problem_names = []

    # Set default values for options.
    default_options = {
        'ptype': 'ubln',
        'mindim': 1,
        'maxdim': np.inf,
        'minb': 0,
        'maxb': np.inf,
        'minlcon': 0,
        'maxlcon': np.inf,
        'minnlcon': 0, 
        'maxnlcon': np.inf,
        'mincon': 0,
        'maxcon': np.inf,
        'excludelist': []
    }
    for key in default_options:
        options.setdefault(key, default_options[key])

    # Load problem info from CSV file and filter problems based on options.
    current_dir = os.path.dirname(os.path.abspath(__file__))
    probinfo_path = os.path.join(current_dir, 'probinfo_python.csv')
    probinfo = pd.read_csv(probinfo_path)
    for _, row in probinfo.iterrows():
        name = row['name']
        ptype = row['ptype']
        dim = row['dim']
        mb = row['mb']
        mlcon = row['mlcon']
        mnlcon = row['mnlcon']
        mcon = row['mcon']

        type_match = ptype in options['ptype']
        dim_match = (dim >= options['mindim'] and dim <= options['maxdim'])
        mb_match = (mb >= options['minb'] and mb <= options['maxb'])
        mlcon_match = (mlcon >= options['minlcon'] and mlcon <= options['maxlcon'])
        mnlcon_match = (mnlcon >= options['minnlcon'] and mnlcon <= options['maxnlcon'])
        mcon_match = (mcon >= options['mincon'] and mcon <= options['maxcon'])
        exclude_match = name not in options['excludelist']

        problem_match = type_match and dim_match and mb_match and mlcon_match and mnlcon_match and mcon_match and exclude_match

        if problem_match:
            problem_names.append(name)

    return problem_names


def custom_get_info():
    """
    This is a toy example to show one way of collecting information about the problems in the custom problem library. We will store the information in a CSV file so that `custom_select` can use it to select problems.
    """
    current_dir = os.path.dirname(os.path.abspath(__file__))
    problems_dir = os.path.join(current_dir, 'python_problems')
    # Collect all the py file names in the directory `python_problems`.
    problem_files = [f for f in os.listdir(problems_dir) if f.endswith('.py') and not f.startswith('__')]
    problem_names = [os.path.splitext(f)[0] for f in problem_files]
    problem_names = sorted(problem_names)  # Sort the file names.

    n_problems = len(problem_names)
    problem_info = {
        'name': problem_names,
        'ptype': ['u'] * n_problems,
        'dim': [2] * n_problems,
        'mb': [0] * n_problems,
        'mlcon': [0] * n_problems,
        'mnlcon': [0] * n_problems,
        'mcon': [0] * n_problems,
    }
    for i in range(n_problems):
        problem_name = problem_names[i]
        p = custom_load(problem_name)
        problem_info['name'][i] = p.name
        problem_info['ptype'][i] = p.ptype
        problem_info['dim'][i] = p.n
        problem_info['mb'][i] = p.mb
        problem_info['mlcon'][i] = p.mlcon
        problem_info['mnlcon'][i] = p.mnlcon
        problem_info['mcon'][i] = p.mcon

    df = pd.DataFrame(problem_info)
    df.to_csv(os.path.join(current_dir, 'probinfo_python.csv'), index=False)
    return df
