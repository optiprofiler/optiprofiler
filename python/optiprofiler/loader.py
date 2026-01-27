import os
import re
import glob
import datetime
import numpy as np
import h5py
import pickle
from typing import Dict, List, Any, Tuple, Optional, Union, Callable
from enum import Enum
from pathlib import Path
from optiprofiler.utils import FeatureName, ProfileOption, FeatureOption, ProblemOption, ProblemError, get_logger

def search_in_dir(current_path: str, pattern: str, max_depth: int, current_depth: int, files: List[Dict]) -> List[Dict]:
    """Recursive function to search for files with a specified pattern within a limited directory depth."""
    # Check if the current depth exceeds the maximum depth
    if current_depth > max_depth:
        return files  # Stop searching deeper
    
    # Get files in the current directory matching the pattern
    matched_files = glob.glob(os.path.join(current_path, pattern))
    for file_path in matched_files:
        if os.path.isfile(file_path):
            file_info = {
                'name': os.path.basename(file_path),
                'folder': os.path.dirname(file_path),
                'date': datetime.datetime.fromtimestamp(os.path.getmtime(file_path)),
                'bytes': os.path.getsize(file_path),
                'isdir': False,
                'datenum': os.path.getmtime(file_path)
            }
            files.append(file_info)
    
    # Get all subdirectories in the current directory
    for item in os.listdir(current_path):
        item_path = os.path.join(current_path, item)
        # Skip '.' and '..' directories
        if os.path.isdir(item_path) and item != '.' and item != '..':
            # Recursively search in the subdirectory
            files = search_in_dir(item_path, pattern, max_depth, current_depth + 1, files)
    
    return files


def truncate_solvers(results_plib: Dict[str, Any], solvers_to_load: List[int]) -> Dict[str, Any]:
    """Truncate the loaded data by the 'solvers_to_load' field."""
    # Create a copy to avoid modifying the original
    results_plib = results_plib.copy()
    
    # Truncate solver-related fields
    results_plib['solver_names'] = [results_plib['solver_names'][i] for i in solvers_to_load]
    results_plib['fun_histories'] = results_plib['fun_histories'][:, solvers_to_load, :, :]
    results_plib['maxcv_histories'] = results_plib['maxcv_histories'][:, solvers_to_load, :, :]
    results_plib['fun_outs'] = results_plib['fun_outs'][:, solvers_to_load, :]
    results_plib['maxcv_outs'] = results_plib['maxcv_outs'][:, solvers_to_load, :]
    results_plib['n_evals'] = results_plib['n_evals'][:, solvers_to_load, :]
    results_plib['computation_times'] = results_plib['computation_times'][:, solvers_to_load, :]
    results_plib['solvers_successes'] = results_plib['solvers_successes'][:, solvers_to_load, :]
    results_plib['merit_histories'] = results_plib['merit_histories'][:, solvers_to_load, :, :]
    results_plib['merit_outs'] = results_plib['merit_outs'][:, solvers_to_load, :]
    
    return results_plib


def truncate_problems(results_plib: Dict[str, Any], problem_options: Dict[str, Any]) -> Dict[str, Any]:
    """Truncate the loaded data by the 'problem_options' field."""
    # Create a copy to avoid modifying the original
    results_plib = results_plib.copy()
    
    # Initialize the mask for problems to load
    p_to_load = np.ones(len(results_plib['problem_names']), dtype=bool)
    
    # Filter problems based on problem type
    if ProblemOption.PTYPE.value in problem_options:
        if hasattr(results_plib, 'ptype'):
            results_plib['ptype'] = list(set(results_plib['ptype']).intersection(problem_options[ProblemOption.PTYPE.value]))
        for i, ptype in enumerate(results_plib['problem_types']):
            if ptype not in problem_options[ProblemOption.PTYPE.value]:
                p_to_load[i] = False
    
    # Filter problems based on minimum dimension
    if ProblemOption.MINDIM.value in problem_options:
        if hasattr(results_plib, 'mindim'):
            results_plib['mindim'] = max(results_plib['mindim'], problem_options[ProblemOption.MINDIM.value])
        for i, dim in enumerate(results_plib['problem_dims']):
            if dim < problem_options[ProblemOption.MINDIM.value]:
                p_to_load[i] = False
    
    # Filter problems based on maximum dimension
    if ProblemOption.MAXDIM.value in problem_options:
        if hasattr(results_plib, 'maxdim'):
            results_plib['maxdim'] = min(results_plib['maxdim'], problem_options[ProblemOption.MAXDIM.value])
        for i, dim in enumerate(results_plib['problem_dims']):
            if dim > problem_options[ProblemOption.MAXDIM.value]:
                p_to_load[i] = False
    
    # Filter problems based on minimum number of bound constraints
    if ProblemOption.MINB.value in problem_options:
        if hasattr(results_plib, 'minb'):
            results_plib['minb'] = max(results_plib['minb'], problem_options[ProblemOption.MINB.value])
        for i, mb in enumerate(results_plib['problem_mbs']):
            if mb < problem_options[ProblemOption.MINB.value]:
                p_to_load[i] = False
    
    # Filter problems based on maximum number of bound constraints
    if ProblemOption.MAXB.value in problem_options:
        if hasattr(results_plib, 'maxb'):
            results_plib['maxb'] = min(results_plib['maxb'], problem_options[ProblemOption.MAXB.value])
        for i, mb in enumerate(results_plib['problem_mbs']):
            if mb > problem_options[ProblemOption.MAXB.value]:
                p_to_load[i] = False
    
    # Filter problems based on minimum number of linear constraints
    if ProblemOption.MINLCON.value in problem_options:
        if hasattr(results_plib, 'minlcon'):
            results_plib['minlcon'] = max(results_plib['minlcon'], problem_options[ProblemOption.MINLCON.value])
        for i, mlcon in enumerate(results_plib['problem_mlcons']):
            if mlcon < problem_options[ProblemOption.MINLCON.value]:
                p_to_load[i] = False
    
    # Filter problems based on maximum number of linear constraints
    if ProblemOption.MAXLCON.value in problem_options:
        if hasattr(results_plib, 'maxlcon'):
            results_plib['maxlcon'] = min(results_plib['maxlcon'], problem_options[ProblemOption.MAXLCON.value])
        for i, mlcon in enumerate(results_plib['problem_mlcons']):
            if mlcon > problem_options[ProblemOption.MAXLCON.value]:
                p_to_load[i] = False
    
    # Filter problems based on minimum number of nonlinear constraints
    if ProblemOption.MINNLCON.value in problem_options:
        if hasattr(results_plib, 'minnlcon'):
            results_plib['minnlcon'] = max(results_plib['minnlcon'], problem_options[ProblemOption.MINNLCON.value])
        for i, mnlcon in enumerate(results_plib['problem_mnlcons']):
            if mnlcon < problem_options[ProblemOption.MINNLCON.value]:
                p_to_load[i] = False
    
    # Filter problems based on maximum number of nonlinear constraints
    if ProblemOption.MAXNLCON.value in problem_options:
        if hasattr(results_plib, 'maxnlcon'):
            results_plib['maxnlcon'] = min(results_plib['maxnlcon'], problem_options[ProblemOption.MAXNLCON.value])
        for i, mnlcon in enumerate(results_plib['problem_mnlcons']):
            if mnlcon > problem_options[ProblemOption.MAXNLCON.value]:
                p_to_load[i] = False
    
    # Filter problems based on minimum number of constraints
    if ProblemOption.MINCON.value in problem_options:
        if hasattr(results_plib, 'mincon'):
            results_plib['mincon'] = max(results_plib['mincon'], problem_options[ProblemOption.MINCON.value])
        for i, mcon in enumerate(results_plib['problem_mcons']):
            if mcon < problem_options[ProblemOption.MINCON.value]:
                p_to_load[i] = False
    
    # Filter problems based on maximum number of constraints
    if ProblemOption.MAXCON.value in problem_options:
        if hasattr(results_plib, 'maxcon'):
            results_plib['maxcon'] = min(results_plib['maxcon'], problem_options[ProblemOption.MAXCON.value])
        for i, mcon in enumerate(results_plib['problem_mcons']):
            if mcon > problem_options[ProblemOption.MAXCON.value]:
                p_to_load[i] = False
    
    # Filter problems based on problem names
    if ProblemOption.PROBLEM_NAMES.value in problem_options:
        for i, name in enumerate(results_plib['problem_names']):
            if name not in problem_options[ProblemOption.PROBLEM_NAMES.value]:
                p_to_load[i] = False
    
    # Filter problems based on exclude list
    if ProblemOption.EXCLUDELIST.value in problem_options:
        for i, name in enumerate(results_plib['problem_names']):
            if name in problem_options[ProblemOption.EXCLUDELIST.value]:
                p_to_load[i] = False
    
    # Apply the mask to truncate the data
    p_to_load = np.where(p_to_load)[0]  # Convert boolean mask to indices
    
    # Truncate problem-related fields
    results_plib['fun_histories'] = results_plib['fun_histories'][p_to_load, :, :, :]
    results_plib['maxcv_histories'] = results_plib['maxcv_histories'][p_to_load, :, :, :]
    results_plib['fun_outs'] = results_plib['fun_outs'][p_to_load, :, :]
    results_plib['maxcv_outs'] = results_plib['maxcv_outs'][p_to_load, :, :]
    if results_plib['fun_inits'].ndim == 1:
        results_plib['fun_inits'] = results_plib['fun_inits'][p_to_load]
    else:
        results_plib['fun_inits'] = results_plib['fun_inits'][p_to_load, :]
    if results_plib['maxcv_inits'].ndim == 1:
        results_plib['maxcv_inits'] = results_plib['maxcv_inits'][p_to_load]
    else:
        results_plib['maxcv_inits'] = results_plib['maxcv_inits'][p_to_load, :]
    results_plib['n_evals'] = results_plib['n_evals'][p_to_load, :, :]
    results_plib['problem_names'] = [results_plib['problem_names'][i] for i in p_to_load]
    results_plib['problem_types'] = [results_plib['problem_types'][i] for i in p_to_load]
    results_plib['problem_dims'] = results_plib['problem_dims'][p_to_load]
    results_plib['problem_mbs'] = results_plib['problem_mbs'][p_to_load]
    results_plib['problem_mlcons'] = results_plib['problem_mlcons'][p_to_load]
    results_plib['problem_mnlcons'] = results_plib['problem_mnlcons'][p_to_load]
    results_plib['problem_mcons'] = results_plib['problem_mcons'][p_to_load]
    results_plib['computation_times'] = results_plib['computation_times'][p_to_load, :, :]
    results_plib['solvers_successes'] = results_plib['solvers_successes'][p_to_load, :, :]
    results_plib['merit_histories'] = results_plib['merit_histories'][p_to_load, :, :, :]
    results_plib['merit_outs'] = results_plib['merit_outs'][p_to_load, :, :]
    if results_plib['merit_inits'].ndim == 1:
        results_plib['merit_inits'] = results_plib['merit_inits'][p_to_load]
    else:
        results_plib['merit_inits'] = results_plib['merit_inits'][p_to_load, :]
    
    return results_plib


def merit_fun_compute(merit_fun: Callable, fun_values: np.ndarray, maxcv_values: np.ndarray, maxcv_inits: np.ndarray) -> np.ndarray:
    """Compute merit values from function values and constraint violations."""
    return merit_fun(fun_values, maxcv_values, maxcv_inits)


def recompute_merits(results_plib: Dict[str, Any], merit_fun: Callable) -> Dict[str, Any]:
    """Recompute the merit values for the loaded data if the 'merit_fun' is provided."""
    # Create a copy to avoid modifying the original
    results_plib = results_plib.copy()
    
    try:
        results_plib['merit_histories'] = merit_fun_compute(
            merit_fun, 
            results_plib['fun_histories'], 
            results_plib['maxcv_histories'], 
            results_plib['maxcv_inits']
        )
        results_plib['merit_outs'] = merit_fun_compute(
            merit_fun, 
            results_plib['fun_outs'], 
            results_plib['maxcv_outs'], 
            results_plib['maxcv_inits']
        )
        results_plib['merit_inits'] = merit_fun_compute(
            merit_fun, 
            results_plib['fun_inits'], 
            results_plib['maxcv_inits'], 
            results_plib['maxcv_inits']
        )
    except Exception as e:
        raise ValueError(f"The merit function provided in the options is not valid: {str(e)}")
    
    return results_plib


def load_results_from_h5(file_path: str) -> List[Dict[str, Any]]:
    """
    Load results from an HDF5 file.
    
    This function reads the HDF5 file created by `save_results_to_h5` and reconstructs the list of problem library results.
    It handles the decoding of strings (which are stored as bytes in HDF5) and the unpickling of complex objects.
    
    Parameters
    ----------
    file_path : str
        The path to the HDF5 file to load.
        
    Returns
    -------
    results_plibs : List[Dict[str, Any]]
        The loaded results.
    """
    results_plibs = []
    
    with h5py.File(file_path, 'r') as f:
        # Iterate through all problem libraries in the file
        for plib_key in f.keys():
            plib_group = f[plib_key]
            plib_data = {}
            
            # Load all datasets in this problem library
            for key in plib_group.keys():
                if isinstance(plib_group[key], h5py.Dataset):
                    # Check if it is a pickled dataset
                    if key.endswith('_pickled'):
                        data_void = plib_group[key][...]
                        data = pickle.loads(data_void.tobytes())
                        real_key = key[:-8]
                        plib_data[real_key] = data
                        continue

                    data = plib_group[key][...]
                    
                    # Convert string datasets back to lists or scalars
                    if data.dtype.kind == 'O' or data.dtype.kind == 'S' or data.dtype.kind == 'U':
                        if data.ndim == 1:
                            data = [s.decode('utf-8') if isinstance(s, bytes) else s for s in data]
                        elif data.ndim == 0:
                            s = data.item()
                            data = s.decode('utf-8') if isinstance(s, bytes) else s
                        
                    plib_data[key] = data
                elif isinstance(plib_group[key], h5py.Group):
                    # Handle nested groups
                    nested_data = {}
                    for nested_key in plib_group[key].keys():
                        if nested_key.endswith('_pickled'):
                            data_void = plib_group[key][nested_key][...]
                            data = pickle.loads(data_void.tobytes())
                            real_nested_key = nested_key[:-8]
                            nested_data[real_nested_key] = data
                            continue

                        nested_val = plib_group[key][nested_key][...]
                        
                        if nested_val.dtype.kind == 'O' or nested_val.dtype.kind == 'S' or nested_val.dtype.kind == 'U':
                            if nested_val.ndim == 1:
                                nested_val = [s.decode('utf-8') if isinstance(s, bytes) else s for s in nested_val]
                            elif nested_val.ndim == 0:
                                s = nested_val.item()
                                nested_val = s.decode('utf-8') if isinstance(s, bytes) else s
                        
                        nested_data[nested_key] = nested_val
                    plib_data[key] = nested_data
            
            results_plibs.append(plib_data)
    
    return results_plibs


def load_results(problem_options: Dict[str, Any], profile_options: Dict[str, Any]) -> Tuple[List[Dict[str, Any]], Dict[str, Any]]:
    """
    Load the results by the given options.
    
    This function restores the state of a previous experiment by loading the data from the 'data_for_loading.h5' file.
    It also filters the loaded data based on the current problem and profile options (e.g., selecting specific solvers or problems).
    
    Parameters
    ----------
    problem_options : Dict[str, Any]
        The problem options used to filter the loaded problems.
    profile_options : Dict[str, Any]
        The profile options used to specify which experiment to load (via 'load' option) and how to filter solvers.
        
    Returns
    -------
    results_plibs : List[Dict[str, Any]]
        The loaded and filtered results.
    profile_options : Dict[str, Any]
        The updated profile options (e.g., solver names might be updated based on loaded data).
    """
    results_plibs = []
    
    # Check if 'load' option is provided
    if ProfileOption.LOAD.value not in profile_options or not profile_options[ProfileOption.LOAD.value]:
        return results_plibs, profile_options
    
    # Check if 'load' is 'latest' or a timestamp
    time_stamp_pattern = r'^\d{4}(\d{2})(\d{2})_(\d{2})(\d{2})(\d{2})$'
    load_value = profile_options[ProfileOption.LOAD.value]
    
    if load_value != 'latest' and not re.match(time_stamp_pattern, load_value):
        raise ValueError("The option `load` should be either 'latest' or a time stamp in the format of 'yyyyMMdd_HHmmss'.")
    
    # Set the path to search for the data
    if ProfileOption.BENCHMARK_ID.value in profile_options:
        search_path = os.path.join(os.getcwd(), profile_options[ProfileOption.BENCHMARK_ID.value])
    else:
        search_path = os.getcwd()
        profile_options[ProfileOption.BENCHMARK_ID.value] = '.'
    
    # Find the path of the data to load
    time_stamp_files = []
    
    if load_value == 'latest':
        # Find all timestamp files and get the latest one
        time_stamp_files = search_in_dir(search_path, 'time_stamp_*.txt', 5, 0, [])
        
        if not time_stamp_files:
            raise FileNotFoundError(f"Failed to load data since no time_stamp files are found in the directory '{search_path}'. "
                                    f"Note that the search is limited to 5 levels of subdirectories.")
        
        # Sort files by date (newest first)
        def key_func(x):
            return x['date']
        time_stamp_files.sort(key=key_func, reverse=True)
        latest_time_stamp_file = time_stamp_files[0]
        path_data = latest_time_stamp_file['folder']
    else:
        # Search for the specific timestamp file
        pattern = f'time_stamp_{load_value}.txt'
        time_stamp_files = search_in_dir(search_path, pattern, 5, 0, [])
        
        if not time_stamp_files:
            raise FileNotFoundError(f"Failed to load data since no time_stamp named 'time_stamp_{load_value}.txt' is found in the directory '{search_path}'. "
                                    f"Note that the search is limited to 5 levels of subdirectories.")
        
        path_data = time_stamp_files[0]['folder']
    
    # Load data from the 'data_for_loading.h5' file
    data_file_path = os.path.join(path_data, 'data_for_loading.h5')
    
    try:
        results_plibs = load_results_from_h5(data_file_path)
    except Exception as e:
        raise FileNotFoundError(f"Failed to load data from the 'data_for_loading.h5' file in the directory '{path_data}': {str(e)}")
    
    if not results_plibs:
        raise ValueError(f"Failed to load any problem libraries from the 'data_for_loading.h5' file in the directory '{path_data}'.")
    
    # Process the 'solvers_to_load' field
    n_solvers_loaded = results_plibs[0]['fun_histories'].shape[1]
    
    if ProfileOption.SOLVERS_TO_LOAD.value in profile_options:
        solvers_to_load = profile_options[ProfileOption.SOLVERS_TO_LOAD.value]
        
        if any(idx >= n_solvers_loaded for idx in solvers_to_load) or (
                ProfileOption.SOLVER_NAMES.value in profile_options and 
                len(profile_options[ProfileOption.SOLVER_NAMES.value]) != len(solvers_to_load)):
            raise ValueError("The option `solvers_to_load` should be a list of different integers selected from 0 to the "
                             "total number of solvers minus 1 in the loaded data, and at least two indices should be provided.")
    else:
        solvers_to_load = list(range(n_solvers_loaded))
    
    # Assign values to 'profile_options' by the loaded 'results_plibs'
    if ProfileOption.SOLVER_NAMES.value not in profile_options:
        profile_options[ProfileOption.SOLVER_NAMES.value] = [results_plibs[0]['solver_names'][i] for i in solvers_to_load]
    
    if ProfileOption.FEATURE_STAMP.value not in profile_options and 'feature_stamp' in results_plibs[0]:
        profile_options[ProfileOption.FEATURE_STAMP.value] = results_plibs[0]['feature_stamp']
    
    # Select problem libraries to load by the 'plibs' field
    if ProblemOption.PLIBS.value in problem_options:
        valid_plib_names = [plib['plib'] for plib in results_plibs]
        plibs_to_load = problem_options[ProblemOption.PLIBS.value]
        
        if not all(isinstance(plib, str) for plib in plibs_to_load) or not all(plib in valid_plib_names for plib in plibs_to_load):
            raise ValueError(f"When you use the option `load`, the option `plibs` should be a list of strings, "
                             f"and each string should be one of the problem libraries in the loaded data ('{', '.join(valid_plib_names)}').")
        
        # Select the problem libraries to load
        results_plibs = [plib for plib in results_plibs if plib['plib'] in plibs_to_load]
        
        if not results_plibs:
            raise ValueError("No problem libraries are selected from the loaded data.")
    
    # Truncate data by the 'solvers_to_load' field
    for i, results_plib in enumerate(results_plibs):
        results_plibs[i] = truncate_solvers(results_plib, solvers_to_load)
        if 'results_plib_plain' in results_plib:
            results_plibs[i]['results_plib_plain'] = truncate_solvers(results_plib['results_plib_plain'], solvers_to_load)
    
    # Truncate data by 'problem_options'
    for i, results_plib in enumerate(results_plibs):
        results_plibs[i] = truncate_problems(results_plib, problem_options)
        if 'results_plib_plain' in results_plib:
            results_plibs[i]['results_plib_plain'] = truncate_problems(results_plib['results_plib_plain'], problem_options)
    
    # Recompute merit values if profile_options.merit_fun is provided
    if ProfileOption.MERIT_FUN.value in profile_options:
        merit_fun = profile_options[ProfileOption.MERIT_FUN.value]
        for i, results_plib in enumerate(results_plibs):
            results_plibs[i] = recompute_merits(results_plib, merit_fun)
            if 'results_plib_plain' in results_plib:
                results_plibs[i]['results_plib_plain'] = recompute_merits(results_plib['results_plib_plain'], merit_fun)
    
    return results_plibs, profile_options


def save_results_to_h5(results_plibs: List[Dict[str, Any]], file_path: str) -> None:
    """
    Save results to an HDF5 file.
    
    This function iterates through the list of problem library results and saves them into an HDF5 file.
    It handles different data types including numpy arrays, lists of strings, and nested dictionaries.
    For data types not natively supported by HDF5 (like lists of numbers or mixed types), it uses pickle serialization.
    
    Parameters
    ----------
    results_plibs : List[Dict[str, Any]]
        A list of dictionaries, where each dictionary contains the results for a problem library.
    file_path : str
        The path to the HDF5 file where the results will be saved.
    """
    with h5py.File(file_path, 'w') as f:
        # Save each problem library
        for i, plib in enumerate(results_plibs):
            plib_group = f.create_group(f'plib_{i}')
            
            # Save all fields in this problem library
            for key, value in plib.items():
                if isinstance(value, np.ndarray):
                    # Save numpy arrays directly
                    # Note: HDF5 does not support compression for scalar datasets (ndim=0).
                    if value.ndim == 0:
                        plib_group.create_dataset(key, data=value)
                    else:
                        plib_group.create_dataset(key, data=value, compression='gzip')
                elif isinstance(value, list):
                    # Convert lists to appropriate datasets
                    if all(isinstance(item, str) for item in value):
                        # List of strings
                        string_dt = h5py.special_dtype(vlen=str)
                        plib_group.create_dataset(key, data=np.array(value, dtype=string_dt))
                    else:
                        # Other lists (assuming numeric)
                        try:
                            arr = np.array(value)
                            # Handle scalar arrays resulting from single-element lists
                            if arr.ndim == 0:
                                plib_group.create_dataset(key, data=arr)
                            else:
                                plib_group.create_dataset(key, data=arr, compression='gzip')
                        except:
                            # If conversion fails (e.g., mixed types), store as pickled object
                            plib_group.create_dataset(f"{key}_pickled", data=np.void(pickle.dumps(value)))
                elif isinstance(value, dict):
                    # Handle nested dictionaries
                    nested_group = plib_group.create_group(key)
                    for nested_key, nested_value in value.items():
                        if isinstance(nested_value, np.ndarray):
                            if nested_value.ndim == 0:
                                nested_group.create_dataset(nested_key, data=nested_value)
                            else:
                                nested_group.create_dataset(nested_key, data=nested_value, compression='gzip')
                        elif isinstance(nested_value, list) and all(isinstance(item, str) for item in nested_value):
                            string_dt = h5py.special_dtype(vlen=str)
                            nested_group.create_dataset(nested_key, data=np.array(nested_value, dtype=string_dt))
                        else:
                            try:
                                arr = np.array(nested_value)
                                if arr.ndim == 0:
                                    nested_group.create_dataset(nested_key, data=arr)
                                else:
                                    nested_group.create_dataset(nested_key, data=arr, compression='gzip')
                            except:
                                nested_group.create_dataset(f"{nested_key}_pickled", data=np.void(pickle.dumps(nested_value)))
                else:
                    # Scalar values
                    if isinstance(value, str):
                        string_dt = h5py.special_dtype(vlen=str)
                        plib_group.create_dataset(key, data=np.array(value, dtype=string_dt))
                    else:
                        try:
                            plib_group.create_dataset(key, data=value)
                        except:
                            plib_group.create_dataset(f"{key}_pickled", data=np.void(pickle.dumps(value)))


def save_options(options: Dict[str, Any], file_path: Union[str, Path]) -> None:
    """
    Save the options dictionary to a pickle file.
    
    This is necessary because options dictionaries often contain Python objects 
    (like function handles or Enum members) that cannot be saved to HDF5 or JSON.
    
    Parameters
    ----------
    options : Dict[str, Any]
        The options dictionary to save.
    file_path : Union[str, Path]
        The path to the pickle file.
    """
    try:
        with open(file_path, 'wb') as f:
            pickle.dump(options, f)
    except Exception as e:
        get_logger().warning(f"Failed to save options to {file_path}: {e}")

