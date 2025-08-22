import sys, os, io, re, importlib
from contextlib import redirect_stdout
import numpy as np
import pandas as pd

from ..utils import add_optiprofiler

add_optiprofiler()

from optiprofiler.problems import Problem

def s2mpj_load(problem_name, *args):
    """
    Load the S2MPJ problem.

    Parameters
    ----------
    problem_name : str
        Name of the problem in S2MPJ.
    *args : tuple
        Additional arguments to pass to the problem class constructor

    Returns
    -------
    `Problem`
        An instance of the Problem class.
    """

    # Add the path of the problem to the system path.
    try:
        current_dir = os.path.dirname(os.path.abspath(__file__))
        src_dir = os.path.join(current_dir, 'src')
        if src_dir not in sys.path:
            # Add the subfolder 'src' under the current directory to the system path (containing 's2mpjlib.py').
            sys.path.insert(0, src_dir)
        python_problems_dir = os.path.join(src_dir, 'python_problems')
        if python_problems_dir not in sys.path:
            # Add the subfolder 'python_problems' under the 'src' directory to the system path.
            sys.path.insert(0, python_problems_dir)
    except Exception as err:
        print(f"Failed to add the path of Python problems of S2MPJ to the system path: {err}")
        raise

    # Check if 'problem_name' has the pattern '_n_m' or 'n'. If it has, find the position of the pattern and return the dimension 'n' and the number of constraints 'm'.
    # Note that when 'm' is 0, the pattern is '_n' instead of '_n_0'.
    pattern = r'_(\d+)_(\d+)$|_(\d+)$'
    match = re.search(pattern, problem_name)
    if match:
        problem_name = problem_name[:match.start()]
        if match.group(1) and match.group(2):
            dim = int(match.group(1))
            mcon = int(match.group(2))
        elif match.group(3):
            dim = int(match.group(3))
            mcon = 0
        else:
            raise ValueError(f"Invalid problem name format: {problem_name}")
    else:
        dim = None
        mcon = None

    # Check csv file for the corresponding arguments.
    csv_file = os.path.join(current_dir, 'probinfo_python.csv')
    if match:
        df = pd.read_csv(csv_file, index_col=0)
        if problem_name in df.index:
            row = df.loc[problem_name]
            if dim is not None:
                dims_str = str(row['dims']).strip()
                mcons_str = str(row['mcons']).strip()
                argins_str = str(row['argins']).strip()
                
                if dims_str and mcons_str and argins_str:
                    dims_list = [int(x) for x in dims_str.split()]
                    mcons_list = [int(x) for x in mcons_str.split()]
                    
                    if problem_name == 'NUFFIELD' and argins_str.startswith('{5.0}{'):
                        variable_part = argins_str[len('{5.0}{'):-1]
                        variable_args = [float(x) for x in variable_part.split()]
                        fixed_args = [5.0]
                    elif problem_name == 'TRAINF' and argins_str.startswith('{1.5}{2}{'):
                        variable_part = argins_str[len('{1.5}{2}{'):-1]
                        variable_args = [float(x) for x in variable_part.split()]
                        fixed_args = [1.5, 2]
                    else:
                        variable_args = [float(x) for x in argins_str.split()]
                        fixed_args = []
                    
                    matching_indices = []
                    for i, (d, m) in enumerate(zip(dims_list, mcons_list)):
                        if d == dim and m == mcon:
                            matching_indices.append(i)
                    
                    if matching_indices:
                        index = matching_indices[0]
                        if index < len(variable_args):
                            selected_arg = variable_args[index]

                            if fixed_args:
                                args = args + tuple(fixed_args) + (selected_arg,)
                            else:
                                args = args + (selected_arg,)

    # Import the problem module dynamically.
    problem_module = importlib.import_module(f'python_problems.{problem_name}')
    # Get the class from the module.
    problem_class = getattr(problem_module, problem_name)
    # Create an instance of the class from the module.
    p = problem_class(*args)

    # List of feasibility problems.
    feasibility_list = [
        'AIRCRFTA', 'ARGAUSS', 'ARGLALE', 'ARGLBLE', 'ARGTRIG', 'ARTIF', 'BAmL1SP', 'BARDNE', 'BEALENE', 'BENNETT5', 'BIGGS6NE', 'BOOTH', 'BOXBOD', 'BRATU2D', 'BRATU2DT', 'BRATU3D', 'BROWNBSNE', 'BROWNDENE', 'BROYDN3D', 'CBRATU2D', 'CBRATU3D', 'CHANDHEQ', 'CHEMRCTA', 'CHWIRUT2', 'CLUSTER', 'COOLHANS', 'CUBENE', 'CYCLIC3', 'CYCLOOCF', 'CYCLOOCT', 'DANIWOOD', 'DANWOOD', 'DECONVBNE', 'DENSCHNBNE', 'DENSCHNDNE', 'DENSCHNFNE', 'DEVGLA1NE', 'DEVGLA2NE', 'DRCAVTY1', 'DRCAVTY2', 'DRCAVTY3', 'ECKERLE4', 'EGGCRATENE', 'EIGENA', 'EIGENB', 'ELATVIDUNE', 'ENGVAL2NE', 'ENSO', 'ERRINROSNE', 'ERRINRSMNE', 'EXP2NE', 'EXTROSNBNE', 'FLOSP2HH', 'FLOSP2HL', 'FLOSP2HM', 'FLOSP2TH', 'FLOSP2TL', 'FLOSP2TM', 'FREURONE', 'GENROSEBNE', 'GOTTFR', 'GROWTH', 'GULFNE', 'HAHN1', 'HATFLDANE', 'HATFLDBNE', 'HATFLDCNE', 'HATFLDDNE', 'HATFLDENE', 'HATFLDFLNE', 'HATFLDF', 'HATFLDG', 'HELIXNE', 'HIMMELBA', 'HIMMELBC', 'HIMMELBD', 'HIMMELBFNE', 'HS1NE', 'HS25NE', 'HS2NE', 'HYDCAR20', 'HYDCAR6', 'HYPCIR', 'INTEGREQ', 'INTEQNE', 'KOEBHELBNE', 'KOWOSBNE', 'KSS', 'LANCZOS1', 'LANCZOS2', 'LANCZOS3', 'LEVYMONE10', 'LEVYMONE5', 'LEVYMONE6', 'LEVYMONE7', 'LEVYMONE8', 'LEVYMONE9', 'LEVYMONE', 'LIARWHDNE', 'LIN', 'LINVERSENE', 'LSC1', 'LSC2', 'LUKSAN11', 'LUKSAN12', 'LUKSAN13', 'LUKSAN14', 'LUKSAN17', 'LUKSAN21', 'LUKSAN22', 'MANCINONE', 'METHANB8', 'METHANL8', 'MEYER3NE', 'MGH09', 'MGH10', 'MISRA1A', 'MISRA1B', 'MISRA1C', 'MISRA1D', 'MODBEALENE', 'MSQRTA', 'MSQRTB', 'MUONSINE', 'n10FOLDTR', 'NELSON', 'NONSCOMPNE', 'NYSTROM5', 'OSBORNE1', 'OSBORNE2', 'OSCIGRNE', 'OSCIPANE', 'PALMER1ANE', 'PALMER1BNE', 'PALMER1ENE', 'PALMER1NE', 'PALMER2ANE', 'PALMER2BNE', 'PALMER2ENE', 'PALMER3ANE', 'PALMER3BNE', 'PALMER3ENE', 'PALMER4ANE', 'PALMER4BNE', 'PALMER4ENE', 'PALMER5ANE', 'PALMER5BNE', 'PALMER5ENE', 'PALMER6ANE', 'PALMER6ENE', 'PALMER7ANE', 'PALMER7ENE', 'PALMER8ANE', 'PALMER8ENE', 'PENLT1NE', 'PENLT2NE', 'POROUS1', 'POROUS2', 'POWELLBS', 'POWELLSQ', 'POWERSUMNE', 'PRICE3NE', 'PRICE4NE', 'QINGNE', 'QR3D', 'RAT42', 'RAT43', 'RECIPE', 'REPEAT', 'RES', 'ROSZMAN1', 'RSNBRNE', 'SANTA', 'SEMICN2U', 'SEMICON1', 'SEMICON2', 'SPECANNE', 'SSBRYBNDNE', 'SSINE', 'THURBER', 'TQUARTICNE', 'VANDERM1', 'VANDERM2', 'VANDERM3', 'VANDERM4', 'VARDIMNE', 'VESUVIA', 'VESUVIO', 'VESUVIOU', 'VIBRBEAMNE', 'WATSONNE', 'WAYSEA1NE', 'WAYSEA2NE', 'YATP1CNE', 'YATP2CNE', 'YFITNE', 'ZANGWIL3'
    ]
    is_feasibility = problem_name in feasibility_list

    # Collect the problem parameters for transforming into a Problem instance.
    name = problem_name

    fun = lambda x: _getfun(p, is_feasibility, x)
    grad = lambda x: _getgrad(p, is_feasibility, x)
    hess = lambda x: _gethess(p, is_feasibility, x)

    x0 = p.x0
    xl = p.xlower
    xu = p.xupper

    if not hasattr(p, 'lincons'):
        p.lincons = None
    if not hasattr(p, 'nle'):
        p.nle = 0
    if not hasattr(p, 'neq'):
        p.neq = 0
    if not hasattr(p, 'nge'):
        p.nge = 0

    buf = io.StringIO()
    with redirect_stdout(buf):
        try:
            cx, jx = p.cJx(x0)[:2]
            if hasattr(cx, 'toarray'):
                cx = cx.toarray()
            if hasattr(jx, 'toarray'):
                jx = jx.toarray()
            bx = jx @ x0 - cx
        except Exception:
            jx = None
            bx = None

    nonlincons = np.setdiff1d(np.arange(p.m), p.lincons)
    idx_le = np.arange(p.nle)
    idx_eq = np.arange(p.nle, p.nle + p.neq)
    idx_ge = np.arange(p.nle + p.neq, p.nle + p.neq + p.nge)
    idx_aeq = np.intersect1d(idx_eq, p.lincons)
    idx_aub_le = np.intersect1d(idx_le, p.lincons)
    idx_aub_ge = np.intersect1d(idx_ge, p.lincons)
    idx_cle = np.intersect1d(idx_le, nonlincons)
    idx_ceq = np.intersect1d(idx_eq, nonlincons)
    idx_cge = np.intersect1d(idx_ge, nonlincons)
    aeq = jx[idx_aeq, :] if jx is not None else None
    aub = np.vstack([jx[idx_aub_le, :], -jx[idx_aub_ge, :]]) if jx is not None else None
    beq = bx[idx_aeq] if bx is not None else None
    bub = np.concatenate([bx[idx_aub_le], -bx[idx_aub_ge]]) if bx is not None else None

    getidx = lambda y, idx: y[idx] if y is not None else None
    ceq = lambda x: getidx(_getcx(p, x), idx_ceq)
    cub = lambda x: np.concatenate([getidx(_getcx(p, x), idx_cle),
                                    -getidx(_getcx(p, x), idx_cge)]) if getidx(_getcx(p, x), idx_cle) is not None else None

    getidx_list = lambda y, idx: [y[i] for i in idx] if y is not None else []
    hceq = lambda x: getidx_list(_getHx(p, x), idx_ceq)
    hcub = lambda x: getidx_list(_getHx(p, x), idx_cle) + getidx_list(_getHx(p, x), idx_cge)

    getidx_mat = lambda y, idx: y[idx, :] if y is not None else None
    jceq = lambda x: getidx_mat(_getJx(p, x), idx_ceq)
    jcub = lambda x: np.vstack([getidx_mat(_getJx(p, x), idx_cle),
                                -getidx_mat(_getJx(p, x), idx_cge)]) if getidx_mat(_getJx(p, x), idx_cle) is not None else None

    # Construct the Problem instance.
    problem = Problem(fun, x0, name=name, xl=xl, xu=xu, aub=aub, bub=bub, aeq=aeq, beq=beq, cub=cub, ceq=ceq, grad=grad, hess=hess, jcub=jcub, jceq=jceq, hcub=hcub, hceq=hceq)

    return problem

def s2mpj_select(options):
    """
    Select problems from the S2MPJ collection that satisfy given criteria.
    
    Parameters
    ----------
    options : dict
        A dictionary containing selection criteria:
        - ptype: problem type, string containing any of 'u', 'b', 'l', 'n' 
                 (default: 'ubln')
        - mindim: minimum dimension (default: 1)
        - maxdim: maximum dimension (default: inf)
        - minb: minimum number of bound constraints (default: 0)
        - maxb: maximum number of bound constraints (default: inf)
        - minlcon: minimum number of linear constraints (default: 0)
        - maxlcon: maximum number of linear constraints (default: inf)
        - minnlcon: minimum number of nonlinear constraints (default: 0)
        - maxnlcon: maximum number of nonlinear constraints (default: inf)
        - mincon: minimum total number of constraints (default: 0)
        - maxcon: maximum total number of constraints (default: inf)
        - oracle: oracle level, 0=zeroth-order, 1=first-order, 2=second-order 
                 (default: 0)
        - excludelist: list of problems to exclude (default: [])
    
    Returns
    -------
    list
        A list of problem names that satisfy the criteria.
    """
    # Set variable_size option
    variable_size = 'default'
    # Check if 'variable_size.txt' exists under the same directory as this script
    current_dir = os.path.dirname(os.path.abspath(__file__))
    variable_size_path = os.path.join(current_dir, 'variable_size.txt')
    if os.path.exists(variable_size_path):
        try:
            with open(variable_size_path, 'r') as f:
                variable_size = f.readline().strip()
        except:
            variable_size = 'default'
    
    if variable_size not in ['default', 'min', 'max', 'all']:
        raise ValueError("Invalid `variable_size` in the file `variable_size.txt`. Please set it to 'default', 'min', 'max', or 'all'.")
    
    # Initialize result lists
    problem_names = []
    
    # Set default values for options
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
        'oracle': 0,
        'excludelist': []
    }
    
    # Update default options with provided options
    for key in options:
        if key not in default_options:
            raise ValueError(f"Invalid option: {key}")
        default_options[key] = options[key]
    options = default_options
    
    # Add known problematic problems to exclude list
    exclude_problems = ['DANWOODLS', 'MISRA1CLS', 'ROSSIMP1', 'ROSSIMP2', 'ROSSIMP3']
    if options['oracle'] != 0:
        exclude_problems.append('NOZZLEfp')
    
    options['excludelist'] = list(set(options['excludelist'] + exclude_problems))
    
    # Get the directory where this script is located
    current_dir = os.path.dirname(os.path.abspath(__file__))
    
    # Load problem info from CSV file
    try:
        probinfo_path = os.path.join(current_dir, 'probinfo_python.csv')
        probinfo = pd.read_csv(probinfo_path)
    except:
        raise FileNotFoundError(f"Could not find or load problem info file at {probinfo_path}")
    
    # Helper function to safely convert values
    def safe_convert(value, default=0):
        if pd.isna(value) or value == 'unknown':
            return default
        try:
            return int(value)
        except (ValueError, TypeError):
            try:
                return float(value)
            except (ValueError, TypeError):
                return default
    
    # Process each problem
    for _, row in probinfo.iterrows():
        problem_name = row['problem_name']
        
        # Skip problems in exclude list
        if problem_name in options['excludelist']:
            continue
        
        # Check if problem type matches
        ptype = row['ptype']
        if not any(pt in options['ptype'] for pt in ptype):
            continue
        
        # Get problem attributes
        dim = safe_convert(row['dim'])
        mb = safe_convert(row['mb'])
        mlcon = safe_convert(row['mlcon'])
        mnlcon = safe_convert(row['mnlcon'])
        mcon = safe_convert(row['mcon'])
        
        # Check if argins (variable sizes) exist
        argins = row['argins'] if pd.notna(row['argins']) else ''
        dims = row['dims'].split() if pd.notna(row['dims']) else []
        mbs = row['mbs'].split() if pd.notna(row['mbs']) else []
        mlcons = row['mlcons'].split() if pd.notna(row['mlcons']) else []
        mnlcons = row['mnlcons'].split() if pd.notna(row['mnlcons']) else []
        mcons = row['mcons'].split() if pd.notna(row['mcons']) else []
        
        # Convert string arrays to numeric
        if dims:
            dims = [safe_convert(d) for d in dims]
        if mbs:
            mbs = [safe_convert(m) for m in mbs]
        if mlcons:
            mlcons = [safe_convert(m) for m in mlcons]
        if mnlcons:
            mnlcons = [safe_convert(m) for m in mnlcons]
        if mcons:
            mcons = [safe_convert(m) for m in mcons]
        
        # Check if default dimension and constraints satisfy criteria
        default_satisfy = (
            dim >= options['mindim'] and dim <= options['maxdim'] and
            mb >= options['minb'] and mb <= options['maxb'] and
            mlcon >= options['minlcon'] and mlcon <= options['maxlcon'] and
            mnlcon >= options['minnlcon'] and mnlcon <= options['maxnlcon'] and
            mcon >= options['mincon'] and mcon <= options['maxcon']
        )
        
        # If default satisfies and (no variable sizes or we only want default), add the problem
        if default_satisfy and (not dims or variable_size == 'default'):
            problem_names.append(problem_name)
        
        # Skip variable size processing if we only want default
        if variable_size == 'default':
            continue
        
        # Process variable sizes if they exist
        if dims:
            # Create mask for configurations that satisfy criteria
            mask = []
            for i in range(len(dims)):
                satisfies = (
                    dims[i] >= options['mindim'] and dims[i] <= options['maxdim'] and
                    (i >= len(mbs) or mbs[i] >= options['minb']) and
                    (i >= len(mbs) or mbs[i] <= options['maxb']) and
                    (i >= len(mlcons) or mlcons[i] >= options['minlcon']) and
                    (i >= len(mlcons) or mlcons[i] <= options['maxlcon']) and
                    (i >= len(mnlcons) or mnlcons[i] >= options['minnlcon']) and
                    (i >= len(mnlcons) or mnlcons[i] <= options['maxnlcon']) and
                    (i >= len(mcons) or mcons[i] >= options['mincon']) and
                    (i >= len(mcons) or mcons[i] <= options['maxcon']) and
                    (dims[i] != dim or (i < len(mcons) and mcons[i] != mcon))
                )
                mask.append(satisfies)
            
            idxs = [i for i, m in enumerate(mask) if m]
            
            if not idxs:
                if default_satisfy:
                    problem_names.append(problem_name)
                continue
            
            if variable_size == 'min':
                # Find minimum dimension among configurations that satisfy
                min_dim = min(dims[i] for i in idxs)
                idxs_dim = [i for i in idxs if dims[i] == min_dim]
                
                # Find minimum constraints among those
                if mcons and idxs_dim:
                    min_mcon = min(mcons[i] if i < len(mcons) else 0 for i in idxs_dim)
                    idxs = [i for i in idxs_dim if i >= len(mcons) or mcons[i] == min_mcon]
                    idxs = [idxs[0]]  # Take just the first
                    
                    # Compare with default
                    if default_satisfy and (dim < min_dim or (dim == min_dim and mcon < min_mcon)):
                        problem_names.append(problem_name)
                        continue
            
            elif variable_size == 'max':
                # Find maximum dimension among configurations that satisfy
                max_dim = max(dims[i] for i in idxs)
                idxs_dim = [i for i in idxs if dims[i] == max_dim]
                
                # Find maximum constraints among those
                if mcons and idxs_dim:
                    max_mcon = max(mcons[i] if i < len(mcons) else 0 for i in idxs_dim)
                    idxs = [i for i in idxs_dim if i < len(mcons) and mcons[i] == max_mcon]
                    idxs = [idxs[0]]  # Take just the first
                    
                    # Compare with default
                    if default_satisfy and (dim > max_dim or (dim == max_dim and mcon > max_mcon)):
                        problem_names.append(problem_name)
                        continue
            
            # For 'all' we just process all idxs
            
            # Add problems with specific dimensions/constraints
            for idx in idxs:
                if idx < len(mcons) and mcons[idx] > 0:
                    variant_name = f"{problem_name}_{dims[idx]}_{mcons[idx]}"
                else:
                    variant_name = f"{problem_name}_{dims[idx]}"
                
                problem_names.append(variant_name)
    
    return problem_names

def _getfun(p, is_feasibility, x):
    if is_feasibility:
        f = 0
    else:
        buf = io.StringIO()
        with redirect_stdout(buf):
            try:
                f = p.fx(x)
            except Exception:
                f = np.empty(0)
    return f

def _getgrad(p, is_feasibility, x):
    if is_feasibility:
        g = np.zeros_like(x)
    else:
        buf = io.StringIO()
        with redirect_stdout(buf):
            try:
                _, g = p.fgx(x)
                g = g.toarray() if hasattr(g, 'toarray') else g
            except Exception:
                g = None
    return g

def _gethess(p, is_feasibility, x):
    if is_feasibility:
        h = np.zeros((len(x), len(x)))
    else:
        buf = io.StringIO()
        with redirect_stdout(buf):
            try:
                _, _, h = p.fgHx(x)
                h = h.toarray() if hasattr(h, 'toarray') else h
            except Exception:
                h = None
    return h

def _getcx(p, x):
    buf = io.StringIO()
    with redirect_stdout(buf):
        try:
            c = p.cx(x)
            c = c.toarray() if hasattr(c, 'toarray') else c
        except Exception:
            c = None
    return c

def _getJx(p, x):
    buf = io.StringIO()
    with redirect_stdout(buf):
        try:
            _, j = p.cJx(x)[:2]
            j = j.toarray() if hasattr(j, 'toarray') else j
        except Exception:
            j = None
    return j

def _getHx(p, x):
    buf = io.StringIO()
    with redirect_stdout(buf):
        try:
            _, _, h = p.cJHx(x)
            for i in range(len(h)):
                if hasattr(h[i], 'toarray'):
                    h[i] = h[i].toarray()
        except Exception:
                h = None
        return h