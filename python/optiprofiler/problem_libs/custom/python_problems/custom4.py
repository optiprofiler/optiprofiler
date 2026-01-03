import numpy as np

def custom4():
    """
    This is a toy example to show how to construct a .py file that returns a dictionary that describes an nonlinearly constrained optimization problem.
    """
    p_dict = {}
    p_dict['name'] = 'custom4'
    p_dict['x0'] = [1.0, 1.0]
    p_dict['fun'] = fun
    p_dict['xl'] = [-5.0, -5.0]
    p_dict['xu'] = [5.0, 5.0]
    p_dict['aub'] = [1.0, 1.0]
    p_dict['bub'] = 1
    p_dict['aeq'] = [1.0, -1.0]
    p_dict['beq'] = 0
    p_dict['cub'] = cub
    p_dict['ceq'] = ceq

    return p_dict

def fun(x):
    return sum(xi**2 for xi in x)

def cub(x):
    return -np.sin(x)

def ceq(x):
    return sum(xi**4 for xi in x)