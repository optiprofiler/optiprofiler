def custom3():
    """
    This is a toy example to show how to construct a .py file that returns a dictionary that describes a linearly constrained optimization problem.
    """
    p_dict = {}
    p_dict['name'] = 'custom3'
    p_dict['x0'] = [1.0, 1.0]
    p_dict['fun'] = fun
    p_dict['xl'] = [-5.0, -5.0]
    p_dict['xu'] = [5.0, 5.0]
    p_dict['aub'] = [1.0, 1.0]
    p_dict['bub'] = 1
    p_dict['aeq'] = [1.0, -1.0]
    p_dict['beq'] = 0

    return p_dict


def fun(x):
    return sum(xi**2 for xi in x)