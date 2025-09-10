def custom1():
    """
    This is a toy example to show how to construct a .py file that returns a
    dictionary that describes an unconstrained optimization problem.
    """
    p_dict = {}
    p_dict['name'] = 'custom1'
    p_dict['x0'] = [1.0, 1.0]
    p_dict['fun'] = fun

    return p_dict

def fun(x):
    return sum(xi**2 for xi in x)