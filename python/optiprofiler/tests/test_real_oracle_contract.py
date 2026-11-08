"""Real optimization inputs must not be changed by complex-to-real casting."""
import numpy as np
import pytest
from optiprofiler import Problem


@pytest.mark.parametrize('channel', ['cub', 'ceq', 'grad', 'hess', 'jcub', 'jceq'])
def test_complex_oracle_arrays_are_rejected(channel):
    options = dict(fun=lambda x: float(x@x), x0=[1.],
                   cub=lambda x: np.array([1.]), ceq=lambda x: np.array([1.]))
    shape = (1, 1) if channel in ('hess', 'jcub', 'jceq') else (1,)
    # Constraint dimensions are discovered at x0. Keep that valid, then
    # exercise a later trial point outside the oracle's real-valued domain.
    options[channel] = lambda x: np.ones(shape) if x[0] == 1 else np.full(shape, 1+2j)
    problem = Problem(**options)
    with pytest.raises(ValueError, match='real'):
        getattr(problem, channel)([2.])


@pytest.mark.parametrize('channel', ['cub', 'ceq'])
def test_complex_initial_constraint_is_rejected_during_dimension_discovery(channel):
    with pytest.raises(ValueError):
        Problem(fun=lambda x: float(x@x), x0=[1.],
                **{channel: lambda x: np.array([1+2j])})


@pytest.mark.parametrize('method', ['fun', 'cub', 'ceq', 'grad', 'hess', 'maxcv'])
def test_complex_point_cannot_silently_move_to_real_part(method):
    problem = Problem(fun=lambda x: float(x@x), x0=[1.])
    with pytest.raises(ValueError, match='real'):
        getattr(problem, method)(np.array([1+2j]))


@pytest.mark.parametrize('key,value', [
    ('x0', np.array([1+2j])), ('xl', np.array([1+2j])),
    ('aub', np.array([[1+2j]])), ('bub', np.array([1+2j])),
])
def test_complex_problem_data_is_not_silently_cast(key, value):
    options = dict(fun=lambda x: float(x@x), x0=[1.], aub=[[1.]], bub=[2.])
    options[key] = value
    with pytest.raises(ValueError, match='real'):
        Problem(**options)
