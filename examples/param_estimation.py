from simulation import simulate

import numpy as np
from scipy.optimize import fmin

__author__ = "Elisa Tonello"
__copyright__ = "Copyright (c) 2016, Elisa Tonello"
__license__ = "BSD"
__version__ = "0.0.1"


def distance(x, y):
    """Eucleadian distance."""
    return np.linalg.norm(x - y)


def score_distance(sol, data, indices):
    """Return the total distance on the points identified by indices."""
    return sum([distance(sol[s][indices], data[s][indices]) for s in sol])


def estimate_params(crn, start_t, end_t, incr, sample_times, data, start_params, **kw):
    """Find the values of the parameters that minimize the distance
    to the points at sample_times. Use fmin from scipy."""
    # find indices for comparison
    t = np.arange(start_t, end_t, incr)
    indices = map(lambda x: np.where(t >= x)[0][0], sample_times)
    param_names = start_params.keys()
    param_values = start_params.values()

    # function to minimize
    def score(par):
        scores = []
        params = dict(zip(param_names, par))
        for dataset in data:
            initial_cond = dict((species, dataset[species][0]) for species in dataset)
            sol, _ = simulate(crn, params, initial_cond, start_t, end_t, incr)
            scores.append(score_distance(sol, dataset, indices))
        return sum(scores)

    return fmin(score, param_values, full_output = True, **kw)
