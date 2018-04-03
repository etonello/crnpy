from crnpy.crn import CRN

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import numpy as np
from scipy import integrate
import sympy as sp

__author__ = "Elisa Tonello, Eitan Lees"
__copyright__ = "Copyright (c) 2016, Elisa Tonello"
__license__ = "BSD"
__version__ = "0.0.1"


def plot_simulation(filename, data, t, colors, title):
    """Create a png file with the plot of the simulation."""
    plt.figure()
    ax = plt.subplot(1, 1, 1)
    ax.spines["right"].set_visible(False)
    ax.spines["top"].set_visible(False)
    ax.yaxis.set_ticks_position('left')
    ax.xaxis.set_ticks_position('bottom')
    plt.tick_params(axis = 'both', which = 'major', labelsize = 8)
    for species in data.keys():
        plt.plot(t, data[species], color = colors[species])
    plt.legend(list(data.keys()), loc = 'upper right', prop = {'size': 9})
    plt.title(title, fontsize = 12)
    plt.savefig(filename)
    plt.close()


def simulate(crn, par, initial_cond, start_t, end_t, incr):
    """Simulate the deterministic dynamics."""
    # time
    times = np.arange(start_t, end_t, incr)

    # inserting rate constants in derivatives
    eqs = [e.subs(par.items()) for e in crn.equations()]

    # turning sympy equations into lambda functions
    lam_funcs = list(map(lambda eq: sp.lambdify(crn.species, eq), eqs))

    # integration
    sol = integrate.odeint(lambda x, t: list(map(lambda func: func(*x), lam_funcs)),
                           [initial_cond[s] for s in crn.species],
                           times)
    return dict(zip(crn.species, np.transpose(sol))), times
