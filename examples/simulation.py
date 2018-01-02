from crnpy.crn import CRN

import matplotlib.pyplot as plt
import numpy as np
from scipy import integrate
import sympy as sp

__author__ = "Elisa Tonello"
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
    plt.legend(data.keys(), loc = 'upper right', prop = {'size': 9})
    plt.title(title, fontsize = 12)
    plt.savefig(filename)
    plt.close()


def odes(x, t, species, eqs):
    """System of differential equations of the CRN."""
    vals = dict()
    for index, s in enumerate(species):
        vals[s] = x[index]
    vals[t] = t

    return [eq.evalf(subs = vals) for eq in eqs]


def simulate(crn, par, initial_cond, start_t, end_t, incr):
    """Simulate the deterministic dynamics."""
    # time
    times = np.arange(start_t, end_t, incr)

    # derivatives
    eqs = [e.subs(par.items()) for e in crn.equations()]

    # integration
    sol = integrate.odeint(lambda x, t: odes(x, t, [sp.Symbol(s) for s in crn.species], eqs),
                           [initial_cond[s] for s in crn.species],
                           times)
    return dict(zip(crn.species, np.transpose(sol))), times
