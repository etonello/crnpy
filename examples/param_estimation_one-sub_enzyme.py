import sys
import os
sys.path.insert(0, os.path.join(
    os.path.dirname(os.path.realpath(__file__)), '..'))

from crnpy.crn import CRN, from_react_strings

from simulation import simulate
from param_estimation import score_distance, estimate_params

import matplotlib.pyplot as plt
import numpy as np


def create_plot(ni, position, t, data, filename, add_times = None, add_data = None, title = None, ylabel = None):
    colors = {'s': 'darkgreen', 'e': 'red', 'es': 'yellow', 'p': 'blue'}
    ax = plt.subplot(ni, 3, position)
    ax.spines["right"].set_visible(False)
    ax.spines["top"].set_visible(False)
    ax.yaxis.set_ticks_position('left')
    ax.xaxis.set_ticks_position('bottom')
    if ylabel: ax.set_ylabel(ylabel, fontsize = 6)
    plt.tick_params(axis = 'both', which = 'major', labelsize = 4)
    map(lambda species: plt.plot(t, data[species], color = colors[species]), data.keys())
    if add_times is not None: map(lambda species: plt.plot(add_times, add_data[species], 'D', markersize = 2, color = colors[species]), data.keys())
    plt.legend(data.keys(), loc = 'upper right', prop = {'size': 6})
    if title: plt.title(title, fontsize = 6)
    plt.savefig(filename)


def plot_comparison():
    ni = len(initials)
    plt.figure()
    filename = 'one-sub_enzyme_comparison.png'
    qss_scores = [score_distance(qss_sols[i], sols[i], indices) for i in range(ni)]
    estimated_scores = [score_distance(estimate_sols[i], sols[i], indices) for i in range(ni)]
    print("Qss scores: {}".format(qss_scores))
    print("Qss form, estimated parameters scores:".format(estimated_scores))

    for i in range(ni):
        ylabel = '$s_0 = %s$, $e_0 = %s$'%(initials[i]['s'], initials[i]['e'])
        title = "Original model\n$k_{1} = %s$, $k_{-1} = %s$, $k_{2} = %s$"%\
                (round(params['k1']), round(params['k_1']), round(params['k2']))
        if i == 0: create_plot(ni, 1 + i * 3, times, sols[i], filename, title = title, ylabel = ylabel)
        else: create_plot(ni, 1 + i * 3, times, sols[i], filename, ylabel = ylabel)

        title = "Qss model\ntotal err = $%s$"%(round(sum(qss_scores)))
        sample_data = [dict((s, sol[s][indices]) for s in sol) for sol in sols]
        if i == 0: create_plot(ni, 2 + i * 3, times, qss_sols[i], filename, add_times = sample_times, add_data = sample_data[i], title = title)
        else: create_plot(ni, 2 + i * 3, times, qss_sols[i], filename, add_times = sample_times, add_data = sample_data[i])

        title = "Estimated params model, total err = $%s$\n$k_{1} = %s$, $k_{-1} = %s$, $k_{2} = %s$"%\
                (round(sum(estimated_scores)), round(estimated_params['k1']), \
                round(estimated_params['k_1']), round(estimated_params['k2']))
        if i == 0: create_plot(ni, 3 + i * 3, times, estimate_sols[i], filename, add_times = sample_times, add_data = sample_data[i], title = title)
        else: create_plot(ni, 3 + i * 3, times, estimate_sols[i], filename, add_times = sample_times, add_data = sample_data[i])

    plt.close()


if __name__ == "__main__":
    # Given the structure and form of the rates of the network
    # identified by the qss reduction, what are the best values
    # for the parameters, to approximate the original network?
    # Wse the initial network to create "experimental data",
    # and find the parameters that minimize the distance.

    crn = from_react_strings(['s + e (k_1)<->(k1) es', 'es ->(k2) p + e'])

    start_t, end_t, incr = 0, 1, 0.01
    params = {'k1': 15, 'k_1': 10, 'k2': 6}
    initials = ({'s': 5, 'e': 3, 'es': 0, 'p': 0},
                {'s': 2, 'e': 3, 'es': 0, 'p': 0},
                {'s': 5, 'e': 1, 'es': 0, 'p': 0})

    # times where solution is sampled
    sample_times = np.arange(0, 1, 0.2)
    times = np.arange(start_t, end_t, incr)
    indices = map(lambda x: np.where(times >= x)[0][0], sample_times)

    # similate the original dynamics
    sols = [simulate(crn, params, initial, start_t, end_t, incr)[0] for initial in initials]

    # reduce network and simulate the dynamics
    crn.qss('es')
    qss_sols = [simulate(crn, params, initial, start_t, end_t, incr)[0] for initial in initials]

    # find the parameters that better approximate the original dynamics in the sample_times
    # using initial parameters as initial guesses
    # then simulate the system with the new parameters
    estimate = estimate_params(crn, start_t, end_t, incr, sample_times, sols, params, maxiter = 1000, disp = True)
    estimated_params = dict(zip(params.keys(), estimate[0]))
    estimate_sols = [simulate(crn, estimated_params, initial, start_t, end_t, incr)[0] for initial in initials]

    # plot the results
    plot_comparison()
