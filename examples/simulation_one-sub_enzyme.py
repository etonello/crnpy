import sys
import os
sys.path.insert(0, os.path.join(
    os.path.dirname(os.path.realpath(__file__)), '..'))

from crnpy.crn import CRN, from_react_strings

from simulation import simulate, plot_simulation


if __name__ == "__main__":
    # Simulation of one-substrate enzyme kinetics,
    # before and after removal of the intermediate via qss.

    crn = from_react_strings(['s + e (k_1)<->(k1) es', 'es ->(k2) p + e'])

    start_t, end_t, incr = 0, 1, 0.01
    initial = {'s': 5, 'e': 3, 'es': 0, 'p': 0}
    params = {'comp': 1, 'k1': 15, 'k_1': 10, 'k2': 6}

    colors = {'s': 'darkgreen', 'e': 'red', 'es': 'yellow', 'p': 'blue'}

    data, t = simulate(crn, params, initial, start_t, end_t, incr)
    title = "One-substrate enzyme kinetics\n $k_{1} = %s$, $k_{-1} = %s$, $k_{2} = %s$"%\
            (params['k1'], params['k_1'], params['k2'])
    plot_simulation("one-sub_enzyme_plt.png", data, t, colors, title)

    crn.qss('es')

    data, t = simulate(crn, params, initial, start_t, end_t, incr)
    title = "One-substrate enzyme kinetics, after qss reduction\n $k_{1} = %s$, $k_{-1} = %s$, $k_{2} = %s$"%\
            (params['k1'], params['k_1'], params['k2'])
    plot_simulation("one-sub_enzyme_reduced_plt.png", data, t, colors, title)
