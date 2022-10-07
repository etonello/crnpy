import numpy as np
import matplotlib.pyplot as plt


def weighted_quantile(values, quantiles, sample_weight=None,
                      values_sorted=False, old_style=False):
    ''' Very close to numpy.percentile, but supports weights.
        Note: quantiles should be in [0, 1]!
    :param values: numpy.array with data
    :param quantiles: array-like with many quantiles needed
    :param sample_weight: array-like of the same length as `array`
    :param values_sorted: bool, if True, then will avoid sorting of
        initial array
    :param old_style: if True, will correct output to be consistent
        with numpy.percentile.
    :return: numpy.array with computed quantiles.
    '''
    values = np.array(values)
    quantiles = np.array(quantiles)
    if sample_weight is None:
        sample_weight = np.ones(len(values))
    sample_weight = np.array(sample_weight)
    assert np.all(quantiles >= 0) and np.all(quantiles <= 1), \
        'quantiles should be in [0, 1]'

    if not values_sorted:
        sorter = np.argsort(values)
        values = values[sorter]
        sample_weight = sample_weight[sorter]

    weighted_quantiles = np.cumsum(sample_weight) - 0.5 * sample_weight
    if old_style:
        # To be convenient with numpy.percentile
        weighted_quantiles -= weighted_quantiles[0]
        weighted_quantiles /= weighted_quantiles[-1]
    else:
        weighted_quantiles /= np.sum(sample_weight)
    return np.interp(quantiles, weighted_quantiles, values)


def plot_samples(trajs, t, ref_traj=None, lower_q_bound=1/3, upper_q_bound=2/3, alpha=0.2, restraints=None, weights=None, crn=None, sim_incr=0.001):
    if weights is None:
        w = np.ones(trajs.shape[0])
    else:
        w = weights
    w /= np.sum(w)
    x = range(trajs.shape[1])
    qtrajs = np.apply_along_axis(lambda x: weighted_quantile(
        x, [lower_q_bound, 1/2, upper_q_bound], sample_weight=w), 0, trajs)
    mtrajs = np.sum(trajs * w[:, np.newaxis, np.newaxis], axis=0)
    qtrajs[0, :, :] = qtrajs[0, :, :] - qtrajs[1, :, :] + mtrajs
    qtrajs[2, :, :] = qtrajs[2, :, :] - qtrajs[1, :, :] + mtrajs
    qtrajs[1, :, :] = mtrajs
    fig, ax = plt.subplots(dpi=100)
    for i in range(trajs.shape[-1]):
        ax.plot(t, qtrajs[1, :, i], color=f'C{i}', label=crn.species[i])
        ax.fill_between(t, qtrajs[0, :, i],
                        qtrajs[2, :, i], color=f'C{i}', alpha=alpha)
        if ref_traj is not None:
            ax.plot(t, ref_traj[:, i], '--',
                    color=f'C{i}', label=crn.species[i])
    ax.set_xlabel('Time (s)')
    ax.set_ylabel('Species mass fraction')
    ax.set_ylim([0,1])
    handles, labels = ax.get_legend_handles_labels()
    unique = [(h, l) for i, (h, l) in enumerate(
        zip(handles, labels)) if l not in labels[:i]]
    ax.legend(*zip(*unique), loc='upper left', bbox_to_anchor=(1.1, 0.9))
    if restraints is not None:
        for r in restraints:
            ax.plot(r[2]*sim_incr, r[0], marker='o', color='k')

