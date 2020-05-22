from datetime import datetime
from matplotlib import pyplot as plt
import matplotlib.dates as mdates
import os
import seaborn as sns

import snowav.framework.figures


def basin_total(plotorder, labels, end_date, barcolors, swi, swe, wy, vollbl,
                figsize, figs_path, fig_name, dpi=200, flight_flag=False,
                flight_dates=None, logger=None):
    """Basin total SWE and SWI.

    Args
    ------
    plotorder {list}: basins
    labels {list}: basin labels
    end_date {datetime}: end date
    barcolors {list}: colors
    swi {DataFrame}: daily swi
    swe {DataFrame}: daily swe
    wy {int}: water year
    vollbl {str}: volume label
    figsize {list}: figure dimensions
    figs_path {str}: base path for figures
    fig_name {str}: figure file name
    dpi {int}: figure dpi
    flight_flag {bool}: flights
    flight_dates {list}: flight dates
    """

    barcolors.insert(0, 'black')
    swe_title = 'Basin SWE'
    swi_title = 'Basin SWI'

    sns.set_style('darkgrid')
    sns.set_context("notebook")

    plt.close(8)
    fig, (ax, ax1) = plt.subplots(num=8, figsize=figsize, dpi=dpi, nrows=1,
                                  ncols=2)

    for iters, name in enumerate(plotorder):
        swe[name].plot(ax=ax, color=barcolors[iters], label=labels[name])
        swi[name].plot(ax=ax1, color=barcolors[iters], label='_nolegend_')

    del barcolors[0]

    if flight_flag:
        for i, d in enumerate(flight_dates):
            if i == 0:
                lb = 'flight update'.format(wy)
            else:
                lb = '__nolabel__'

            ax.axvline(x=d, linestyle=':', linewidth=0.75, color='k', label=lb)

    x_end_date = end_date
    ax1.yaxis.set_label_position("right")
    ax1.set_xlim((datetime(wy - 1, 10, 1), x_end_date))
    ax.set_xlim((datetime(wy - 1, 10, 1), x_end_date))
    ax1.tick_params(axis='y')
    ax1.yaxis.tick_right()
    ax.legend(loc='upper left', fontsize=8)
    swey = ax.get_ylim()
    swiy = ax1.get_ylim()

    if swey[1] < swiy[1]:
        ax1.set_ylim((-0.1, swiy[1]))
        ax.set_ylim((-0.1, swiy[1]))

    if swey[1] >= swiy[1]:
        ax1.set_ylim((-0.1, swey[1]))
        ax.set_ylim((-0.1, swey[1]))

    for tick, tick1 in zip(ax.get_xticklabels(), ax1.get_xticklabels()):
        tick.set_rotation(45)
        tick1.set_rotation(45)

    ax1.set_ylabel('[{}]'.format(vollbl))
    ax1.set_xlabel('')
    ax.set_xlabel('')
    ax.axes.set_title(swe_title)
    ax1.axes.set_title(swi_title)
    ax.set_ylabel('[{}]'.format(vollbl))

    ax.tick_params(axis='x', labelsize=8)
    ax1.tick_params(axis='x', labelsize=8)
    ax.xaxis.set_major_formatter(mdates.DateFormatter('%b-%d'))
    ax1.xaxis.set_major_formatter(mdates.DateFormatter('%b-%d'))

    fig_name = os.path.join(os.path.abspath(figs_path), fig_name)
    snowav.framework.figures.save_fig(fig, fig_name)

    if logger is not None:
        logger.info(' Saved: {}'.format(fig_name))

