
import numpy as np
from matplotlib import pyplot as plt
import seaborn as sns
from datetime import datetime
import pandas as pd
from matplotlib.dates import DateFormatter
import dateutil.parser
import copy
import snowav.framework.figures


def basin_total(args, logger = None):
    '''
    Basin total SWE and SWI.

    Args
    ----------
    args : dict
        dictionary with required inputs, see swi() figure for more information.

    logger : dict
        snowav logger

    '''

    plotorder = args['plotorder']
    labels = args['labels']
    barcolors = args['barcolors']
    barcolors.insert(0,'black')
    swi_summary = args['swi_summary']
    swe_summary = args['swe_summary']
    swi_end_val = copy.deepcopy(swi_summary)
    swe_title = 'Basin SWE'
    swi_title = 'Basin SWI'

    sns.set_style('darkgrid')
    sns.set_context("notebook")

    plt.close(8)
    fig,(ax,ax1) = plt.subplots(num=8, figsize=args['figsize'],
                                dpi=args['dpi'], nrows = 1, ncols = 2)

    for iters,name in enumerate(plotorder):
        swe_summary[name].plot(ax=ax, color = barcolors[iters],
                               label = labels[name])
        swi_summary[name].plot(ax=ax1,color = barcolors[iters],
                               label='_nolegend_')

    if args['flt_flag']:
        for i,d in enumerate(args['flight_dates']):
            if i == 0:
                lb = 'flight update'.format(args['wy'])

            else:
                lb = '__nolabel__'

            ax.axvline(x=d, linestyle=':', linewidth=0.75, color='k', label=lb)

    x_end_date = args['end_date']

    # forecast
    # if args['forecast_flag']:
        # print('basin_total forecast in progress')
        # ax.axvline(x=snow.for_start_date,
        #            linestyle = ':',
        #            linewidth = 0.75,
        #            color = 'r')
        # ax1.axvline(x=snow.for_start_date,
        #            linestyle = ':',
        #            linewidth = 0.75,
        #            color = 'r')

    ax1.yaxis.set_label_position("right")
    ax1.set_xlim((datetime(args['wy']-1, 10, 1), x_end_date))
    ax.set_xlim((datetime(args['wy']-1, 10, 1), x_end_date))
    ax1.tick_params(axis='y')
    ax1.yaxis.tick_right()
    ax.legend(loc='upper left')

    # Put on the same yaxis
    swey = ax.get_ylim()
    swiy = ax1.get_ylim()

    if swey[1] < swiy[1]:
        ax1.set_ylim((-0.1,swiy[1]))
        ax.set_ylim((-0.1,swiy[1]))

    if swey[1] >= swiy[1]:
        ax1.set_ylim((-0.1,swey[1]))
        ax.set_ylim((-0.1,swey[1]))

    for tick,tick1 in zip(ax.get_xticklabels(),ax1.get_xticklabels()):
        tick.set_rotation(30)
        tick1.set_rotation(30)

    ax1.set_ylabel(r'[{}]'.format(args['vollbl']))
    ax1.set_xlabel('')
    ax.set_xlabel('')
    ax.axes.set_title(swe_title)
    ax1.axes.set_title(swi_title)
    ax.set_ylabel(r'[{}]'.format(args['vollbl']))

    del barcolors[0]

    fig_name_short = 'basin_total_'
    fig_name = '{}{}{}.png'.format(args['figs_path'],fig_name_short,args['directory'])
    if logger is not None:
        logger.info(' saving {}'.format(fig_name))
    snowav.framework.figures.save_fig(fig, fig_name)

    return fig_name_short
