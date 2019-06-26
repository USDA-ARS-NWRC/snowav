
import os
from datetime import datetime
import pandas as pd
import numpy as np
from matplotlib import pyplot as plt
import seaborn as sns
import copy
import matplotlib.dates as mdates
import snowav.framework.figures

def compare_runs(args, logger = None):
    '''
    Basin total values for existing database runs.

    Currently only plotting SWE and SWI.

    Args
    ----------
    args : dict
        dictionary with required inputs, see swi() figure for more information.

    logger : list
        snowav logger

    '''

    start_date = datetime(args['wy'] - 1,10,1)
    end_date = args['end_date']
    bid = args['plotorder'][0]

    sns.set_style('darkgrid')
    sns.set_context("notebook")

    for i, var in enumerate(args['variables']):

        plt.close(i)
        fig = plt.figure(num = i, figsize = args['figsize'], dpi = args['dpi'])
        ax = plt.gca()

        formatter = mdates.DateFormatter('%b')
        ax.xaxis.set_major_formatter(formatter)

        if 'z' in var:
            lbl = 'mm'
        if ('vol' in var) or 'avail' in var:
            lbl = args['vollbl']
        if var == 'density':
            lbl = 'kg/m^3'

        if var == 'swe_vol':
            title = 'Basin Total SWE'
        if var == 'swi_vol':
            title = 'Basin Total SWI'

        for name in args['dict'][var]:
            ax.plot(args['dict'][var][name], label = name)

        if args['flag']:
            for ix,d in enumerate(args['flight_dates']):
                if ix == 0:
                    lb = 'wy{} flight update'.format(args['wy'])
                else:
                    lb = '__nolabel__'

                ax.axvline(x = d, linestyle = ':', linewidth = 0.75,
                           color = 'k', label = lb)

        for tick in ax.get_xticklabels():
            tick.set_rotation(30)

        formatter = mdates.DateFormatter('%b')
        ax.xaxis.set_major_formatter(formatter)
        ax.legend()
        ax.set_ylabel('[{}]'.format(lbl))
        ax.set_xlim((datetime(args['wy']-1,10,1),datetime(args['wy'],8,1)))

        ax.set_title(title)
        plt.tight_layout()

        fig_name_short = 'compare_{}_'.format(var)
        fig_name = '{}{}{}.png'.format(args['figs_path'],fig_name_short,args['directory'])
        if logger is not None:
            logger.info(' saving {}'.format(fig_name))
        snowav.framework.figures.save_fig(fig, fig_name)

    return fig_name_short
