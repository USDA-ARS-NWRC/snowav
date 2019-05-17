
import numpy as np
import matplotlib
from matplotlib import pyplot as plt
import matplotlib.colors as mcolors
from mpl_toolkits.axes_grid1 import make_axes_locatable
import seaborn as sns
import copy
import cmocean
import matplotlib.patches as mpatches
from datetime import datetime
import pandas as pd
from matplotlib.dates import DateFormatter
import pandas as pd
from snowav import database
from snowav.database.tables import Basins
import dateutil.parser

def basin_total(snow, forecast = None):

    '''

    '''
    name_append = snow.name_append
    swe_title = 'Basin SWE'
    swi_title = 'Basin SWI'

    sns.set_style('darkgrid')
    sns.set_context("notebook")

    plt.close(8)
    fig,(ax,ax1) = plt.subplots(num=8, figsize=snow.figsize,
                                dpi=snow.dpi, nrows = 1, ncols = 2)

    snow.barcolors.insert(0,'black')

    # Make df from database
    swe_summary = pd.DataFrame(columns = snow.plotorder)
    swi_summary = pd.DataFrame(columns = snow.plotorder)

    for bid in snow.plotorder:
        r = database.database.query(snow,
                                    datetime(snow.wy-1,10,1),
                                    snow.end_date,
                                    snow.run_name,
                                    bid,
                                    'swe_vol')

        r2 = database.database.query(snow,
                                    datetime(snow.wy-1,10,1),
                                    snow.end_date,
                                    snow.run_name,
                                    bid,
                                    'swi_vol')

        v = r[(r['elevation'] == 'total')]
        v2 = r2[(r2['elevation'] == 'total')]

        for iter,d in enumerate(v['date_time'].values):
            swe_summary.loc[d,bid] = v['value'].values[iter]
            swi_summary.loc[d,bid] = v2['value'].values[iter]

    swi_summary.sort_index(inplace=True)
    swe_summary.sort_index(inplace=True)
    swi_summary = swi_summary.cumsum()
    swi_end_val = copy.deepcopy(swi_summary)

    if snow.basin == 'LAKES' or snow.basin == 'RCEW':
        plotorder = [snow.plotorder[0]]
    else:
        plotorder = snow.plotorder

    for iters,name in enumerate(plotorder):
        swe_summary[name].plot(ax=ax, color = snow.barcolors[iters])
        swi_summary[name].plot(ax=ax1,color = snow.barcolors[iters], label='_nolegend_')

    if snow.flt_flag:
        for i,d in enumerate(snow.flight_diff_dates):
            if i == 0:
                lb = 'flight update'.format(snow.wy)
            else:
                lb = '__nolabel__'
            ax.axvline(x=d,linestyle=':',linewidth=0.75,color='k',label=lb)
            # ax1.axvline(x=d,linestyle=':',linewidth=0.75,color='k',label=lb)

    # add in other years
    x_end_date = snow.end_date

    # forecast
    if forecast is not None:
        start_date = snow.for_start_date
        end_date = snow.for_end_date
        run_name = snow.for_run_name
        name_append = snow.name_append + '_forecast'
        swe_title = 'Forecast Basin SWE'
        swi_title = 'Forecast Basin SWI'
        x_end_date = snow.for_end_date

        # Make df from database
        swe_summary = pd.DataFrame(columns = snow.plotorder)
        swi_summary = pd.DataFrame(columns = snow.plotorder)

        for bid in snow.plotorder:
            r = database.database.query(snow,
                                        start_date,
                                        end_date,
                                        run_name,
                                        bid,
                                        'swe_vol')

            r2 = database.database.query(snow,
                                        start_date,
                                        end_date,
                                        run_name,
                                        bid,
                                        'swi_vol')

            v = r[(r['elevation'] == 'total')]
            v2 = r2[(r2['elevation'] == 'total')]

            for iter,d in enumerate(v['date_time'].values):
                swe_summary.loc[d,bid] = v['value'].values[iter]
                swi_summary.loc[d,bid] = v2['value'].values[iter]

        swi_summary.sort_index(inplace=True)

        # as a starting spot, add actual run
        swi_summary.iloc[0,:] = swi_summary.iloc[0,:] + swi_end_val.iloc[-1,:].values
        swi_summary = swi_summary.cumsum()

        if snow.basin == 'LAKES' or snow.basin == 'RCEW':
            plotorder = [snow.plotorder[0]]
        else:
            plotorder = snow.plotorder

        for iters,name in enumerate(plotorder):
            swe_summary[name].plot(ax=ax,
                                   color = snow.barcolors[iters],
                                   linestyle = ':',
                                   label='_nolegend_')
            swi_summary[name].plot(ax=ax1,
                                   color = snow.barcolors[iters],
                                   linestyle = ':',
                                   label='_nolegend_')

        ax.axvline(x=snow.for_start_date,
                   linestyle = ':',
                   linewidth = 0.75,
                   color = 'r')
        ax1.axvline(x=snow.for_start_date,
                   linestyle = ':',
                   linewidth = 0.75,
                   color = 'r')

    ax1.yaxis.set_label_position("right")
    ax1.set_xlim((datetime(snow.wy -1 , 10, 1),x_end_date))
    ax.set_xlim((datetime(snow.wy - 1, 10, 1),x_end_date))
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

    ax1.set_ylabel(r'[{}]'.format(snow.vollbl))
    ax.axes.set_title(swe_title)
    ax1.axes.set_title(swi_title)
    ax.set_ylabel(r'[{}]'.format(snow.vollbl))

    # plt.tight_layout()
    del snow.barcolors[0]

    snow._logger.info(' saving {}basin_total_{}.png'.format(snow.figs_path,
                                                          name_append))
    plt.savefig('{}basin_total_{}.png'.format(snow.figs_path,name_append))
