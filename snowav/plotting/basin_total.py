
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

def basin_total(snow):

    '''

    '''

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

    swi_summary = swi_summary.cumsum()

    sns.set_style('darkgrid')
    sns.set_context("notebook")

    plt.close(8)
    fig,(ax,ax1) = plt.subplots(num=8, figsize=snow.figsize,
                                dpi=snow.dpi, nrows = 1, ncols = 2)

    snow.barcolors.insert(0,'black')

    if snow.basin == 'LAKES' or snow.basin == 'RCEW':
        plotorder = [snow.plotorder[0]]
    else:
        plotorder = snow.plotorder

    for iters,name in enumerate(plotorder):
        swe_summary[name].plot(ax=ax, color = snow.barcolors[iters])
        swi_summary[name].plot(ax=ax1,color = snow.barcolors[iters], label='_nolegend_')

    if snow.flight_dates is not None:
        for d in snow.flight_dates:
            ax.axvline(x=d,linestyle = ':',linewidth = 0.75, color = 'k')
            ax1.axvline(x=d,linestyle = ':',linewidth = 0.75, color = 'k')

    ax1.yaxis.set_label_position("right")
    ax1.set_xlim((datetime(snow.wy -1 , 10, 1),snow.end_date))
    ax.set_xlim((datetime(snow.wy - 1, 10, 1),snow.end_date))
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

    ax1.set_ylabel(r'[%s]'%(snow.vollbl))
    ax.axes.set_title('Basin SWE')
    ax1.axes.set_title('Accumulated Basin SWI')
    ax.set_ylabel(r'[%s]'%(snow.vollbl))

    # plt.tight_layout()
    del snow.barcolors[0]

    snow._logger.info('saving figure to %sbasin_total_%s.png'%(snow.figs_path,
                                                  snow.name_append))
    plt.savefig('%sbasin_total_%s.png'%(snow.figs_path,snow.name_append))

    ########################################
    #         Second figure
    ########################################

    # This needs to be improved...
    if snow.basin_total_flag is False:
        snow._logger.info('No basin total summary files specified in config file, '
              + 'skipping...')
        return

    if snow.basin == 'BRB':
        main = 'Boise River Basin'
        multiswe = pd.read_csv(snow.summary_swe)
        multiswi = pd.read_csv(snow.summary_swi)
        multiswe['date'] = pd.to_datetime(multiswe['date'])
        multiswi['date'] = pd.to_datetime(multiswi['date'])

        multiswi.wy17 = multiswi.wy17.cumsum()
        multiswe.loc[range(0,len(swe_summary[main])),'wy18'] = swe_summary[main].values.copy()
        multiswi.loc[range(0,len(swi_summary[main])),'wy18'] = swi_summary[main].values.copy()

        if snow.units == 'SI':
            multiswe.wy17 = np.multiply(multiswe.wy17,0.00123348)
            multiswi.wy17 = np.multiply(multiswi.wy17,0.00123348)
            multiswe.wy16 = np.multiply(multiswe.wy16,0.00123348)
            multiswi.wy16 = np.multiply(multiswi.wy16,0.00123348)
            multiswe.wy15 = np.multiply(multiswe.wy15,0.00123348)
            multiswi.wy15 = np.multiply(multiswi.wy15,0.00123348)

        plt.close(8)
        fig,(ax,ax1) = plt.subplots(num=8, figsize=snow.figsize,
                                    dpi=snow.dpi, nrows = 1, ncols = 2)

        ax.plot(multiswe['date'], multiswe['wy13'], color = 'c',label = 'wy2013')
        ax.plot(multiswe['date'], multiswe['wy15'], color = 'g',label = 'wy2015')
        ax.plot(multiswe['date'], multiswe['wy16'], color = 'r',label = 'wy2016')
        ax.plot(multiswe['date'], multiswe['wy17'], color = 'k',label = 'wy2017')
        ax.plot(multiswe['date'], multiswe['wy18'], color = 'b', label = 'wy2018')
        #
        ax1.plot(multiswi['date'], multiswi['wy13'], color = 'c',label = 'wy2013')
        ax1.plot(multiswi['date'], multiswi['wy15'], color = 'g',label = 'wy2015')
        ax1.plot(multiswi['date'], multiswi['wy16'], color = 'r',label = 'wy2016')
        ax1.plot(multiswi['date'], multiswi['wy17'], color = 'k',label = 'wy2017')
        ax1.plot(multiswi['date'], multiswi['wy18'], color = 'b', label = 'wy2018')

        formatter = DateFormatter('%b')
        ax.xaxis.set_major_formatter(formatter)
        ax1.xaxis.set_major_formatter(formatter)

        ax1.yaxis.set_label_position("right")
        ax1.set_xlim((datetime(snow.wy-1, 10, 1),datetime(snow.wy, 8, 1)))
        ax.set_xlim((datetime(snow.wy-1, 10, 1),datetime(snow.wy, 8, 1)))
        ax1.tick_params(axis='y')
        ax1.yaxis.tick_right()
        ax.legend(loc='upper left')

        for tick,tick1 in zip(ax.get_xticklabels(),ax1.get_xticklabels()):
            tick.set_rotation(30)
            tick1.set_rotation(30)

        ax.set_ylabel('[%s]'%(snow.vollbl))
        ax1.set_ylabel('[%s]'%(snow.vollbl))
        ax.axes.set_title('Water Year SWE')
        ax1.axes.set_title('Accumulated Basin SWI')

        ax.set_ylim((-0.1,ax1.get_ylim()[1]))
        plt.tight_layout()

        snow._logger.info('saving figure to %sbasin_total_multiyr_%s.png'%(snow.figs_path,snow.name_append))
        plt.savefig('%sbasin_total_multiyr_%s.png'%(snow.figs_path,snow.name_append))
