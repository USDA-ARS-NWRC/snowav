
import numpy as np
import matplotlib
from matplotlib import pyplot as plt
import matplotlib.colors as mcolors
from mpl_toolkits.axes_grid1 import make_axes_locatable
import seaborn as sns
import copy
import cmocean
import matplotlib.patches as mpatches
from snowav import database
from snowav.database.tables import Basins
import pandas as pd

def state_by_elev(snow):
    '''
    Plots SWE by elevation.

    Currently no longer separating by melt/nonmelt, consider adding subplot
    if moving back to that scheme.

    '''

    # Fill df from database
    melt = pd.DataFrame(index = snow.edges, columns = snow.plotorder)
    nonmelt = pd.DataFrame(index = snow.edges, columns = snow.plotorder)

    for bid in snow.plotorder:
        r = database.database.query(snow, snow.start_date, snow.end_date,
                                    snow.run_name, bid, 'swe_avail')

        r2 = database.database.query(snow, snow.start_date, snow.end_date,
                                    snow.run_name, bid, 'swe_unavail')

        for elev in snow.edges:
            v = r[(r['elevation'] == str(elev))
                  & (r['date_time'] == snow.end_date)
                  & (r['basin_id'] == Basins.basins[bid]['basin_id'])]
            v2 = r2[(r2['elevation'] == str(elev))
                  & (r2['date_time'] == snow.end_date)
                  & (r2['basin_id'] == Basins.basins[bid]['basin_id'])]
            melt.loc[elev,bid] = v['value'].values[0]
            nonmelt.loc[elev,bid] = v2['value'].values[0]

    lim = np.max(melt[snow.plotorder[0]]) + np.max(nonmelt[snow.plotorder[0]])
    ylim = np.max(lim) + np.max(lim)*0.2
    # colors = ['xkcd:rose red','xkcd:cool blue']

    sns.set_style('darkgrid')
    sns.set_context("notebook")

    plt.close(7)
    fig = plt.figure(num=7, figsize=snow.figsize, dpi=snow.dpi)
    ax = plt.gca()

    for iters,name in enumerate(snow.plotorder):
        if iters == 0:
            clr = 'k'
        else:
            clr = snow.barcolors[iters - 1]

        if snow.dplcs == 0:
            kaf = str(np.int(sum(melt[name] + nonmelt[name])))
        else:
            kaf = str(np.round(sum(melt[name]+ nonmelt[name]),snow.dplcs))

        # Currently only plotting melt + nonmelt
        # If separating nonmelt, consider subplots
        ax.plot(range(0,len(snow.edges)),
                melt[name] + nonmelt[name],
                color = clr,
                label = name + ': {} {}'.format(kaf,snow.vollbl))

    xts = ax.get_xticks()
    if len(xts) < 6:
        dxt = xts[1] - xts[0]
        xts = np.arange(xts[0],xts[-1] + 1 ,dxt/2)
        ax.set_xticks(xts)

    edges_lbl = []
    for i in xts[0:len(xts)-1]:
        edges_lbl.append(str(int(snow.edges[int(i)])))

    ax.set_xticklabels(str(i) for i in edges_lbl)
    ax.set_xlabel('elevation [{}]'.format(snow.elevlbl))
    ax.legend(loc = 2, fontsize = 10)
    ax.set_ylim((0,ylim))
    ax.set_ylabel(snow.units)

    for tick in ax.get_xticklabels():
        tick.set_rotation(30)

    ax.set_xlim((snow.xlims[0]-0.5,snow.xlims[1]-1))
    fig.tight_layout()
    fig.subplots_adjust(top=0.92,wspace = 0.1)
    fig.suptitle('SWE, %s'%snow.report_date.date().strftime("%Y-%-m-%-d"))

    snow._logger.info('saving figure to {}swe_elev_{}.png'.format(snow.figs_path,snow.name_append))
    plt.savefig('{}swe_elev_{}.png'.format(snow.figs_path,snow.name_append))
