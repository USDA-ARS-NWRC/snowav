
import numpy as np
from matplotlib import pyplot as plt
import matplotlib.colors as mcolors
from mpl_toolkits.axes_grid1 import make_axes_locatable
import seaborn as sns
import copy
import cmocean
import matplotlib.patches as mpatches
from snowav import database
from snowav.database.tables import BASINS
import pandas as pd

def state_by_elev(snow):
    '''
    Plots SWE by elevation, delineated by melt/nonmelt

    '''

    # snow.melt
    # snow.nonmelt
    # Calculate accumulated swi during the specified period
    ixs = np.where(snow.outputs['dates'] == snow.start_date)[0][0]
    ixe = np.where(snow.outputs['dates'] == snow.end_date)[0][0]
    accum = np.zeros((snow.nrows,snow.ncols))

    for n in range(ixs,ixe):
        accum = accum + snow.outputs['swi_z'][n]

    accum = np.multiply(accum,snow.depth_factor)

    # Make df from database
    snow.melt = pd.DataFrame(index = snow.edges, columns = snow.plotorder)
    snow.nonmelt = pd.DataFrame(index = snow.edges, columns = snow.plotorder)

    loc = snow.database
    start_date = snow.end_date
    end_date = snow.end_date
    bid = snow.plotorder

    for bid in snow.plotorder:
        r = database.database.query_basin_value(loc, start_date,
                                                end_date, bid, 'swe_avail')
        r2 = database.database.query_basin_value(loc, start_date,
                                                end_date, bid, 'swe_unavail')

        for elev in snow.edges:
            v = r[r['elevation'] == str(elev)]
            v2 = r2[r2['elevation'] == str(elev)]
            snow.melt.loc[elev,bid] = np.nansum(v['value'].values)
            snow.nonmelt.loc[elev,bid] = np.nansum(v2['value'].values)

    lim = np.max(snow.melt[snow.plotorder[0]]) + np.max(snow.nonmelt[snow.plotorder[0]])
    ylim = np.max(lim) + np.max(lim)*0.3
    colors = ['xkcd:rose red','xkcd:cool blue']
    fs = list(snow.figsize)
    fs[0] = fs[0]*0.9
    fs = tuple(fs)

    sns.set_style('darkgrid')
    sns.set_context("notebook")

    nf = len(snow.masks)

    if nf > 1 and nf < 5:
        nr = 2
        nc = 2
    elif nf < 7:
        nr = 2
        nc = 3
    if nf == 1:
        nr = 1
        nc = 1

    plt.close(7)
    fig,ax = plt.subplots(num=7, figsize=fs, dpi=snow.dpi,
                           nrows = nr, ncols = nc)
    ax = np.array(ax)
    axs = ax.ravel()

    if nf == 5:
        fig.delaxes(axs[5])

    for iters,name in enumerate(snow.plotorder):
        # iters = 0
        # name = snow.plotorder[iters]
        axs[iters].bar(range(0,len(snow.edges)),snow.melt[name],
                       color = colors[0], bottom = snow.nonmelt[name])
        axs[iters].bar(range(0,len(snow.edges)),snow.nonmelt[name],
                       color = colors[1], label = 'unavail ')

        xts = axs[iters].get_xticks()
        if len(xts) < 6:
            dxt = xts[1] - xts[0]
            xts = np.arange(xts[0],xts[-1] + 1 ,dxt/2)
            axs[iters].set_xticks(xts)

        edges_lbl = []
        for i in xts[0:len(xts)-1]:
            edges_lbl.append(str(int(snow.edges[int(i)])))

        axs[iters].set_xticklabels(str(i) for i in edges_lbl)

        if iters > nc - 1:
            axs[iters].set_xlabel('elevation [%s]'%(snow.elevlbl))

        # Put yaxis on right
        if iters == nc - 1 or iters == nf -1 :
            axs[iters].yaxis.set_label_position("right")
            axs[iters].yaxis.tick_right()

        if iters == 1 and nc == 3:
            axs[iters].set_yticklabels([])

        if iters <= nc - 1 and nf != 1:
            axs[iters].set_xticklabels([])

        axs[iters].tick_params(axis='x')
        axs[iters].tick_params(axis='y')

        # Get basin total storage in strings for label
        kaf = str(np.int(sum(snow.melt[name]) + sum(snow.nonmelt[name])))

        if snow.dplcs == 0:
            kaf = str(np.int(sum(snow.melt[name]) + sum(snow.nonmelt[name])))
        else:
            kaf = str(np.round(sum(snow.melt[name])
                               + sum(snow.nonmelt[name]),snow.dplcs))

        axs[iters].text(0.5,0.92,'%s - %s %s'
                        %(snow.masks[name]['label'],kaf,snow.vollbl),
                        horizontalalignment='center',
                        transform=axs[iters].transAxes, fontsize = 10)

        if iters == 1 and nc == 3:
            axs[iters].set_yticklabels([])
        else:
            axs[iters].set_ylabel(snow.units)

        lbl = []
        for n in (0,1):
            if n == 0:
                if snow.dplcs == 0:
                    kafa = str(np.int(sum(snow.melt[name])))
                else:
                    kafa = str(np.round(sum(snow.melt[name]),snow.dplcs))

                tmpa = (r'avail = %s')%(kafa)
                lbl.append(tmpa)

            if snow.dplcs == 0:
                kafna = str(np.int(sum(snow.nonmelt[name])))
            else:
                kafna = str(np.round(sum(snow.nonmelt[name]),snow.dplcs))

            tmpna = ('unavail = %s')%(kafna)
            lbl.append(tmpna)

        axs[iters].legend(lbl, loc = (0.025, 0.65),fontsize = 9)
        axs[iters].set_ylim((0,ylim))

        for tick in axs[iters].get_xticklabels():
            tick.set_rotation(30)

        axs[iters].set_xlim((snow.xlims[0]-0.5,snow.xlims[1]+0.5))

    fig.tight_layout()
    # for n in range(0,len(axs)):
    #     axs[n].set_xticks(xts)
    #     axs[n].set_xlim(snow.xlims)

    fig.subplots_adjust(top=0.92,wspace = 0.1)
    fig.suptitle('SWE, %s'%snow.end_date.date().strftime("%Y-%-m-%-d"))

    snow._logger.info('saving figure to %sswe_elev%s.png'%(snow.figs_path,snow.name_append))
    plt.savefig('%sswe_elev%s.png'%(snow.figs_path,snow.name_append))
