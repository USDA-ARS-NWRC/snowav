
import numpy as np
from matplotlib import pyplot as plt
import matplotlib.colors as mcolors
from mpl_toolkits.axes_grid1 import make_axes_locatable
import seaborn as sns
import copy
import cmocean
import matplotlib.patches as mpatches
import pandas as pd
from snowav import database
from snowav.database.tables import BASINS

def pixel_swe(snow):

    '''

    '''

    # Start and end indices
    ixs = np.where(snow.outputs['dates'] == snow.start_date)[0][0]
    ixe = np.where(snow.outputs['dates'] == snow.end_date)[0][0]

    # delta_swe = snow.outputs['swe_z'][ixe] - snow.outputs['swe_z'][ixs]
    # delta_swe = np.multiply(delta_swe,snow.depth_factor)

    # Make df from database
    swe_byelev = pd.DataFrame(index = snow.edges, columns = snow.plotorder)
    depth_byelev = pd.DataFrame(index = snow.edges, columns = snow.plotorder)

    for bid in snow.plotorder:
        r = database.database.query(snow.database,
                                                snow.end_date,
                                                snow.end_date,
                                                bid,
                                                'swe_z')
        r2 = database.database.query(snow.database,
                                                snow.end_date,
                                                snow.end_date,
                                                bid,
                                                'depth')

        for iter,elev in enumerate(snow.edges):
            v = r[r['elevation'] == str(elev)]
            v2 = r2[r2['elevation'] == str(elev)]
            swe_byelev.loc[elev,bid] = r['value'].values[iter]
            depth_byelev.loc[elev,bid] = r2['value'].values[iter]

    sns.set_style('darkgrid')
    sns.set_context('notebook')

    plt.close(2)
    fig,(ax,ax1) = plt.subplots(num=2, figsize=snow.figsize,
                                dpi=snow.dpi, nrows = 1, ncols = 2)

    sumorder = snow.plotorder[1:]
    if snow.basin == 'LAKES' or snow.basin == 'RCEW':
        sumorder = [snow.plotorder[0]]
        swid = 0.45
    else:
        sumorder = snow.plotorder[1::]
        swid = 0.25

    wid = np.linspace(-0.3,0.3,len(sumorder))

    for iters,name in enumerate(sumorder):
        ax.bar(range(0,len(snow.edges))-wid[iters],
                depth_byelev[name],
                color = snow.barcolors[iters], width = swid, edgecolor = 'k',label = name)

        ax1.bar(range(0,len(snow.edges))-wid[iters],
                swe_byelev[name],
                color = snow.barcolors[iters], width = swid, edgecolor = 'k',label = name)

    plt.tight_layout()
    xts = ax1.get_xticks()
    edges_lbl = []
    for i in xts[0:len(xts)-1]:
        edges_lbl.append(str(int(snow.edges[int(i)])))

    ax.set_xticklabels(str(i) for i in edges_lbl)
    ax1.set_xticklabels(str(i) for i in edges_lbl)
    for tick,tick1 in zip(ax.get_xticklabels(),ax1.get_xticklabels()):
        tick.set_rotation(30)
        tick1.set_rotation(30)

    ax.set_xlim((snow.xlims[0]-0.5,snow.xlims[1]+0.5))
    ax1.set_xlim((snow.xlims[0]-0.5,snow.xlims[1]+0.5))

    ylims = ax1.get_ylim()
    ax1.set_ylim((ylims[0],ylims[1]+ylims[1]*0.2))

    ylims = ax.get_ylim()
    ax.set_ylim((ylims[0],ylims[1]+ylims[1]*0.2))

    ax.set_ylabel('mean depth [%s]'%(snow.depthlbl))
    ax.set_xlabel('elevation [ft]')
    ax1.set_ylabel('mean SWE [%s]'%(snow.depthlbl))
    ax1.set_xlabel('elevation [ft]')
    ax.legend(loc='upper left')

    ax1.yaxis.set_label_position("right")
    ax1.yaxis.tick_right()

    ax.set_title('Mean Depth, %s'%(snow.end_date.date().strftime("%Y-%-m-%-d")))
    ax1.set_title('Mean SWE, %s'%(snow.end_date.date().strftime("%Y-%-m-%-d")))

    plt.tight_layout()
    snow._logger.info('saving figure to %smean_swe_depth%s.png'%(snow.figs_path,snow.name_append))
    plt.savefig('%smean_swe_depth%s.png'%(snow.figs_path,snow.name_append))
