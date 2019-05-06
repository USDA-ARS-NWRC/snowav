
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
from snowav.plotting.plotlims import plotlims as plotlims
from datetime import timedelta
import os
from snowav.utils.MidpointNormalize import MidpointNormalize



def swe_difference(day):
    '''

    '''

    outputs = day.outputs
    daily_outputs = day.daily_outputs
    results = day.results
    depth_factor = day.depth_factor
    masks = day.masks
    plotorder = day.plotorder
    dpi = 200
    figsize = (10,5)
    depthlbl = day.depthlbl
    vollbl = day.vollbl
    elevlbl = day.elevlbl
    basin = day.basin
    dplcs = 1
    xlims = day.xlims
    figs_path = day.figs_path
    name_append = 'day'
    edges = day.edges
    barcolors = day.barcolors

    if day.basin == 'LAKES':
        dplcs = 2

    if day.value in ['swe_z','swe_vol','swe_avail','swe_unavail']:
        v = 'swe_z'
    else:
        v = day.value

    title = '{}\n{} -\n{}'.format(v,day.nc_path[1],day.nc_path[0])

    swe_a = np.multiply(outputs[v][0], depth_factor)
    swe_b = np.multiply(outputs[v][1], depth_factor)
    delta_swe = swe_b - swe_a

    delta_swe_byelev = pd.DataFrame(index = edges, columns = plotorder)
    delta_swe_byelev = day.results[outputs['dates'][1]] - day.results[outputs['dates'][0]]
    delta_swe_byelev.fillna(value=0, inplace=True)
    delta_swe_byelev.drop(index='total', inplace=True)

    vMin,vMax = np.nanpercentile(delta_swe,[0.5,99.5])
    colorsbad = plt.cm.Set1_r(np.linspace(0., 1, 1))
    colors1 = cmocean.cm.matter_r(np.linspace(0., 1, 127))
    colors2 = plt.cm.Blues(np.linspace(0, 1, 128))
    colors = np.vstack((colorsbad,colors1, colors2))
    mymap = mcolors.LinearSegmentedColormap.from_list('my_colormap', colors)

    lims = plotlims(basin, plotorder)

    ixf = delta_swe == 0
    pmask = masks[plotorder[0]]['mask']
    ixo = pmask == 0

    delta_swe[ixf] = -100000
    delta_swe[ixo] = np.nan

    cmap = copy.copy(mymap)
    cmap.set_bad('white',1.)

    sns.set_style('darkgrid')
    sns.set_context("notebook")

    fig,(ax,ax1) = plt.subplots(num=0, figsize=figsize,
                                facecolor = 'white', dpi=dpi,
                                nrows = 1, ncols = 2)

    h = ax.imshow(delta_swe, interpolation='none',
        cmap = cmap, norm=MidpointNormalize(midpoint=0,
                                            vmin = vMin-0.01,vmax=vMax+0.01))

    # Basin boundaries
    for name in masks:
        ax.contour(masks[name]['mask'],cmap = "Greys",linewidths = 1)

    # Do pretty stuff for the left plot
    h.axes.get_xaxis().set_ticks([])
    h.axes.get_yaxis().set_ticks([])
    divider = make_axes_locatable(ax)
    cax = divider.append_axes("right", size="5%", pad=0.2)
    cbar = plt.colorbar(h, cax = cax)
    cbar.set_label('[{}]'.format(depthlbl))

    h.axes.set_title(title,fontsize=8)

    if dplcs == 0:
        tlbl = '{}: {} {}'.format(plotorder[0],
                             str(int(delta_swe_byelev[plotorder[0]].sum())),
                             vollbl)
    else:
        tlbl = '{}: {} {}'.format(plotorder[0],
                             str(np.round(delta_swe_byelev[plotorder[0]].sum(),
                                          dplcs)),vollbl)

    for iters,name in enumerate(lims.sumorder):

        if dplcs == 0:
            lbl = '{}: {} {}'.format(name,
                                str(int(delta_swe_byelev[name].sum())),
                                vollbl)
        else:
            lbl = '{}: {} {}'.format(name,
                                str(np.round(delta_swe_byelev[name].sum(),
                                dplcs)),vollbl)

        if iters == 0:
            ax1.bar(range(0,len(edges)),delta_swe_byelev[name],
                    color = barcolors[iters],
                    edgecolor = 'k',label = lbl)

        else:
            ax1.bar(range(0,len(edges)),delta_swe_byelev[name],
                    bottom = pd.DataFrame(delta_swe_byelev[lims.sumorder[0:iters]]).sum(axis = 1).values,
                    color = barcolors[iters], edgecolor = 'k',label = lbl)


    ax1.xaxis.set_ticks(range(0,len(edges)))
    plt.tight_layout()
    ax1.set_xlim((xlims[0]-0.5,xlims[1]+0.5))

    edges_lbl = []
    for i in range(0,len(edges)):
        edges_lbl.append(str(int(edges[int(i)])))

    ax1.set_xticklabels(str(i) for i in edges_lbl)
    for tick in ax1.get_xticklabels():
        tick.set_rotation(30)

    ax1.set_xlim((xlims[0]-0.5,xlims[1]+0.5))

    ylims = ax1.get_ylim()
    max = ylims[1] + abs(ylims[1] - ylims[0])
    min = ylims[0] - abs(ylims[1] - ylims[0])*0.1
    ax1.set_ylim((min, max))

    ax1.set_ylabel('{} - per elevation band'.format(vollbl))
    ax1.set_xlabel('elevation [{}]'.format(elevlbl))
    ax1.axes.set_title('Change in {}'.format(day.value))

    ax1.yaxis.set_label_position("right")
    ax1.yaxis.tick_right()

    # basin total and legend
    if len(plotorder) > 1:
        ax1.legend(loc=(lims.legx,lims.legy))

    ax1.text(lims.btx,lims.bty,tlbl,horizontalalignment='center',
             transform=ax1.transAxes,fontsize = 10)

    plt.tight_layout()
    plt.savefig('{}swe_difference_{}.png'.format(figs_path,name_append))

    if (day is not None) and day.show:
        plt.show()
