
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


def swe_volume(snow=None, forecast=None, day=None):
    '''

    '''

    if day is not None:
        outputs = day.outputs
        title = 'SWE \n {}'.format(day.nc_path)
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
        name_append = ''

        if day.basin == 'LAKES':
            imgx = day.imgx
            imgy = day.imgy
            dplcs = 2

    elif forecast is None:
        run_name = snow.run_name
        outputs = snow.outputs
        start_date = snow.start_date
        end_date = snow.end_date
        name_append = snow.name_append
        title = 'SWE \n {}'.format(snow.report_date.date().strftime("%Y-%-m-%-d"))

    else:
        run_name = snow.for_run_name
        outputs = snow.for_outputs
        start_date = snow.for_start_date
        end_date = snow.for_end_date
        name_append = snow.name_append + '_forecast'
        title = 'Forecast SWE \n {}'.format(end_date.date().strftime("%Y-%-m-%-d"))

    if snow is not None:
        depth_factor = snow.depth_factor
        masks = snow.masks
        plotorder = snow.plotorder
        figsize = snow.figsize
        dpi = snow.dpi
        depthlbl = snow.depthlbl
        vollbl = snow.vollbl
        elevlbl = snow.elevlbl
        basin = snow.basin
        dplcs = 1
        xlims = snow.xlims
        figs_path = snow.figs_path

        if snow.basin == 'LAKES':
            imgx = snow.imgx
            imgy = snow.imgy
            dplcs = 2

    state = copy.deepcopy(outputs['swe_z'][-1])
    state = np.multiply(state, depth_factor)

    qMin,qMax = np.nanpercentile(state,[0,99.8])
    clims = (qMin,qMax)
    pmask = masks[plotorder[0]]['mask']
    ixo = pmask == 0

    # Prepare no-snow and outside of the basin for the colormaps
    ixz = state == 0
    state[ixz] = -1

    # Colormap for snow.state
    colorsbad = plt.cm.Set2_r(np.linspace(0., 1, 1))
    colors1 = cmocean.cm.haline_r(np.linspace(0., 1, 254))
    colors = np.vstack((colorsbad,colors1))
    mymap = mcolors.LinearSegmentedColormap.from_list('my_colormap', colors)

    state[ixo] = np.nan
    mymap.set_bad('white',1.)

    sns.set_style('darkgrid')
    sns.set_context("notebook")

    plt.close(21)
    fig,(ax,ax1) = plt.subplots(num=21, figsize=figsize,
                                facecolor = 'white', dpi=dpi,
                                nrows = 1, ncols = 2)

    h = ax.imshow(state, clim=clims, cmap = mymap)

    if snow.basin == 'LAKES':
        ax.set_xlim(imgx)
        ax.set_ylim(imgy)

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

    h.axes.set_title(title)

    # Get basin-specific lims
    lims = plotlims(basin, plotorder)

    patches = [mpatches.Patch(color='grey', label='snow free')]

    # If there is meaningful snow-free area, include path and label
    if sum(sum(ixz)) > 1000:
        patches = [mpatches.Patch(color='grey', label='snow free')]
        ax.legend(handles=patches, bbox_to_anchor=(lims.pbbx, 0.05),
                  loc=2, borderaxespad=0. )

    # basin total and legend
    ax1.legend(loc=(lims.legx,lims.legy),markerscale = 0.5)

    # ax1.text(lims.btx,lims.bty,tlbl,horizontalalignment='center',
    #          transform=ax1.transAxes,fontsize = 10)

    if snow is not None:
        swe = pd.DataFrame(index = edges, columns = plotorder)

        for bid in snow.plotorder:
            r2 = database.database.query(snow, start_date, end_date,
                                        run_name, bid, 'swe_vol')

            for elev in snow.edges:
                v2 = r2[(r2['elevation'] == str(elev))
                      & (r2['date_time'] == end_date)
                      & (r2['basin_id'] == Basins.basins[bid]['basin_id'])]
                swe.loc[elev,bid] = v2['value'].values[0]

    else:
        swe =

    ix = len(barcolors)
    # Total basin label

    if len(plotorder) > 1:
        sumorder = plotorder[1::]
    else:
        sumorder = plotorder

    for iters,name in enumerate(sumorder):

        if snow.dplcs == 0:
            ukaf = str(np.int(np.nansum(swe[name])))
        else:
            ukaf = str(np.round(np.nansum(swe[name]),dplcs))

        if iters == 0:
            ax1.bar(range(0,len(edges)),swe[name],
                    color = barcolors[iters],
                    edgecolor = 'k',label = name + ': {} {}'.format(ukaf,vollbl))

        else:
            ax1.bar(range(0,len(snow.edges)),swe[name],
                    bottom = pd.DataFrame(swe[sumorder[0:iters]]).sum(axis = 1).values,
                    color = barcolors[iters], edgecolor = 'k',label = name + ': {} {}'.format(ukaf,vollbl))

    ylims = ax1.get_ylim()
    max = ylims[1] + ylims[1]*0.5
    min = 0
    ax1.set_ylim((min, max))

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

    ax1.set_xlabel('elevation [{}]'.format(elevlbl))
    ax1.set_ylabel('{} - per elevation band'.format(vollbl))
    ax1.yaxis.set_label_position("right")
    ax1.yaxis.tick_right()
    ax1.legend(loc = 2, fontsize = 10)

    for tick in ax1.get_xticklabels():
        tick.set_rotation(30)

    fig.tight_layout()
    fig.subplots_adjust(top=0.92,wspace = 0.2)

    if snow is not None:
        snow._logger.info('saving {}swe_volume_{}.png'.format(figs_path,name_append))
    plt.savefig('{}swe_volume_{}.png'.format(figs_path,name_append))
