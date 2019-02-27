
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


def swe_volume(snow):
    '''

    '''

    state = copy.deepcopy(snow.outputs['swe_z'][snow.ixe])
    state = np.multiply(state,snow.depth_factor)

    qMin,qMax = np.nanpercentile(state,[0,99.8])
    clims = (qMin,qMax)
    pmask = snow.masks[snow.plotorder[0]]['mask']
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
    fig,(ax,ax1) = plt.subplots(num=21, figsize=snow.figsize,
                                facecolor = 'white', dpi=snow.dpi,
                                nrows = 1, ncols = 2)

    h = ax.imshow(state, clim=clims, cmap = mymap)

    if snow.basin == 'LAKES':
        ax.set_xlim(snow.imgx)
        ax.set_ylim(snow.imgy)

    # Basin boundaries
    for name in snow.masks:
        ax.contour(snow.masks[name]['mask'],cmap = "Greys",linewidths = 1)

    # Do pretty stuff for the left plot
    h.axes.get_xaxis().set_ticks([])
    h.axes.get_yaxis().set_ticks([])
    divider = make_axes_locatable(ax)
    cax = divider.append_axes("right", size="5%", pad=0.2)
    cbar = plt.colorbar(h, cax = cax)
    cbar.set_label('[%s]'%(snow.depthlbl))

    h.axes.set_title('SWE \n %s'%(snow.report_date.date().strftime("%Y-%-m-%-d")))

    # Get basin-specific lims
    lims = plotlims(snow.basin, snow.plotorder)

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

    swe = pd.DataFrame(index = snow.edges, columns = snow.plotorder)

    for bid in snow.plotorder:
        r2 = database.database.query(snow, snow.start_date, snow.end_date,
                                    snow.run_name, bid, 'swe_vol')

        for elev in snow.edges:
            v2 = r2[(r2['elevation'] == str(elev))
                  & (r2['date_time'] == snow.end_date)
                  & (r2['basin_id'] == Basins.basins[bid]['basin_id'])]
            swe.loc[elev,bid] = v2['value'].values[0]

    ix = len(snow.barcolors)
    # Total basin label

    if len(snow.plotorder) > 1:
        sumorder = snow.plotorder[1::]
    else:
        sumorder = snow.plotorder

    for iters,name in enumerate(sumorder):

        if snow.dplcs == 0:
            ukaf = str(np.int(np.nansum(swe[name])))
        else:
            ukaf = str(np.round(np.nansum(swe[name]),snow.dplcs))

        if iters == 0:
            ax1.bar(range(0,len(snow.edges)),swe[name],
                    color = snow.barcolors[iters],
                    edgecolor = 'k',label = name + ': {} {}'.format(ukaf,snow.vollbl))

        else:
            ax1.bar(range(0,len(snow.edges)),swe[name],
                    bottom = pd.DataFrame(swe[sumorder[0:iters]]).sum(axis = 1).values,
                    color = snow.barcolors[iters], edgecolor = 'k',label = name + ': {} {}'.format(ukaf,snow.vollbl))

    ylims = ax1.get_ylim()
    max = ylims[1] + ylims[1]*0.5
    min = 0
    ax1.set_ylim((min, max))

    ax1.xaxis.set_ticks(range(0,len(snow.edges)))
    plt.tight_layout()
    ax1.set_xlim((snow.xlims[0]-0.5,snow.xlims[1]+0.5))

    edges_lbl = []
    for i in range(0,len(snow.edges)):
        edges_lbl.append(str(int(snow.edges[int(i)])))

    ax1.set_xticklabels(str(i) for i in edges_lbl)
    for tick in ax1.get_xticklabels():
        tick.set_rotation(30)

    ax1.set_xlim((snow.xlims[0]-0.5,snow.xlims[1]+0.5))

    ax1.set_xlabel('elevation [{}]'.format(snow.elevlbl))
    ax1.set_ylabel('{} - per elevation band'.format(snow.vollbl))
    ax1.yaxis.set_label_position("right")
    ax1.yaxis.tick_right()
    ax1.legend(loc = 2, fontsize = 10)

    for tick in ax1.get_xticklabels():
        tick.set_rotation(30)

    fig.tight_layout()
    fig.subplots_adjust(top=0.92,wspace = 0.2)

    snow._logger.info('saving figure to %sswe_volume_%s.png'%(snow.figs_path,snow.name_append))
    plt.savefig('%sswe_volume_%s.png'%(snow.figs_path,snow.name_append))
