
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

    patches = [mpatches.Patch(color='grey', label='snow free')]

    # Make legend/box defaults and adjust as needed
    pbbx = 0.15
    pbby = 0.05

    if snow.basin in ['KAWEAH', 'RCEW']:
        pbbx = 0.1

    # snow-free
    ax.legend(handles=patches, bbox_to_anchor=(pbbx, pbby),
              loc=2, borderaxespad=0. )

    swe = pd.DataFrame(index = snow.edges, columns = snow.plotorder)

    for bid in snow.plotorder:
        r2 = database.database.query(snow, snow.start_date, snow.end_date,
                                    snow.run_name, bid, 'swe_vol')

        for elev in snow.edges:
            v2 = r2[(r2['elevation'] == str(elev))
                  & (r2['date_time'] == snow.end_date)
                  & (r2['basin_id'] == Basins.basins[bid]['basin_id'])]
            swe.loc[elev,bid] = v2['value'].values[0]

    lim = np.max(np.max(swe[snow.plotorder[0]]))
    if lim <= 1:
        lim = 1

    ylim = np.max(lim) + np.max(lim)*0.2
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

    ax1.set_xlim((snow.xlims[0]-0.5,snow.xlims[1]))
    plt.tight_layout()

    xts = ax1.get_xticks()
    if len(xts) < 6:
        dxt = xts[1] - xts[0]
        xts = np.arange(xts[0],xts[-1] + 1 ,dxt/2)
        ax1.set_xticks(xts)

    edges_lbl = []
    for i in xts[0:len(xts)-1]:
        edges_lbl.append(str(int(snow.edges[int(i)])))

    ax1.set_xticklabels(str(i) for i in edges_lbl)
    ax1.set_xlabel('elevation [{}]'.format(snow.elevlbl))
    ax1.set_ylabel('{} - per elevation band'.format(snow.vollbl))
    ax1.yaxis.set_label_position("right")
    ax1.yaxis.tick_right()
    ax1.legend(loc = 2, fontsize = 10)
    ax1.set_ylim((0,ylim))

    for tick in ax1.get_xticklabels():
        tick.set_rotation(30)

    fig.tight_layout()
    fig.subplots_adjust(top=0.92,wspace = 0.2)

    snow._logger.info('saving figure to %sswe_volume_%s.png'%(snow.figs_path,snow.name_append))
    plt.savefig('%sswe_volume_%s.png'%(snow.figs_path,snow.name_append))
