
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
from snowav.plotting.plotlims import plotlims as plotlims
import pandas as pd
from matplotlib.ticker import FormatStrFormatter


def accumulated(snow):
    '''
    Figure: 0

    outputs['swi_z'] (trimmed to ixe,ixs)
    depth_factor
    edges
    plotorder
    figsize
    start_date
    end_date
    report_date
    run_name
    session (to database)
    masks
    imgx (lakes)
    imgy (lakes)
    depthlbl
    vollbl
    figs_path
    name_append
    xlims

    '''

    # Calculate accumulated swi during the specified period
    accum = np.zeros((len(snow.outputs['swi_z'][0][:,0]),
                      len(snow.outputs['swi_z'][0][0,:])))

    for n in range(snow.ixs,snow.ixe):
        accum = accum + np.multiply(snow.outputs['swi_z'][n],snow.depth_factor)

    # Make df from database
    accum_byelev = pd.DataFrame(index = snow.edges, columns = snow.plotorder)

    for bid in snow.plotorder:
        r = database.database.query(snow, snow.start_date, snow.end_date,
                                    snow.run_name, bid, 'swi_vol')

        for elev in snow.edges:
            v = r[r['elevation'] == str(elev)]
            accum_byelev.loc[elev,bid] = np.nansum(v['value'].values)

    #######################################
    #            plot set up              #
    #######################################

    qMin,qMax = np.nanpercentile(accum,[0,99.8])
    clims = (0,qMax)
    colors1 = cmocean.cm.dense(np.linspace(0., 1, 255))
    colors2 = plt.cm.binary(np.linspace(0, 1, 1))
    colors = np.vstack((colors2, colors1))
    mymap = mcolors.LinearSegmentedColormap.from_list('my_colormap', colors)

    # White background
    pmask = snow.masks[snow.plotorder[0]]['mask']
    ixo = pmask == 0
    accum[ixo] = np.nan
    mymap.set_bad('white',1.)

    # Now set SWI-free to some color
    zs = ~np.isnan(accum)
    zs[zs] &= accum[zs] < 0.02
    accum[zs] = -1
    mymap.set_under('grey',1.)

    # Get basin-specific lims
    lims = plotlims(snow.basin, snow.plotorder)

    sns.set_style('darkgrid')
    sns.set_context("notebook")

    plt.close(0)
    fig,(ax,ax1) = plt.subplots(num=0, figsize=snow.figsize,
                                dpi=snow.dpi, nrows=1, ncols=2)
    h = ax.imshow(accum, clim=clims, cmap=mymap)

    # Basin boundaries
    for name in snow.masks:
        ax.contour(snow.masks[name]['mask'], cmap='Greys',linewidths=1)

    if snow.basin == 'LAKES':
        ax.set_xlim(snow.imgx)
        ax.set_ylim(snow.imgy)

    # Do pretty stuff
    h.axes.get_xaxis().set_ticks([])
    h.axes.get_yaxis().set_ticks([])
    divider = make_axes_locatable(ax)
    cax = divider.append_axes("right", size="4%", pad=0.2)
    cbar = plt.colorbar(h, cax = cax)
    cbar.set_label('[{}]'.format(snow.depthlbl))
    h.axes.set_title('Accumulated SWI \n {} to {}'.format(
                                snow.report_start.date().strftime("%Y-%-m-%-d"),
                                snow.report_date.date().strftime("%Y-%-m-%-d")))

    # Decimal places
    if accum_byelev[snow.plotorder[0]].sum() < 9.9:
        tlbl = '{}: {} {}'.format(snow.plotorder[0],
                             str(int(accum_byelev[snow.plotorder[0]].sum())),
                             snow.vollbl)
    else:
        tlbl = '{}: {} {}'.format(snow.plotorder[0],
                             str(np.round(accum_byelev[snow.plotorder[0]].sum(),
                            snow.dplcs)),snow.vollbl)

    # Plot the bars
    for iters,name in enumerate(lims.sumorder):
        if accum_byelev[snow.plotorder[0]].sum() < 9.9:
            lbl = '{}: {} {}'.format(name,str(int(accum_byelev[name].sum())),
                                snow.vollbl)
        else:
            lbl = '{}: {} {}'.format(name,str(np.round(accum_byelev[name].sum(),
                                snow.dplcs)),snow.vollbl)

        if iters == 0:
            ax1.bar(range(0,len(snow.edges)),accum_byelev[name],
                    color = snow.barcolors[iters],
                    edgecolor = 'k',label = lbl)

        else:
            ax1.bar(range(0,len(snow.edges)),accum_byelev[name],
                    bottom = pd.DataFrame(accum_byelev[lims.sumorder[0:iters]]).sum(axis = 1).values,
                    color = snow.barcolors[iters], edgecolor = 'k',label = lbl)

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

    ax1.set_ylabel('{} - per elevation band'.format(snow.vollbl))
    ax1.set_xlabel('elevation [{}]'.format(snow.elevlbl))
    #ax1.set_xlim((snow.xlims[0]-0.5,snow.xlims[1]+0.5))
    ax1.yaxis.set_label_position("right")
    ax1.yaxis.tick_right()

    ylims = ax1.get_ylim()
    max = ylims[1] + ylims[1]*0.6
    min = 0
    ax1.set_ylim((min, max))

    #plt.tight_layout()
    fig.subplots_adjust(top=0.88)

    # If there is meaningful snow-free area, include path and label
    if sum(sum(zs)) > 1000:
        patches = [mpatches.Patch(color='grey', label='no SWI')]
        ax.legend(handles=patches, bbox_to_anchor=(lims.pbbx, 0.05),
                  loc=2, borderaxespad=0. )

    # fix so that legend isn't obscured
    if accum_byelev.iloc[0,:].sum() > np.nanmax(accum_byelev.sum())*0.9:

        if len(snow.plotorder) > 1:
            ax1.legend(loc=1,markerscale = 0.5)

    else:

        if len(snow.plotorder) > 1:
            ax1.legend(loc=(lims.legx,lims.legy),markerscale = 0.5)

        ax1.text(lims.btx,lims.bty,tlbl,horizontalalignment='center',
                 transform=ax1.transAxes,fontsize = 10)

    snow._logger.info('saving figure to{}swi_{}.png'.format(snow.figs_path,
                                                             snow.name_append))
    plt.tight_layout()
    plt.savefig('{}swi_{}.png'.format(snow.figs_path,snow.name_append))
