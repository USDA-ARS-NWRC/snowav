
from snowav.utils.MidpointNormalize import MidpointNormalize
import numpy as np
from matplotlib import pyplot as plt
import matplotlib.colors as mcolors
from mpl_toolkits.axes_grid1 import make_axes_locatable
import seaborn as sns
import cmocean
import matplotlib.patches as mpatches
import pandas as pd
import copy
import snowav.framework.figures

def precip_depth(args, logger = None):
    '''
    Depth of SWI, precipitation, and rain during the report period.

    Args
    ----------
    args : dict
        dictionary with required inputs, see swi() figure for more information.

    logger : list
        snowav logger

    '''

    accum = args['swi_image']
    precip = args['precip_image']
    rain = args['rain_image']
    accum_byelev = args['swi_df']
    precip_byelev = args['precip_df']
    rain_byelev = args['rain_df']
    plotorder = args['plotorder']
    barcolors = args['barcolors']
    lims = args['lims']
    masks = args['masks']
    edges = args['edges']
    labels = args['labels']

    mm = [args['percent_min'],args['percent_max']]

    if np.nanpercentile(accum,mm)[1] > np.nanpercentile(precip,mm)[1]:
        z,qMax = np.nanpercentile(accum,mm)
    else:
        z,qMax = np.nanpercentile(precip,mm)

    if np.nanpercentile(accum,mm)[0] < np.nanpercentile(precip,mm)[0]:
        qMin,z = np.nanpercentile(accum,mm)
    else:
        qMin,z = np.nanpercentile(precip,mm)

    clims = (0,qMax)

    # Get bar plot ylims
    if np.nanmax(accum_byelev.values) > np.nanmax(precip_byelev.values):
        if len(plotorder) < 5:
            yMax = np.nanmax(accum_byelev.values) + np.nanmax(accum_byelev.values)*0.4
        else:
            yMax = np.nanmax(accum_byelev.values) + np.nanmax(accum_byelev.values)*0.6
    else:
        if len(plotorder) < 5:
            yMax = np.nanmax(precip_byelev.values) + np.nanmax(precip_byelev.values)*0.4
        else:
            yMax = np.nanmax(precip_byelev.values) + np.nanmax(precip_byelev.values)*0.6

    colors1 = cmocean.cm.dense(np.linspace(0, 1, 255))
    colors2 = plt.cm.Set1_r(np.linspace(0, 1, 1))
    colors = np.vstack((colors2, colors1))
    cmap = mcolors.LinearSegmentedColormap.from_list('my_colormap', colors)

    colors1 = plt.cm.nipy_spectral_r(np.linspace(0, 1, 255))
    colors2 = plt.cm.Set1_r(np.linspace(0, 1, 1))
    colors = np.vstack((colors2, colors1))
    cmap1 = mcolors.LinearSegmentedColormap.from_list('my_colormap', colors)

    sns.set_style('darkgrid')
    sns.set_context("notebook")

    plt.close(0)
    fig, ax = plt.subplots(num = 0, figsize = (12,10),dpi=args['dpi'],
                           nrows = 3, ncols = 2)

    ################################################
    #           SWI                                #
    ################################################
    mymap = copy.deepcopy(cmap)
    pmask = masks[plotorder[0]]['mask']
    ixo = pmask == 0
    accum[ixo] = np.nan
    mymap.set_bad('white',1.)

    r = ~np.isnan(accum)
    r[r] &= accum[r] < 0.001
    accum[r] = 0

    h = ax[0,0].imshow(accum,  cmap = mymap, clim=clims)

    for name in masks:
        ax[0,0].contour(masks[name]['mask'],cmap = 'Greys',linewidths = 1)

    h.axes.get_xaxis().set_ticks([])
    h.axes.get_yaxis().set_ticks([])
    divider = make_axes_locatable(ax[0,0])
    cax = divider.append_axes("right", size="4%", pad=0.2)
    cbar = plt.colorbar(h, cax = cax)
    cbar.set_label('[{}]'.format(args['depthlbl']))
    h.axes.set_title('Accumulated SWI')
    xlims = args['xlims']
    if len(plotorder) == 1:
        sumorder = plotorder
        swid = 0.4
        wid = [-0.1, 0.1]
    elif len(plotorder) <= 4:
        sumorder = plotorder[1::]
        swid = 0.25
        wid = np.linspace(-0.3,0.3,len(sumorder))
    elif len(plotorder) == 5:
        sumorder = plotorder[1::]
        swid = 0.2
        wid = np.linspace(-0.3,0.3,len(sumorder))
    elif len(plotorder) > 5:
        sumorder = plotorder[1::]
        swid = 0.1
        wid = np.linspace(-0.4,0.4,len(sumorder))

    for iters,name in enumerate(sumorder):
        lbl = labels[name]
        if len(plotorder) == 1:
            ax[0,1].bar(range(0,len(edges)),
                    accum_byelev[name],
                    color = barcolors[iters], width = swid, edgecolor = 'k', label = lbl)
        else:
            ax[0,1].bar(range(0,len(edges))-wid[iters],
                    accum_byelev[name],
                    color = barcolors[iters], width = swid, edgecolor = 'k', label = lbl)

    ax[0,1].xaxis.set_ticks(range(0,len(edges)))
    ax[0,1].set_xlim((xlims[0]-0.5,xlims[1]-0.5))

    edges_lbl = []
    for i in range(0,len(edges)):
        edges_lbl.append(str(int(edges[int(i)])))

    ax[0,1].set_xticklabels(str(i) for i in edges_lbl)
    for tick in ax[0,1].get_xticklabels():
        tick.set_rotation(30)

    ax[0,1].set_xlim((xlims[0]-0.5,xlims[1]-0.5))
    ax[0,1].set_ylabel('[{}]'.format(args['depthlbl']))
    ax[0,1].yaxis.set_label_position("right")
    ax[0,1].yaxis.tick_right()
    ax[0,1].set_ylim((0,yMax))

    if sum(sum(r)) > 1000:
        patches = [mpatches.Patch(color='grey', label='no SWI')]
        ax[0,0].legend(handles=patches, bbox_to_anchor=(lims.pbbx, 0.05),
                  loc=2, borderaxespad=0. )

    if len(plotorder) > 1:
        ax[0,1].legend(loc=(lims.legx,lims.legy2),markerscale = 0.5)

    ################################################
    #           Precip                             #
    ################################################

    mymap1 = copy.deepcopy(cmap1)
    pmask = masks[plotorder[0]]['mask']
    ixo = pmask == 0
    precip[ixo] = np.nan
    mymap1.set_bad('white',1.)
    r = ~np.isnan(precip)
    r[r] &= precip[r] < 0.001
    precip[r] = 0

    h2 = ax[1,0].imshow(precip, interpolation='none', cmap = mymap1, clim = clims)

    for name in masks:
        ax[1,0].contour(masks[name]['mask'],cmap = "Greys",linewidths = 1)

    h2.axes.get_xaxis().set_ticks([])
    h2.axes.get_yaxis().set_ticks([])
    divider = make_axes_locatable(ax[1,0])
    cax = divider.append_axes("right", size="5%", pad=0.2)
    cbar = plt.colorbar(h2, cax = cax)
    cbar.set_label(r'[{}]'.format(args['depthlbl']))

    h2.axes.set_title('Precipitation')

    for iters,name in enumerate(sumorder):
        lbl = name

        if len(plotorder) == 1:
            ax[1,1].bar(range(0,len(edges)),
                    precip_byelev[name],
                    color = barcolors[iters], width = swid, edgecolor = 'k', label = lbl)
        else:
            ax[1,1].bar(range(0,len(edges))-wid[iters],
                    precip_byelev[name],
                    color = barcolors[iters], width = swid, edgecolor = 'k', label = lbl)

        ax[1,1].set_xlim((xlims[0]-0.5,xlims[1]-0.5))

    ax[1,1].xaxis.set_ticks(range(0,len(edges)))
    ax[1,1].set_xlim((xlims[0]-0.5,xlims[1]-0.5))

    edges_lbl = []
    for i in range(0,len(edges)):
        edges_lbl.append(str(int(edges[int(i)])))

    ax[1,1].set_xticklabels(str(i) for i in edges_lbl)
    for tick in ax[1,1].get_xticklabels():
        tick.set_rotation(30)

    ax[1,1].set_xlim((xlims[0]-0.5,xlims[1]-0.5))
    ax[1,1].set_ylim((0,yMax))
    ax[1,1].set_ylabel(r'[{}]'.format(args['depthlbl']))
    ax[1,1].yaxis.set_label_position("right")
    ax[1,1].tick_params(axis='x')
    ax[1,1].tick_params(axis='y')
    ax[1,1].yaxis.tick_right()

    patches = [mpatches.Patch(color='grey', label='no precip')]

    ################################################
    #           rain                                #
    ################################################

    ix = rain_byelev < 0
    rain_byelev[ix] = 0

    mymap = copy.deepcopy(cmap1)
    pmask = masks[plotorder[0]]['mask']
    ixo = pmask == 0
    rain[ixo] = np.nan
    mymap.set_bad('white',1.)
    r = ~np.isnan(rain)
    r[r] &= rain[r] < 0.001
    rain[r] = 0

    h2 = ax[2,0].imshow(rain, interpolation='none', cmap = mymap, clim = clims)

    for name in masks:
        ax[2,0].contour(masks[name]['mask'],cmap = "Greys",linewidths = 1)

    h2.axes.get_xaxis().set_ticks([])
    h2.axes.get_yaxis().set_ticks([])
    divider = make_axes_locatable(ax[2,0])
    cax = divider.append_axes("right", size="5%", pad=0.2)
    cbar = plt.colorbar(h2, cax = cax)
    cbar.set_label(r'[{}]'.format(args['depthlbl']))

    h2.axes.set_title('Rain')

    for iters,name in enumerate(sumorder):
        lbl = name
        if len(plotorder) == 1:
            ax[2,1].bar(range(0,len(edges)),
                    rain_byelev[name],
                    color = barcolors[iters], width = swid, edgecolor = 'k', label = lbl)
        else:
            ax[2,1].bar(range(0,len(edges))-wid[iters],
                    rain_byelev[name],
                    color = barcolors[iters], width = swid, edgecolor = 'k', label = lbl)

        ax[2,1].set_xlim((xlims[0]-0.5,xlims[1]-0.5))

    ax[2,1].xaxis.set_ticks(range(0,len(edges)))
    ax[2,1].set_xlim((xlims[0]-0.5,xlims[1]-0.5))

    edges_lbl = []
    for i in range(0,len(edges)):
        edges_lbl.append(str(int(edges[int(i)])))

    ax[2,1].set_xticklabels(str(i) for i in edges_lbl)
    for tick in ax[2,1].get_xticklabels():
        tick.set_rotation(30)

    ax[2,1].set_xlim((xlims[0]-0.5,xlims[1]-0.5))
    ax[2,1].set_ylim((0,yMax))
    ax[2,1].set_ylabel(r'[{}]'.format(args['depthlbl']))
    ax[2,1].set_xlabel('elevation [{}]'.format(args['elevlbl']))
    ax[2,1].yaxis.set_label_position("right")
    ax[2,1].tick_params(axis='x')
    ax[2,1].tick_params(axis='y')
    ax[2,1].yaxis.tick_right()

    plt.suptitle(args['title'])
    fig.subplots_adjust(top=0.92)

    fig_name_short = 'precip_depth_'
    fig_name = '{}{}{}.png'.format(args['figs_path'],fig_name_short,args['directory'])
    if logger is not None:
        logger.info(' saving {}'.format(fig_name))
    snowav.framework.figures.save_fig(fig, fig_name)

    return fig_name_short
