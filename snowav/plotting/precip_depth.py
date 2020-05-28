import cmocean
import copy
from matplotlib import pyplot as plt
import matplotlib.colors as mcolors
import matplotlib.patches as mpatches
from mpl_toolkits.axes_grid1 import make_axes_locatable
import numpy as np
import os
import seaborn as sns

import snowav.framework.figures


def precip_depth(swi, precip, rain, swi_byelev, precip_byelev, rain_byelev,
                 plotorder, barcolors, lims, masks, edges, labels,
                 clims_percent, depthlbl, elevlbl, title, xlims,
                 figs_path, fig_name, logger=None):
    """ Depth of SWI, precipitation, and rain during the report period.

    Args
    ------
    swi {}: swi image
    precip {}: precip image
    rain {}: rain image
    swi_byelev {DataFrame}: by elevation dataframe
    precip_byelev {DataFrame}: by elevation dataframe
    rain_byelev {DataFrame}:by elevation dataframe
    plotorder {list}: basins
    masks {list}: basin lookup dict with masks
    edges {list}: elevation bands
    labels {list}: basin labels
    barcolors {list}: plot colors
    title {str}: plot title
    depthlbl {str}: depth label
    elevlbl {str}: elevation label
    xlims {list}: x lims
    figs_path {str}: figure base path
    fig_name {str}: figure name
    """

    xtick_fontsize = 10
    xtick_rotation = 45
    legend_fontsize = 8

    if np.nanpercentile(swi, clims_percent)[1] > np.nanpercentile(precip, clims_percent)[1]:
        z, q_max = np.nanpercentile(swi, clims_percent)
    else:
        z, q_max = np.nanpercentile(precip, clims_percent)

    clims = (-0.1, q_max)
    y_min = -0.02

    # Get bar plot ylims
    if np.nanmax(swi_byelev.values) > np.nanmax(precip_byelev.values):
        if len(plotorder) < 5:
            y_max = np.nanmax(swi_byelev.values) + np.nanmax(swi_byelev.values) * 0.4
        else:
            y_max = np.nanmax(swi_byelev.values) + np.nanmax(swi_byelev.values) * 0.6
    else:
        if len(plotorder) < 5:
            y_max = np.nanmax(precip_byelev.values) + np.nanmax(precip_byelev.values) * 0.4
        else:
            y_max = np.nanmax(precip_byelev.values) + np.nanmax(precip_byelev.values) * 0.6

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
    fig, ax = plt.subplots(num=0, figsize=(12, 10), dpi=200, nrows=3, ncols=2)

    ################################################
    #           SWI                                #
    ################################################
    mymap = copy.deepcopy(cmap)
    pmask = masks[plotorder[0]]['mask']
    ixo = pmask == 0
    swi[ixo] = np.nan
    mymap.set_bad('white', 1.)

    r = ~np.isnan(swi)
    r[r] &= swi[r] < 0.001
    swi[r] = 0

    h = ax[0, 0].imshow(swi, cmap=mymap, clim=clims)

    for name in masks:
        ax[0, 0].contour(masks[name]['mask'], cmap='Greys', linewidths=1)

    h.axes.get_xaxis().set_ticks([])
    h.axes.get_yaxis().set_ticks([])
    divider = make_axes_locatable(ax[0, 0])
    cax = divider.append_axes("right", size="2.5%", pad=0.1)
    cbar = plt.colorbar(h, cax=cax)
    cbar.set_label('SWI [{}]'.format(depthlbl))

    if len(plotorder) == 1:
        sumorder = plotorder
        swid = 0.4
        wid = [-0.1, 0.1]
    elif len(plotorder) <= 4:
        sumorder = plotorder[1::]
        swid = 0.25
        wid = np.linspace(-0.3, 0.3, len(sumorder))
    elif len(plotorder) == 5:
        sumorder = plotorder[1::]
        swid = 0.2
        wid = np.linspace(-0.3, 0.3, len(sumorder))
    elif len(plotorder) > 5:
        sumorder = plotorder[1::]
        swid = 0.1
        wid = np.linspace(-0.4, 0.4, len(sumorder))

    for iters, name in enumerate(sumorder):
        lbl = labels[name]
        if len(plotorder) == 1:
            ax[0, 1].bar(range(0, len(edges)),
                         swi_byelev[name],
                         color=barcolors[iters], width=swid, edgecolor='k',
                         label=lbl)
        else:
            ax[0, 1].bar(range(0, len(edges)) - wid[iters],
                         swi_byelev[name],
                         color=barcolors[iters], width=swid, edgecolor='k',
                         label=lbl)

    ax[0, 1].xaxis.set_ticks(range(0, len(edges)))
    ax[0, 1].set_xticklabels('')
    ax[0, 1].set_xlim((xlims[0] - 0.5, xlims[1] - 0.5))
    ax[0, 1].set_ylabel('[{}]'.format(depthlbl))
    ax[0, 1].yaxis.set_label_position("right")
    ax[0, 1].yaxis.tick_right()
    ax[0, 1].set_ylim((y_min, y_max))

    if sum(sum(r)) > 1000:
        patches = [mpatches.Patch(color='grey', label='snow free')]
        ax[0, 0].legend(handles=patches, bbox_to_anchor=(lims.pbbx, 0.05),
                        loc=2, borderaxespad=0.)

    if len(plotorder) > 1:
        ax[0, 1].legend(loc=2, markerscale=0.5, fontsize=legend_fontsize)

    ################################################
    #           Precip                             #
    ################################################
    mymap1 = copy.deepcopy(cmap1)
    pmask = masks[plotorder[0]]['mask']
    ixo = pmask == 0
    precip[ixo] = np.nan
    mymap1.set_bad('white', 1.)
    r = ~np.isnan(precip)
    r[r] &= precip[r] < 0.001
    precip[r] = 0

    h2 = ax[1, 0].imshow(precip, interpolation='none', cmap=mymap1, clim=clims)

    for name in masks:
        ax[1, 0].contour(masks[name]['mask'], cmap="Greys", linewidths=1)

    h2.axes.get_xaxis().set_ticks([])
    h2.axes.get_yaxis().set_ticks([])
    divider = make_axes_locatable(ax[1, 0])
    cax = divider.append_axes("right", size="2.5%", pad=0.1)
    cbar = plt.colorbar(h2, cax=cax)
    cbar.set_label('precipitation [{}]'.format(depthlbl))

    for iters, name in enumerate(sumorder):
        lbl = name

        if len(plotorder) == 1:
            ax[1, 1].bar(range(0, len(edges)),
                         precip_byelev[name],
                         color=barcolors[iters], width=swid, edgecolor='k',
                         label=lbl)
        else:
            ax[1, 1].bar(range(0, len(edges)) - wid[iters],
                         precip_byelev[name],
                         color=barcolors[iters], width=swid, edgecolor='k',
                         label=lbl)

        ax[1, 1].set_xlim((xlims[0] - 0.5, xlims[1] - 0.5))

    ax[1, 1].xaxis.set_ticks(range(0, len(edges)))
    ax[1, 1].set_xlim((xlims[0] - 0.5, xlims[1] - 0.5))
    ax[1, 1].set_xticklabels('')
    ax[1, 1].set_ylim((y_min, y_max))
    ax[1, 1].set_ylabel('[{}]'.format(depthlbl))
    ax[1, 1].yaxis.set_label_position("right")
    ax[1, 1].yaxis.tick_right()

    ################################################
    #           rain                                #
    ################################################
    ix = rain_byelev < 0
    rain_byelev[ix] = 0

    mymap = copy.deepcopy(cmap1)
    pmask = masks[plotorder[0]]['mask']
    ixo = pmask == 0
    rain[ixo] = np.nan
    mymap.set_bad('white', 1.)
    r = ~np.isnan(rain)
    r[r] &= rain[r] < 0.001
    rain[r] = 0

    h2 = ax[2, 0].imshow(rain, interpolation='none', cmap=mymap, clim=clims)

    for name in masks:
        ax[2, 0].contour(masks[name]['mask'], cmap="Greys", linewidths=1)

    h2.axes.get_xaxis().set_ticks([])
    h2.axes.get_yaxis().set_ticks([])
    divider = make_axes_locatable(ax[2, 0])
    cax = divider.append_axes("right", size="2.5%", pad=0.1)
    cbar = plt.colorbar(h2, cax=cax)
    cbar.set_label('rain [{}]'.format(depthlbl))

    for iters, name in enumerate(sumorder):
        lbl = name
        if len(plotorder) == 1:
            ax[2, 1].bar(range(0, len(edges)),
                         rain_byelev[name],
                         color=barcolors[iters], width=swid, edgecolor='k',
                         label=lbl)
        else:
            ax[2, 1].bar(range(0, len(edges)) - wid[iters],
                         rain_byelev[name],
                         color=barcolors[iters], width=swid, edgecolor='k',
                         label=lbl)

        ax[2, 1].set_xlim((xlims[0] - 0.5, xlims[1] - 0.5))

    ax[2, 1].xaxis.set_ticks(range(0, len(edges)))

    edges_lbl = []
    for i in range(0, len(edges)):
        edges_lbl.append(str(int(edges[int(i)])))

    ax[2, 1].set_xticklabels(str(i) for i in edges_lbl)
    ax[2, 1].tick_params(axis='x', labelsize=xtick_fontsize)
    for tick in ax[2, 1].get_xticklabels():
        tick.set_rotation(xtick_rotation)

    ax[2, 1].set_ylim((y_min, y_max))
    ax[2, 1].set_ylabel('[{}]'.format(depthlbl))
    ax[2, 1].set_xlabel('elevation [{}]'.format(elevlbl))
    ax[2, 1].yaxis.set_label_position("right")
    ax[2, 1].yaxis.tick_right()

    plt.suptitle(title)
    fig.subplots_adjust(top=0.92)

    fig_name = os.path.join(os.path.abspath(figs_path), fig_name)
    snowav.framework.figures.save_fig(fig, fig_name)

    if logger is not None:
        logger.info(' Saved: {}'.format(fig_name))

