import cmocean
import copy
from matplotlib import pyplot as plt
import matplotlib.colors as mcolors
from mpl_toolkits.axes_grid1 import make_axes_locatable
import numpy as np
import os
import pandas as pd
import seaborn as sns

import snowav.framework.figures
from snowav.utils.MidpointNormalize import MidpointNormalize


def image_change(masks, image, df, plotorder, lims, edges, labels, barcolors,
                 clims_percent, vollbl, elevlbl, depthlbl, xlims, figsize,
                 title, dplcs, figs_path, fig_name, dpi=200, logger=None):
    """ Change in SWE depth and volume.

    Args
    ------
    masks {dict}: basin lookup dict with masks
    image {arr}: change in swe image
    df {DataFrame}: swe change dataframe
    plotorder {list}: basins
    edges {list}: elevation bands
    labels {list}: basin labels
    barcolors {list}: plot colors
    title {str}: plot title
    vollbl {str}: volume label
    elevlbl {str}: elevation label
    figsize {list}: figure size
    dplcs {int}: decimal places
    xlims {list}: x lims
    figs_path {str}: figure base path
    fig_name {str}: figure name
    """

    xtick_fontsize = 10
    xtick_rotation = 45
    legend_fontsize = 8

    v_min, v_max = np.nanpercentile(image, clims_percent)
    colorsbad = plt.cm.Set1_r(np.linspace(0., 1, 1))
    colors1 = cmocean.cm.matter_r(np.linspace(0., 1, 127))
    colors2 = plt.cm.Blues(np.linspace(0, 1, 128))
    colors = np.vstack((colorsbad, colors1, colors2))
    mymap = mcolors.LinearSegmentedColormap.from_list('my_colormap', colors)

    ixf = image == 0
    image[ixf] = -100000
    pmask = masks[plotorder[0]]['mask']
    ixo = pmask == 0
    image[ixo] = np.nan
    cmap = copy.copy(mymap)
    cmap.set_bad('white', 1.)

    sns.set_style('darkgrid')
    sns.set_context("notebook")

    plt.close(6)
    fig, (ax, ax1) = plt.subplots(num=6, figsize=figsize, dpi=dpi, nrows=1,
                                  ncols=2)
    h = ax.imshow(image,
                  interpolation='none',
                  cmap=cmap,
                  norm=MidpointNormalize(midpoint=0,
                                         vmin=v_min - 0.01,
                                         vmax=v_max + 0.01))

    for name in masks:
        ax.contour(masks[name]['mask'], cmap="Greys", linewidths=1)

    h.axes.get_xaxis().set_ticks([])
    h.axes.get_yaxis().set_ticks([])
    divider = make_axes_locatable(ax)
    cax = divider.append_axes("right", size="2.5%", pad=0.1)
    cbar = plt.colorbar(h, cax=cax)
    cbar.set_label(r'$\Delta$ SWE [{}]'.format(depthlbl))

    h.axes.set_title(title)

    for iters, name in enumerate(lims.sumorder):
        if dplcs == 0:
            lbl = '{}: {} {}'.format(labels[name], str(int(df[name].sum())),
                                     vollbl)
        else:
            lbl = '{}: {} {}'.format(labels[name], str(np.round(df[name].sum(),
                                                                dplcs)), vollbl)

        if iters == 0:
            ax1.bar(range(0, len(edges)), df[name],
                    color=barcolors[iters],
                    edgecolor='k', label=lbl)

        else:
            ax1.bar(range(0, len(edges)), df[name],
                    bottom=pd.DataFrame(df[lims.sumorder[0:iters]]).sum(axis=1).values,
                    color=barcolors[iters], edgecolor='k', label=lbl)

    ax1.xaxis.set_ticks(range(0, len(edges)))

    edges_lbl = []
    for i in range(0, len(edges)):
        edges_lbl.append(str(int(edges[int(i)])))

    ax1.set_xticklabels(str(i) for i in edges_lbl)
    ax1.tick_params(axis='x', labelsize=xtick_fontsize)
    for tick in ax1.get_xticklabels():
        tick.set_rotation(xtick_rotation)

    ax1.set_xlim((xlims[0] - 0.5, xlims[1] - 0.5))

    ylims = ax1.get_ylim()
    ymax = ylims[1] + abs(ylims[1] - ylims[0])
    ymin = ylims[0] - abs(ylims[1] - ylims[0]) * 0.1
    ax1.set_ylim((ymin, ymax))

    ax1.set_ylabel('{}'.format(vollbl))
    ax1.set_xlabel('elevation [{}]'.format(elevlbl))
    ax1.axes.set_title('Change in Volume')
    ax1.yaxis.set_label_position("right")
    ax1.yaxis.tick_right()

    if len(plotorder) > 1:
        ax1.legend(loc=2, fontsize=legend_fontsize)

    plt.tight_layout()
    fig.subplots_adjust(top=0.88)

    fig_name = os.path.join(os.path.abspath(figs_path), fig_name)
    snowav.framework.figures.save_fig(fig, fig_name)

    if logger is not None:
        logger.info(' Saved: {}'.format(fig_name))
