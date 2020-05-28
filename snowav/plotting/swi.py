import cmocean
from matplotlib import pyplot as plt
import matplotlib.colors as mcolors
import matplotlib.patches as mpatches
from mpl_toolkits.axes_grid1 import make_axes_locatable
import numpy as np
import os
import pandas as pd
import seaborn as sns

import snowav.framework.figures


def swi(masks, image, df, plotorder, lims, edges, labels, barcolors, vollbl,
        elevlbl, depthlbl, clim_percents, title, figsize, figs_path, fig_name,
        dplcs, xlims, dpi=200, logger=None):
    """ SWI figure.

    Args
    ------
    masks {dict}: basin lookup dict with masks
    image {arr}: swi array
    df {DataFrame}: swi dataframe
    plotorder {list}: basins
    edges {list}: elevation bands
    labels {list}: basin labels
    barcolors {list}: plot colors
    title {str}: plot title
    vollbl {str}: volume label
    elevlbl {str}: elevation label
    depthlbl {str}: depth label
    figsize {list}: figure size
    clim_percents {list}: percentiles for color lims
    dplcs {int}: decimal places
    xlims {list}: x lims
    figs_path {str}: figure base path
    fig_name {str}: figure name
    """

    xtick_fontsize = 10
    xtick_rotation = 45
    legend_fontsize = 8

    q_min, q_max = np.nanpercentile(image, [0, clim_percents[1]])
    clims = (0, q_max)
    colors1 = cmocean.cm.dense(np.linspace(0, 1, 255))
    colors2 = plt.cm.Set1_r(np.linspace(0, 1, 1))
    colors = np.vstack((colors2, colors1))
    mymap = mcolors.LinearSegmentedColormap.from_list('my_colormap', colors)

    pmask = masks[plotorder[0]]['mask']
    ixo = pmask == 0
    image[ixo] = np.nan
    mymap.set_bad('white', 1.)

    # Now set SWI-free to some color
    r = ~np.isnan(image)
    r[r] &= image[r] < 0.001
    image[r] = 0

    sns.set_style('darkgrid')
    sns.set_context("notebook")

    plt.close(0)
    fig, (ax, ax1) = plt.subplots(num=0, figsize=figsize, dpi=dpi, nrows=1,
                                  ncols=2)

    h = ax.imshow(image, cmap=mymap, clim=clims)

    for name in masks:
        ax.contour(masks[name]['mask'], cmap='Greys', linewidths=1)

    h.axes.get_xaxis().set_ticks([])
    h.axes.get_yaxis().set_ticks([])
    divider = make_axes_locatable(ax)
    cax = divider.append_axes("right", size="2.5%", pad=0.1)
    cbar = plt.colorbar(h, cax=cax)
    cbar.set_label('[{}]'.format(depthlbl))
    h.axes.set_title(title)

    # Plot the bars
    for iters, name in enumerate(lims.sumorder):
        if df[plotorder[0]].sum() < 9.9:
            lbl = '{}: {} {}'.format(labels[name], str(int(df[name].sum())),
                                     vollbl)
        else:
            lbl = '{}: {} {}'.format(labels[name], str(np.round(df[name].sum(),
                                                                dplcs)), vollbl)

        if iters == 0:
            ax1.bar(range(0, len(edges)), df[name], color=barcolors[iters],
                    edgecolor='k', label=lbl)

        else:
            ax1.bar(range(0, len(edges)), df[name],
                    bottom=pd.DataFrame(df[lims.sumorder[0:iters]]).sum(axis=1).values,
                    color=barcolors[iters], edgecolor='k', label=lbl)

    ax1.xaxis.set_ticks(range(0, len(edges)))
    plt.tight_layout()
    ax1.set_xlim((xlims[0] - 0.5, xlims[1] - 0.5))

    edges_lbl = []
    for i in range(0, len(edges)):
        edges_lbl.append(str(int(edges[int(i)])))

    ax1.set_xticklabels(str(i) for i in edges_lbl)
    ax1.tick_params(axis='x', labelsize=xtick_fontsize)
    for tick in ax1.get_xticklabels():
        tick.set_rotation(xtick_rotation)

    ax1.set_ylabel('{}'.format(vollbl))
    ax1.set_xlabel('elevation [{}]'.format(elevlbl))
    ax1.yaxis.set_label_position("right")
    ax1.yaxis.tick_right()

    ylims = ax1.get_ylim()
    ymax = ylims[1] + ylims[1] * 0.6
    ax1.set_ylim((0, ymax))

    fig.subplots_adjust(top=0.88)

    # If there is meaningful snow-free area, include patch and label
    if sum(sum(r)) > 1000:
        patches = [mpatches.Patch(color='grey', label='no SWI')]
        ax.legend(handles=patches, bbox_to_anchor=(lims.pbbx, 0.05),
                  loc=2, borderaxespad=0.)

    plt.tight_layout()

    if len(plotorder) > 1:
        ax1.legend(loc=2, markerscale=0.5, fontsize=legend_fontsize)

    fig_name = os.path.join(os.path.abspath(figs_path), fig_name)
    snowav.framework.figures.save_fig(fig, fig_name)

    if logger is not None:
        logger.info(' Saved: {}'.format(fig_name))

