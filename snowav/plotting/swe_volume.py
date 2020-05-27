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


def swe_volume(masks, image, df, plotorder, lims, edges, labels, barcolors,
               clim_percents, title, depthlbl, vollbl, elevlbl, xlims, dplcs,
               figs_path, fig_name, figsize, dpi=200, show=False, logger=None):
    """
    Current SWE depth and volume.

    Note: scripts/sample_figure.py and command line tool call this
    function. If changes are made here update those calls as necessary.

    masks {dict}: basin lookup dict with masks
    image {arr}: swe array
    df {DataFrame}: swe dataframe
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

    q_min, q_max = np.nanpercentile(image, [0, clim_percents[1]])
    clims = (q_min, q_max)
    pmask = masks[plotorder[0]]['mask']
    ixo = pmask == 0
    ixz = image == 0
    image[ixz] = -1

    colorsbad = plt.cm.Set2_r(np.linspace(0., 1, 1))
    colors1 = cmocean.cm.haline_r(np.linspace(0., 1, 254))
    colors = np.vstack((colorsbad, colors1))
    mymap = mcolors.LinearSegmentedColormap.from_list('my_colormap', colors)

    image[ixo] = np.nan
    mymap.set_bad('white', 1.)

    sns.set_style('darkgrid')
    sns.set_context("notebook")

    plt.close(21)
    fig, (ax, ax1) = plt.subplots(num=21, figsize=figsize, facecolor='white',
                                  dpi=dpi, nrows=1, ncols=2)

    h = ax.imshow(image, clim=clims, cmap=mymap)

    for name in masks:
        ax.contour(masks[name]['mask'], cmap="Greys", linewidths=1)

    h.axes.get_xaxis().set_ticks([])
    h.axes.get_yaxis().set_ticks([])
    divider = make_axes_locatable(ax)
    cax = divider.append_axes("right", size="2.5%", pad=0.1)
    cbar = plt.colorbar(h, cax=cax)
    cbar.set_label('[{}]'.format(depthlbl))
    h.axes.set_title(title)

    # If there is meaningful snow-free area, include path and label
    if sum(sum(ixz)) > 1000:
        patches = [mpatches.Patch(color='grey', label='snow free')]
        ax.legend(handles=patches, bbox_to_anchor=(lims.pbbx, 0.05),
                  loc=2, borderaxespad=0.)

    for iters, name in enumerate(lims.sumorder):
        if dplcs == 0:
            ukaf = str(np.int(np.nansum(df[name].values)))
        else:
            ukaf = str(np.round(np.nansum(df[name].values), dplcs))

        if iters == 0:

            ax1.bar(range(0, len(edges)), df[name],
                    color=barcolors[iters],
                    edgecolor='k',
                    label=labels[name] + ': {} {}'.format(ukaf, vollbl))

        else:
            ax1.bar(range(0, len(edges)), df[name],
                    bottom=pd.DataFrame(df[lims.sumorder[0:iters]]).sum(axis=1).values,
                    color=barcolors[iters],
                    edgecolor='k',
                    label=labels[name] + ': {} {}'.format(ukaf, vollbl))

    ylims = ax1.get_ylim()
    ymax = ylims[1] + ylims[1] * 0.5
    ymin = 0
    ax1.set_ylim((ymin, ymax))

    # Keep these to use in cold_content()
    swe_ylims = (ymin, ymax)

    ax1.xaxis.set_ticks(range(0, len(edges)))
    ax1.set_xlim((xlims[0] - 0.5, xlims[1] - 0.5))

    edges_lbl = []
    for i in range(0, len(edges)):
        edges_lbl.append(str(int(edges[int(i)])))

    ax1.set_xticklabels(str(i) for i in edges_lbl)
    ax1.tick_params(axis='x', labelsize=xtick_fontsize)
    for tick in ax1.get_xticklabels():
        tick.set_rotation(xtick_rotation)

    ax1.set_xlabel('elevation [{}]'.format(elevlbl))
    ax1.set_ylabel('{}'.format(vollbl))
    ax1.yaxis.set_label_position("right")
    ax1.yaxis.tick_right()
    ax1.legend(loc=2, fontsize=legend_fontsize)

    fig.tight_layout()
    fig.subplots_adjust(top=0.92, wspace=0.2)

    fig_name = os.path.join(os.path.abspath(figs_path), fig_name)
    snowav.framework.figures.save_fig(fig, fig_name)

    if logger is not None:
        logger.info(' Saved: {}'.format(fig_name))

    if show:
        plt.show()

    return swe_ylims
