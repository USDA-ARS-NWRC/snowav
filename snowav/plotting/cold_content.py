from matplotlib import pyplot as plt
from mpl_toolkits.axes_grid1 import make_axes_locatable
import numpy as np
import os
import pandas as pd
import seaborn as sns

import snowav.framework.figures


def cold_content(masks, swe, image, df, plotorder, lims, edges, labels,
                 barcolors, title, vollbl, elevlbl, figsize, dplcs, xlims,
                 figs_path, fig_name, ylims=None, dpi=200, logger=None):
    """ SWE unavailable for melt based on snowpack cold content.

    Args
    ------
    masks {dict}: basin lookup dict with masks
    swe {arr}: swe array
    image {arr}: cold content array
    df {DataFrame}: cold content dataframe
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
    ylims {list}: ylims
    """

    xtick_fontsize = 10
    xtick_rotation = 45
    legend_fontsize = 8

    clims2 = (-5, 0)
    pmask = masks[plotorder[0]]['mask']
    ixo = pmask == 0
    ixz = swe == 0
    image[ixz] = 1

    mymap1 = plt.cm.Spectral_r
    image[ixo] = np.nan
    mymap1.set_bad('white')
    mymap1.set_over('white', 1)

    sns.set_style('darkgrid')
    sns.set_context("notebook")

    plt.close(1)
    fig, (ax, ax1) = plt.subplots(num=1, figsize=figsize, facecolor='white',
                                  dpi=dpi, nrows=1, ncols=2)

    h = ax.imshow(image, clim=clims2, cmap=mymap1)

    for name in masks:
        ax.contour(masks[name]['mask'], cmap="Greys", linewidths=1)

    h.axes.get_xaxis().set_ticks([])
    h.axes.get_yaxis().set_ticks([])
    h.axes.set_title(title)
    divider = make_axes_locatable(ax)
    cax2 = divider.append_axes("right", size="2.5%", pad=0.1)
    cbar = plt.colorbar(h, cax=cax2)
    cbar.set_label('[MJ/$m^3$]')
    cbar.ax.tick_params()

    ax1.legend(loc='upper left', markerscale=0.5, fontsize=legend_fontsize)

    for iters, name in enumerate(lims.sumorder):
        if dplcs == 0:
            ukaf = str(np.int(np.nansum(df[name])))
        else:
            ukaf = str(np.round(np.nansum(df[name]), dplcs))

        if iters == 0:
            ax1.bar(range(0, len(edges)), df[name],
                    color=barcolors[iters],
                    edgecolor='k',
                    label=labels[name] + ': {} {}'.format(ukaf, vollbl))

        else:
            ax1.bar(range(0, len(edges)), df[name],
                    bottom=pd.DataFrame(df[lims.sumorder[0:iters]]).sum(axis=1).values,
                    color=barcolors[iters], edgecolor='k',
                    label=labels[name] + ': {} {}'.format(ukaf, vollbl))

    if ylims is not None:
        ax1.set_ylim(ylims)
    else:
        ylims = ax1.get_ylim()
        ymax = ylims[1] + ylims[1] * 0.5
        ax1.set_ylim((0, ymax))

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

    ax1.set_xlabel('elevation [{}]'.format(elevlbl))
    ax1.set_ylabel('{} '.format(vollbl))
    ax1.yaxis.set_label_position("right")
    ax1.yaxis.tick_right()
    ax1.legend(loc=2, fontsize=8)

    fig.tight_layout()
    fig.subplots_adjust(top=0.92, wspace=0.2)

    fig_name = os.path.join(os.path.abspath(figs_path), fig_name)
    snowav.framework.figures.save_fig(fig, fig_name)

    if logger is not None:
        logger.info(' Saved: {}'.format(fig_name))
