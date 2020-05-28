import copy
import matplotlib.patches as mpatches
from matplotlib import pyplot as plt
from mpl_toolkits.axes_grid1 import make_axes_locatable
import numpy as np
import os
import seaborn as sns

import snowav.framework.figures


def density(masks, image, density, plotorder, lims, edges, barcolors,
            clims_percent, title, xlims, figsize, figs_path, fig_name, dpi=200,
            logger=None):
    """ Density figure.

    masks {dict}: basin lookup dict with masks
    image {arr}: density array
    density {dict}: subbasin and elevation band density from process()
    plotorder {list}: basins
    barcolors {list}: plot colors
    clims_percent {list}: colorlims
    title {str}: plot title
    figsize {list}: figure size
    dplcs {int}: decimal places
    xlims {list}: x lims
    figs_path {str}: figure base path
    fig_name {str}: figure name

    """

    xtick_fontsize = 10
    xtick_rotation = 45
    legend_fontsize = 8

    q_min, q_max = np.nanpercentile(image, clims_percent)
    cvalue = copy.deepcopy(image)

    mymap = plt.cm.get_cmap('BuGn', 8)
    sns.set_style('darkgrid')
    sns.set_context("notebook")

    pmask = masks[plotorder[0]]['mask']
    ixo = pmask == 0
    cvalue[ixo] = np.nan
    mymap.set_bad('white', 1.)

    ixf = cvalue == 0
    cvalue[ixf] = -1
    mymap.set_under('lightgrey', 1.)

    plt.close(4)
    fig, (ax, ax1) = plt.subplots(num=4, figsize=figsize, dpi=dpi, nrows=1,
                                  ncols=2)
    h = ax.imshow(cvalue, vmin=100, vmax=q_max, cmap=mymap)

    for name in masks:
        ax.contour(masks[name]['mask'], cmap='Greys', linewidths=1)

    h.axes.get_xaxis().set_ticks([])
    h.axes.get_yaxis().set_ticks([])
    divider = make_axes_locatable(ax)

    cax = divider.append_axes("right", size="2.5%", pad=0.1)
    cbar = plt.colorbar(h, cax=cax)
    cbar.ax.tick_params()
    cbar.set_label('[kg/$m^3$]')

    h.axes.set_title(title)

    wid = 1 / len(plotorder)
    widths = np.arange((-0.5 + wid), (0.5 - wid), wid)

    if len(widths) == 0:
        widths = [0]

    bp = []
    for i, name in enumerate(lims.sumorder):
        for idx, edge in enumerate(edges):
            data = density[name][edge].flatten()
            data = data[~np.isnan(data)]
            if np.nansum(data) > 0:
                bpo = ax1.boxplot(data,
                                  positions=[idx + widths[i]],
                                  widths=wid)

                for element in ['boxes', 'whiskers', 'caps']:
                    plt.setp(bpo[element], color=barcolors[i])

        bp.append(bpo["boxes"][0])

    ax1.set_xticks(np.arange(0, len(edges)))
    ax1.set_xticklabels([str(x) for x in edges])
    ax1.set_xlim(xlims[0] - 1, xlims[1])
    ax1.tick_params(axis='x', labelsize=xtick_fontsize)
    for tick in ax1.get_xticklabels():
        tick.set_rotation(xtick_rotation)

    ax1.set_ylabel('[kg/$m^3$]')

    ax1.yaxis.set_label_position("right")
    ax1.yaxis.tick_right()
    ax1.legend(bp, lims.sumorder, loc=1, fontsize=legend_fontsize)

    plt.tight_layout()
    fig.subplots_adjust(top=0.88)

    if sum(sum(ixf)) > 1000:
        patches = [mpatches.Patch(color='grey', label='snow free')]
        ax.legend(handles=patches, bbox_to_anchor=(lims.pbbx, 0.05),
                  loc=2, borderaxespad=0.)

    fig_name = os.path.join(os.path.abspath(figs_path), fig_name)
    snowav.framework.figures.save_fig(fig, fig_name)

    if logger is not None:
        logger.info(' Saved: {}'.format(fig_name))
