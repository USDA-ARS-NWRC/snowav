
from snowav.utils.MidpointNormalize import MidpointNormalize
import numpy as np
import matplotlib
from matplotlib import pyplot as plt
import matplotlib.colors as mcolors
from mpl_toolkits.axes_grid1 import make_axes_locatable
import seaborn as sns
import copy
import cmocean
import matplotlib.patches as mpatches
import pandas as pd
import snowav.framework.figures

def image_change(args, logger):
    '''
    Makes change in SWE figure, with spatial depth on the left and
    bar volume on right.

    Note: command line tool in snow.py calls this function. If changes are made
    here update that calls as necessary.

    Args
    ----------
    args : dict
        dictionary with required inputs, see swi() figure for more information.

    logger : dict
        snowav logger

    '''

    # print to screen for each figure if desired
    if args['print']:
        print('image_change() figure args:\n','omitting masks and image...\n')
        for name in args.keys():
            if name not in ['masks','image']:
                print(name, ': ', args[name])

    # These get used more and we'll just assign here
    masks = args['masks']
    delta_swe = args['image']
    df = args['df']
    plotorder = args['plotorder']
    lims = args['lims']
    edges = args['edges']
    labels = args['labels']
    barcolors = args['barcolors']

    vMin,vMax = np.nanpercentile(delta_swe,[args['percent_min'],args['percent_max']])
    colorsbad = plt.cm.Set1_r(np.linspace(0., 1, 1))
    colors1 = cmocean.cm.matter_r(np.linspace(0., 1, 127))
    colors2 = plt.cm.Blues(np.linspace(0, 1, 128))
    colors = np.vstack((colorsbad,colors1, colors2))
    mymap = mcolors.LinearSegmentedColormap.from_list('my_colormap', colors)

    ixf = delta_swe == 0
    delta_swe[ixf] = -100000
    pmask = masks[plotorder[0]]['mask']
    ixo = pmask == 0
    delta_swe[ixo] = np.nan
    cmap = copy.copy(mymap)
    cmap.set_bad('white',1.)

    sns.set_style('darkgrid')
    sns.set_context("notebook")

    plt.close(6)
    fig,(ax,ax1) = plt.subplots(num=6, figsize=args['figsize'],
                                dpi=args['dpi'], nrows = 1, ncols = 2)
    h = ax.imshow(delta_swe, interpolation='none',
        cmap = cmap, norm=MidpointNormalize(midpoint=0,
                                            vmin = vMin-0.01,vmax=vMax+0.01))

    # Basin boundaries
    for name in masks:
        ax.contour(masks[name]['mask'],cmap = "Greys",linewidths = 1)

    h.axes.get_xaxis().set_ticks([])
    h.axes.get_yaxis().set_ticks([])
    divider = make_axes_locatable(ax)
    cax = divider.append_axes("right", size="5%", pad=0.2)
    cbar = plt.colorbar(h, cax = cax)
    cbar.set_label(r'$\Delta$ SWE [{}]'.format(args['depthlbl']))

    h.axes.set_title(args['title'])

    if args['dplcs'] == 0:
        tlbl = '{}: {} {}'.format(labels[plotorder[0]],
                             str(int(df[plotorder[0]].sum())),
                             args['vollbl'])
    else:
        tlbl = '{}: {} {}'.format(labels[plotorder[0]],
                             str(np.round(df[plotorder[0]].sum(),
                                          args['dplcs'])),args['vollbl'])

    for iters,name in enumerate(lims.sumorder):

        if args['dplcs'] == 0:
            lbl = '{}: {} {}'.format(labels[name],
                                str(int(df[name].sum())),
                                args['vollbl'])
        else:
            lbl = '{}: {} {}'.format(labels[name],
                                str(np.round(df[name].sum(),
                                args['dplcs'])),args['vollbl'])

        if iters == 0:
            ax1.bar(range(0,len(edges)),df[name],
                    color = barcolors[iters],
                    edgecolor = 'k',label = lbl)

        else:
            ax1.bar(range(0,len(edges)),df[name],
                    bottom = pd.DataFrame(df[lims.sumorder[0:iters]]).sum(axis = 1).values,
                    color = barcolors[iters], edgecolor = 'k',label = lbl)

    ax1.xaxis.set_ticks(range(0,len(edges)))
    plt.tight_layout()
    ax1.set_xlim((args['xlims'][0]-0.5,args['xlims'][1]-0.5))

    edges_lbl = []
    for i in range(0,len(edges)):
        edges_lbl.append(str(int(edges[int(i)])))

    ax1.set_xticklabels(str(i) for i in edges_lbl)
    for tick in ax1.get_xticklabels():
        tick.set_rotation(30)

    ax1.set_xlim((args['xlims'][0]-0.5,args['xlims'][1]-0.5))

    ylims = ax1.get_ylim()
    max = ylims[1] + abs(ylims[1] - ylims[0])
    min = ylims[0] - abs(ylims[1] - ylims[0])*0.1
    ax1.set_ylim((min, max))

    ax1.set_ylabel('{}'.format(args['vollbl']))
    ax1.set_xlabel('elevation [{}]'.format(args['elevlbl']))
    ax1.axes.set_title('Change in Volume')
    ax1.yaxis.set_label_position("right")
    ax1.yaxis.tick_right()

    if len(plotorder) > 1:
        ax1.legend(loc=(lims.legx,lims.legy))

    ax1.text(lims.btx,lims.bty,tlbl,horizontalalignment='center',
             transform=ax1.transAxes,fontsize = 10)

    plt.tight_layout()
    fig.subplots_adjust(top=0.88)

    fig_name_short = 'swe_change_'
    fig_name = '{}{}{}.png'.format(args['figs_path'],fig_name_short,args['directory'])
    if logger is not None:
        logger.info(' saving {}'.format(fig_name))
    snowav.framework.figures.save_fig(fig, fig_name)

    if 'show' in args.keys() and args['show']:
        plt.show()

    return fig_name_short
