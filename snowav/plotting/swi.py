
import numpy as np
from matplotlib import pyplot as plt
import matplotlib.colors as mcolors
from mpl_toolkits.axes_grid1 import make_axes_locatable
import seaborn as sns
import cmocean
import matplotlib.patches as mpatches
import pandas as pd
import snowav.framework.figures

def swi(args, logger=None):
    '''
    SWI figure, with depth image on the left and volume by elevation band on
    the right.

    All other snowav figures reference this function for args description.

    Args
    ----------
    args : dict
        dictionary of values and options for figure. See CoreConfig.ini
            for more information.

            df: standard subbasin by elevation dataframe from
                snowav.database.database.collect()
                       Extended Tuolumne  Tuolumne  Cherry Creek  Eleanor
                3000      0.008              0.008         0.000    0.000
                ...
            image: 2D array of accumulated SWI
            title: figure title
            directory: figure directory name
            figs_path: base directory for saving figures
            edges: array of elevation bins in min, max, step
            plotorder: list of basins from topo.nc
            labels: dictionary of labels to use with plotorder
            lims: more figure specs from snowav.plotting.plotlims(plotorder)
            masks: masks made from topo.nc and masks()
            figsize: figure size
            dpi: figure dpi
            depthlbl: depth label
            vollbl: volume label
            dplcs: decimal places to round outputs
            barcolors: list of colors for bar plots
            xlims: xlims
            elevlbl: elevation label
            depth_clip: lower limit on image depths for plotting
            percent_min: quantile min for image
            percent_max: quantile max for image

    logger : list
        snowav logger

    '''

    # print to screen for each figure if desired
    if args['print']:
        print('swi() figure args:\n','omitting masks and image...\n')
        for name in args.keys():
            if name not in ['masks','image']:
                print(name, ': ', args[name])

    # These get used more and we'll just assign here
    masks = args['masks']
    accum = args['image']
    df = args['df']
    plotorder = args['plotorder']
    lims = args['lims']
    edges = args['edges']
    labels = args['labels']
    barcolors = args['barcolors']

    qMin,qMax = np.nanpercentile(accum,[0,args['percent_max']])
    clims = (0,qMax)
    colors1 = cmocean.cm.dense(np.linspace(0, 1, 255))
    colors2 = plt.cm.Set1_r(np.linspace(0, 1, 1))
    colors = np.vstack((colors2, colors1))
    mymap = mcolors.LinearSegmentedColormap.from_list('my_colormap', colors)

    pmask = masks[plotorder[0]]['mask']
    ixo = pmask == 0
    accum[ixo] = np.nan
    mymap.set_bad('white',1.)

    # Now set SWI-free to some color
    r = ~np.isnan(accum)
    r[r] &= accum[r] < 0.001
    accum[r] = 0

    sns.set_style('darkgrid')
    sns.set_context("notebook")

    plt.close(0)
    fig,(ax,ax1) = plt.subplots(num=0,
                                figsize=args['figsize'],
                                dpi=args['dpi'],
                                nrows=1,
                                ncols=2)

    h = ax.imshow(accum, cmap=mymap, clim=clims)

    for name in masks:
        ax.contour(masks[name]['mask'], cmap='Greys',linewidths=1)

    h.axes.get_xaxis().set_ticks([])
    h.axes.get_yaxis().set_ticks([])
    divider = make_axes_locatable(ax)
    cax = divider.append_axes("right", size="4%", pad=0.2)
    cbar = plt.colorbar(h, cax = cax)
    cbar.set_label('[{}]'.format(args['depthlbl']))
    h.axes.set_title(args['title'])

    # Plot the bars
    for iters,name in enumerate(lims.sumorder):
        if df[plotorder[0]].sum() < 9.9:
            lbl = '{}: {} {}'.format(labels[name],str(int(df[name].sum())),
                                args['vollbl'])
        else:
            lbl = '{}: {} {}'.format(labels[name],str(np.round(df[name].sum(),
                                args['dplcs'])),args['vollbl'])

        if iters == 0:
            ax1.bar(range(0,len(edges)),df[name],color=args['barcolors'][iters],
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

    ax1.set_ylabel('{}'.format(args['vollbl']))
    ax1.set_xlabel('elevation [{}]'.format(args['elevlbl']))
    ax1.yaxis.set_label_position("right")
    ax1.yaxis.tick_right()

    ylims = ax1.get_ylim()
    max = ylims[1] + ylims[1]*0.6
    ax1.set_ylim((0, max))

    fig.subplots_adjust(top=0.88)

    # If there is meaningful snow-free area, include patch and label
    if sum(sum(r)) > 1000:
        patches = [mpatches.Patch(color='grey', label='no SWI')]
        ax.legend(handles=patches, bbox_to_anchor=(lims.pbbx, 0.05),
                  loc=2, borderaxespad=0. )

    plt.tight_layout()

    # fix so that legend isn't obscured
    if len(plotorder) > 1:
        ax1.legend(loc=2,markerscale = 0.5)

    fig_name_short = 'swi_'
    fig_name = '{}{}{}.png'.format(args['figs_path'],fig_name_short,args['directory'])
    if logger is not None:
        logger.info(' saving {}'.format(fig_name))
    snowav.framework.figures.save_fig(fig, fig_name)

    return fig_name_short
