
from mpl_toolkits.axes_grid1 import make_axes_locatable
import numpy as np
from matplotlib import pyplot as plt
import matplotlib.colors as mcolors
import seaborn as sns
import copy
import cmocean
import matplotlib.patches as mpatches
import pandas as pd
import snowav.framework.figures

def density(args, logger = None):
    '''
    Density figure.

    Args
    ----------
    args : dict
        dictionary with required inputs, see swi() figure for more information.

    logger : dict

    '''

    if args['print']:
        print('density() figure args:\n','omitting masks, image, and density...\n')
        for name in args.keys():
            if name not in ['masks','image', 'density']:
                print(name, ': ', args[name])

    masks = args['masks']
    image = args['image']
    density = args['density']
    plotorder = args['plotorder']
    lims = args['lims']
    edges = args['edges']
    labels = args['labels']
    barcolors = args['barcolors']

    qMin,qMax = np.nanpercentile(image,[5,95])
    clims = (qMin,qMax)

    sns.set_style('darkgrid')
    sns.set_context("notebook")

    cvalue = copy.deepcopy(image)

    mymap = plt.cm.get_cmap('BuGn', 8)

    pmask = masks[plotorder[0]]['mask']
    ixo = pmask == 0
    cvalue[ixo] = np.nan
    mymap.set_bad('white',1.)

    ixf = cvalue == 0
    cvalue[ixf] = -1
    mymap.set_under('lightgrey',1.)

    plt.close(4)
    fig,(ax,ax1) = plt.subplots(num=4, figsize = args['figsize'],
                                dpi=args['dpi'], nrows = 1, ncols = 2)
    h = ax.imshow(cvalue, vmin = 200, vmax = 500, cmap = mymap)

    for name in masks:
        ax.contour(masks[name]['mask'],cmap = 'Greys',linewidths = 1)

    h.axes.get_xaxis().set_ticks([])
    h.axes.get_yaxis().set_ticks([])
    divider = make_axes_locatable(ax)

    cax = divider.append_axes("right", size="4%", pad=0.2)
    cbar = plt.colorbar(h, cax = cax)
    cbar.ax.tick_params()
    cbar.set_label('[kg/$m^3$]')

    h.axes.set_title(args['title'])

    sep = 0.05
    wid = 1/len(plotorder)-sep
    widths = np.arange((-1 + wid), (1 - wid), wid)

    for i,name in enumerate(lims.sumorder):
        for iter,edge in enumerate(edges):
            data = density[name][edge].flatten()
            data = data[~np.isnan(data)]

            if np.nansum(data) > 0:
                bp = ax1.boxplot(data,
                                 positions=[iter + widths[i]],
                                 widths=wid)

                for element in ['boxes', 'whiskers', 'caps']:
                    plt.setp(bp[element], color=barcolors[i])

    ax1.set_xticks(np.arange(0,len(edges)))
    ax1.set_xticklabels([str(x) for x in edges])
    ax1.set_xlim(args['xlims'][0],args['xlims'][1])

    for tick in ax1.get_xticklabels():
        tick.set_rotation(30)

    ax1.set_title('Elevation distribution')
    ax1.set_ylabel(r'density [kg/$m^3$]')

    ax1.yaxis.set_label_position("right")
    ax1.yaxis.tick_right()

    plt.tight_layout()
    fig.subplots_adjust(top=0.88)

    if sum(sum(ixf)) > 1000:
        patches = [mpatches.Patch(color='grey', label='snow free')]
        ax.legend(handles=patches, bbox_to_anchor=(lims.pbbx, 0.05),
                  loc=2, borderaxespad=0. )

    fig_name_short = 'density_'
    fig_name = '{}{}{}.png'.format(args['figs_path'],fig_name_short,args['directory'])
    if logger is not None:
        logger.info(' saving {}'.format(fig_name))

    snowav.framework.figures.save_fig(fig, fig_name)

    return fig_name_short

    ###########################################
    #   2nd density figure
    ###########################################
    #
    # plt.close(4)
    # f = plt.figure(num=4, figsize = snow.figsize, dpi=snow.dpi)
    # a = plt.gca()
    #
    # for iters,name in enumerate(snow.plotorder):
    #     if iters == 0:
    #         lname = '8k ft'
    #     else:
    #         lname = '__nolabel__'
    #
    #     a.plot(density_summary[name],
    #            color = snow.barcolors[iters],
    #            label = name)
    #     a.plot(density_summary_8[name],
    #            color = snow.barcolors[iters],
    #            linewidth = 0.75,
    #            linestyle = ':',
    #            label = lname)
    #
    # if snow.flight_dates is not None:
    #     for d in snow.flight_dates:
    #         a.axvline(x=d,linestyle = ':',linewidth = 0.75, color = 'k')
    #
    # a.legend(loc=2)
    # a.set_ylabel('density [kg/$m^3$]')
    # a.set_xlim((snow.start_date, snow.end_date))
    #
    # for tick in a.get_xticklabels():
    #     tick.set_rotation(30)
    #
    # f.tight_layout()
    # f.savefig('{}density_change_{}.png'.format(snow.figs_path,snow.name_append))
