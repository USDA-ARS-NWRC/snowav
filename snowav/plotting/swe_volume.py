
import numpy as np
from matplotlib import pyplot as plt
import matplotlib.colors as mcolors
from mpl_toolkits.axes_grid1 import make_axes_locatable
import seaborn as sns
import cmocean
import matplotlib.patches as mpatches
import pandas as pd
import snowav.framework.figures
from pandas.tools.plotting import table
import six

def swe_volume(args, logger = None):
    '''
    Current SWE depth and volume.

    Args
    ----------
    args : dict
        dictionary of values and options for figure. See CoreConfig.ini
            for more information. Use config option [plots]
            print_args_dict: True to print args to the screen if desired.

            df: standard subbasin by elevation dataframe from
                snowav.database.database.collect()
                       Extended Tuolumne  Tuolumne  Cherry Creek  Eleanor
                3000      0.008              0.008         0.000    0.000
                ...
            image: 2D array of SWE
            title: figure title
            directory: figure directory name
            figs_path: base directory for saving figures
            edges: array of elevation bins in min, max, step
            plotorder: list of basins from topo.nc
            labels: dictionary of labels to use with plotorder
            lims: more figure specs from snowav.plotting.plotlims(plotorder)
            masks: masks made from topo.nc and snowav_masks()
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

        logger : dict
            snowav logger

    '''

    # print to screen for each figure if desired
    if args['print']:
        print('swe_volume() figure args:\n','omitting masks and image...\n')
        for name in args.keys():
            if name not in ['masks','image']:
                print(name, ': ', args[name])

    image = args['image']
    masks = args['masks']
    df = args['df']
    plotorder = args['plotorder']
    lims = args['lims']
    edges = args['edges']
    labels = args['labels']
    barcolors = args['barcolors']

    qMin,qMax = np.nanpercentile(image,[0,args['percent_max']])
    clims = (qMin,qMax)
    pmask = masks[plotorder[0]]['mask']
    ixo = pmask == 0
    ixz = image == 0
    image[ixz] = -1

    colorsbad = plt.cm.Set2_r(np.linspace(0., 1, 1))
    colors1 = cmocean.cm.haline_r(np.linspace(0., 1, 254))
    colors = np.vstack((colorsbad,colors1))
    mymap = mcolors.LinearSegmentedColormap.from_list('my_colormap', colors)

    image[ixo] = np.nan
    mymap.set_bad('white',1.)

    sns.set_style('darkgrid')
    sns.set_context("notebook")

    plt.close(21)
    fig,(ax,ax1) = plt.subplots(num=21, figsize=args['figsize'],
                                facecolor = 'white', dpi=args['dpi'],
                                nrows = 1, ncols = 2)

    h = ax.imshow(image, clim=clims, cmap = mymap)

    for name in masks:
        ax.contour(masks[name]['mask'],cmap = "Greys",linewidths = 1)

    h.axes.get_xaxis().set_ticks([])
    h.axes.get_yaxis().set_ticks([])
    divider = make_axes_locatable(ax)
    cax = divider.append_axes("right", size="5%", pad=0.2)
    cbar = plt.colorbar(h, cax = cax)
    cbar.set_label('[{}]'.format(args['depthlbl']))
    h.axes.set_title(args['title'])

    patches = [mpatches.Patch(color='grey', label='snow free')]

    # If there is meaningful snow-free area, include path and label
    if sum(sum(ixz)) > 1000:
        patches = [mpatches.Patch(color='grey', label='snow free')]
        ax.legend(handles=patches, bbox_to_anchor=(lims.pbbx, 0.05),
                  loc=2, borderaxespad=0. )

    ax1.legend(loc=(lims.legx,lims.legy),markerscale = 0.5)

    for iters,name in enumerate(lims.sumorder):

        if args['dplcs'] == 0:
            ukaf = str(np.int(np.nansum(df[name].values)))
        else:
            ukaf = str(np.round(np.nansum(df[name].values),args['dplcs']))

        if iters == 0:
            ax1.bar(range(0,len(edges)),df[name],
                    color = barcolors[iters],
                    edgecolor='k',
                    label=labels[name] + ': {} {}'.format(ukaf,args['vollbl']))

        else:
            ax1.bar(range(0,len(edges)),df[name],
                    bottom = pd.DataFrame(df[lims.sumorder[0:iters]]).sum(axis=1).values,
                    color = barcolors[iters],
                    edgecolor = 'k',
                    label = labels[name]  + ': {} {}'.format(ukaf,args['vollbl']))

    ylims = ax1.get_ylim()
    max = ylims[1] + ylims[1]*0.5
    min = 0
    ax1.set_ylim((min, max))

    # Keep these to use in cold_content()
    swe_ylims = (min,max)

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
    ax1.set_xlabel('elevation [{}]'.format(args['elevlbl']))
    ax1.set_ylabel('{}'.format(args['vollbl']))
    ax1.yaxis.set_label_position("right")
    ax1.yaxis.tick_right()
    ax1.legend(loc = 2, fontsize = 10)

    for tick in ax1.get_xticklabels():
        tick.set_rotation(30)

    fig.tight_layout()
    fig.subplots_adjust(top=0.92,wspace = 0.2)

    fig_name_short = 'swe_volume_'
    fig_name = '{}{}{}.png'.format(args['figs_path'],fig_name_short,
                                   args['directory'])
    if logger is not None:
        logger.info(' saving {}{}{}.png'.format(args['figs_path'],
                                                fig_name,args['directory']))
    snowav.framework.figures.save_fig(fig, fig_name)

    return fig_name_short, swe_ylims

    # SWE table
    # def render_mpl_table(data, col_width=7.0, row_height=0.625, font_size=14,
    #                      header_color='#40466e', row_colors=['#f1f1f2', 'w'], edge_color='w',
    #                      bbox=[0, 0, 1, 1], header_columns=0,
    #                      ax=None, **kwargs):
    #     if ax is None:
    #         size = (np.array(data.shape[::-1]) + np.array([0, 1])) * np.array([col_width, row_height])
    #         fig, ax = plt.subplots(figsize=size)
    #         ax.axis('off')
    #
    #     mpl_table = ax.table(cellText=data.values, bbox=bbox, colLabels=data.columns, **kwargs)
    #
    #     mpl_table.auto_set_font_size(False)
    #     mpl_table.set_fontsize(font_size)
    #
    #     for k, cell in  six.iteritems(mpl_table._cells):
    #         cell.set_edgecolor(edge_color)
    #         if k[0] == 0 or k[1] < header_columns:
    #             cell.set_text_props(weight='bold', color='w')
    #             cell.set_facecolor(header_color)
    #         else:
    #             cell.set_facecolor(row_colors[k[0]%len(row_colors) ])
    #     return fig, ax
    #
    # ft, at = render_mpl_table(swe_total, header_columns=0, col_width=3.5)
    #
    # if (swe_table_path is not None) and (os.path.exists(swe_table_path)):
    #     fig_names = ['{}{}_swe_table.png'.format(swe_table_path,lims.basin_name),
    #                  '{}{}_swe_table.png'.format(figs_path,lims.basin_name)]
    #     figure.save(ft, fig_names)

    # if (day is not None) and day.show:
    #     plt.show()
