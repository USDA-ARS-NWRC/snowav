
import numpy as np
import matplotlib
from matplotlib import pyplot as plt
import matplotlib.colors as mcolors
from mpl_toolkits.axes_grid1 import make_axes_locatable
import seaborn as sns
import copy
import cmocean
import matplotlib.patches as mpatches
from snowav import database
from snowav.database.tables import Basins
import pandas as pd
from snowav.plotting.plotlims import plotlims as plotlims
from datetime import timedelta
import os
from pandas.tools.plotting import table
import six


def swe_volume(snow=None, forecast=None, day=None):
    '''

    '''

    if day is not None:
        outputs = day.outputs
        title = 'SWE \n{}'.format(day.path)
        depth_factor = day.depth_factor
        masks = day.masks
        plotorder = day.plotorder
        dpi = 200
        figsize = (10,5)
        depthlbl = day.depthlbl
        vollbl = day.vollbl
        elevlbl = day.elevlbl
        basin = day.basin
        dplcs = 1
        xlims = day.xlims
        figs_path = day.figs_path
        name_append = 'day'
        edges = day.edges
        barcolors = day.barcolors

        if day.basin == 'LAKES':
            dplcs = 2

    if snow is not None:
        run_name = snow.run_name
        depth_factor = snow.depth_factor
        masks = snow.masks
        plotorder = snow.plotorder
        figsize = snow.figsize
        dpi = snow.dpi
        depthlbl = snow.depthlbl
        vollbl = snow.vollbl
        elevlbl = snow.elevlbl
        basin = snow.basin
        dplcs = 1
        xlims = snow.xlims
        figs_path = snow.figs_path
        edges = snow.edges
        barcolors = snow.barcolors
        run_name = snow.run_name
        outputs = snow.outputs
        start_date = snow.start_date
        end_date = snow.end_date
        name_append = snow.name_append
        title = 'SWE \n {}'.format(snow.report_date.date().strftime("%Y-%-m-%-d"))

        if snow.basin == 'LAKES':
            dplcs = 2

    if forecast is not None:
        run_name = snow.for_run_name
        outputs = snow.for_outputs
        start_date = snow.for_start_date
        end_date = snow.for_end_date
        name_append = snow.name_append + '_forecast'
        title = 'Forecast SWE \n {}'.format(end_date.date().strftime("%Y-%-m-%-d"))

    if len(plotorder) > 1:
        sumorder = plotorder[1::]
    else:
        sumorder = plotorder


    state = copy.deepcopy(outputs['swe_z'][-1])
    state = np.multiply(state, depth_factor)

    qMin,qMax = np.nanpercentile(state,[0,99.8])
    clims = (qMin,qMax)
    pmask = masks[plotorder[0]]['mask']
    ixo = pmask == 0

    # Prepare no-snow and outside of the basin for the colormaps
    ixz = state == 0
    state[ixz] = -1

    # Colormap for snow.state
    colorsbad = plt.cm.Set2_r(np.linspace(0., 1, 1))
    colors1 = cmocean.cm.haline_r(np.linspace(0., 1, 254))
    colors = np.vstack((colorsbad,colors1))
    mymap = mcolors.LinearSegmentedColormap.from_list('my_colormap', colors)

    state[ixo] = np.nan
    mymap.set_bad('white',1.)

    sns.set_style('darkgrid')
    sns.set_context("notebook")

    plt.close(21)
    fig,(ax,ax1) = plt.subplots(num=21, figsize=figsize,
                                facecolor = 'white', dpi=dpi,
                                nrows = 1, ncols = 2)

    h = ax.imshow(state, clim=clims, cmap = mymap)

    # Basin boundaries
    for name in masks:
        ax.contour(masks[name]['mask'],cmap = "Greys",linewidths = 1)

    # Do pretty stuff for the left plot
    h.axes.get_xaxis().set_ticks([])
    h.axes.get_yaxis().set_ticks([])
    divider = make_axes_locatable(ax)
    cax = divider.append_axes("right", size="5%", pad=0.2)
    cbar = plt.colorbar(h, cax = cax)
    cbar.set_label('[{}]'.format(depthlbl))

    h.axes.set_title(title)

    # Get basin-specific lims
    lims = plotlims(basin, plotorder)

    patches = [mpatches.Patch(color='grey', label='snow free')]

    # If there is meaningful snow-free area, include path and label
    if sum(sum(ixz)) > 1000:
        patches = [mpatches.Patch(color='grey', label='snow free')]
        ax.legend(handles=patches, bbox_to_anchor=(lims.pbbx, 0.05),
                  loc=2, borderaxespad=0. )

    # basin total and legend
    ax1.legend(loc=(lims.legx,lims.legy),markerscale = 0.5)

    datestr = snow.report_date.date().strftime("%Y-%-m-%-d")
    swe = pd.DataFrame(index = edges, columns = plotorder)
    swe_total = pd.DataFrame(index = ['Total'] + sumorder, columns = [datestr, 'SWE [TAF]'])

    if snow is not None:

        for bid in plotorder:
            r2 = database.database.query(snow, start_date, end_date,
                                        run_name, bid, 'swe_vol')

            total = r2[(r2['elevation'] == 'total')
                  & (r2['date_time'] == end_date)
                  & (r2['basin_id'] == Basins.basins[bid]['basin_id'])]

            if bid == plotorder[0]:
                swe_total.loc['Total','SWE [TAF]'] = total['value'].values[0]
                swe_total.loc['Total',datestr] = 'Total'
            else:
                swe_total.loc[bid,'SWE [TAF]'] = total['value'].values[0]
                swe_total.loc[bid,datestr] = bid


            for elev in edges:
                v2 = r2[(r2['elevation'] == str(elev))
                      & (r2['date_time'] == end_date)
                      & (r2['basin_id'] == Basins.basins[bid]['basin_id'])]

                swe.loc[elev,bid] = v2['value'].values[0]

    else:
        # day.daily_outputs['swe_vol'].values
        for bid in plotorder:
            for i,elev in enumerate(edges):
                swe.loc[elev,bid] = day.daily_outputs['swe_vol'][bid][i]

    ix = len(barcolors)
    # Total basin label

    for iters,name in enumerate(sumorder):

        if dplcs == 0:
            ukaf = str(np.int(np.nansum(swe[name].values)))
        else:
            ukaf = str(np.round(np.nansum(swe[name].values),dplcs))

        if iters == 0:
            #print(iters,edges,name,ukaf,vollbl)
            ax1.bar(range(0,len(edges)),swe[name],
                    color = barcolors[iters],
                    edgecolor = 'k',label = name + ': {} {}'.format(ukaf,vollbl))

        else:
            ax1.bar(range(0,len(edges)),swe[name],
                    bottom = pd.DataFrame(swe[sumorder[0:iters]]).sum(axis = 1).values,
                    color = barcolors[iters], edgecolor = 'k',label = name + ': {} {}'.format(ukaf,vollbl))

    ylims = ax1.get_ylim()
    max = ylims[1] + ylims[1]*0.5
    min = 0
    ax1.set_ylim((min, max))

    ax1.xaxis.set_ticks(range(0,len(edges)))
    plt.tight_layout()
    ax1.set_xlim((xlims[0]-0.5,xlims[1]+0.5))

    edges_lbl = []
    for i in range(0,len(edges)):
        edges_lbl.append(str(int(edges[int(i)])))

    ax1.set_xticklabels(str(i) for i in edges_lbl)
    for tick in ax1.get_xticklabels():
        tick.set_rotation(30)

    ax1.set_xlim((xlims[0]-0.5,xlims[1]+0.5))

    ax1.set_xlabel('elevation [{}]'.format(elevlbl))
    ax1.set_ylabel('{} - per elevation band'.format(vollbl))
    ax1.yaxis.set_label_position("right")
    ax1.yaxis.tick_right()
    ax1.legend(loc = 2, fontsize = 10)

    for tick in ax1.get_xticklabels():
        tick.set_rotation(30)

    fig.tight_layout()
    fig.subplots_adjust(top=0.92,wspace = 0.2)

    if snow is not None:
        snow._logger.info(' saving {}swe_volume_{}.png'.format(figs_path,name_append))
    fig.savefig('{}swe_volume_{}.png'.format(figs_path,name_append))

    if day is not None:
        plt.show()

    # SWE table
    def render_mpl_table(data, col_width=7.0, row_height=0.625, font_size=14,
                         header_color='#40466e', row_colors=['#f1f1f2', 'w'], edge_color='w',
                         bbox=[0, 0, 1, 1], header_columns=0,
                         ax=None, **kwargs):
        if ax is None:
            size = (np.array(data.shape[::-1]) + np.array([0, 1])) * np.array([col_width, row_height])
            fig, ax = plt.subplots(figsize=size)
            ax.axis('off')

        mpl_table = ax.table(cellText=data.values, bbox=bbox, colLabels=data.columns, **kwargs)

        mpl_table.auto_set_font_size(False)
        mpl_table.set_fontsize(font_size)

        for k, cell in  six.iteritems(mpl_table._cells):
            cell.set_edgecolor(edge_color)
            if k[0] == 0 or k[1] < header_columns:
                cell.set_text_props(weight='bold', color='w')
                cell.set_facecolor(header_color)
            else:
                cell.set_facecolor(row_colors[k[0]%len(row_colors) ])
        return fig, ax

    tf, ta = render_mpl_table(swe_total, header_columns=0, col_width=3.5)

    tf.savefig('{}swe_table_{}.png'.format(figs_path,name_append))
