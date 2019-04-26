from mpl_toolkits.axes_grid1 import make_axes_locatable
import numpy as np
from matplotlib import pyplot as plt
import matplotlib.colors as mcolors
import seaborn as sns
import copy
import cmocean
import matplotlib.patches as mpatches
import pandas as pd
from snowav import database
from snowav.database.tables import Basins
from snowav.plotting.plotlims import plotlims as plotlims
from datetime import datetime


def density(snow):
    '''

    '''

    density = snow.outputs['density'][snow.ixe]
    qMin,qMax = np.nanpercentile(density,[10,90])
    clims = (qMin,qMax)
    # print('density, ', snow.density)

    # Make df from database
    density_byelev = pd.DataFrame(index = snow.edges, columns = snow.plotorder)
    density_summary = pd.DataFrame(columns = snow.plotorder)
    density_summary_8 = pd.DataFrame(columns = snow.plotorder)
    density_summary_11 = pd.DataFrame(columns = snow.plotorder)
    # delta_density_byelev = pd.DataFrame(index = snow.edges, columns = snow.plotorder)

    for bid in snow.plotorder:
        r = database.database.query(snow, snow.start_date, snow.end_date,
                                    snow.run_name, bid, 'density')

        dr = database.database.query(snow,
                                     datetime(snow.wy-1,10,1),
                                     snow.end_date,
                                     snow.run_name,
                                     bid,
                                     'density')

        vd = dr[(dr['elevation'] == 'total')]
        vd_8 = dr[(dr['elevation'] == '8000')]
        vd_11 = dr[(dr['elevation'] == '11000')]

        for iter,d in enumerate(vd['date_time'].values):
            density_summary.loc[d,bid] = vd['value'].values[iter]
            if (vd_8['value'].values[iter] > 0) and (vd_11['value'].values[iter] > 0):
                density_summary_8.loc[d,bid] = vd_8['value'].values[iter]
                density_summary_11.loc[d,bid] = vd_11['value'].values[iter]

        for elev in snow.edges:
            # a = r[(r['elevation'] == str(elev))
            #       & (r['date_time'] >= snow.start_date)
            #       & (r['date_time'] <= snow.start_date)]

            v = r[(r['elevation'] == str(elev)) & (r['date_time'] == snow.start_date)]
            v2 = r[(r['elevation'] == str(elev)) & (r['date_time'] == snow.end_date)]

            # delta_density_byelev.loc[elev,bid] = np.nansum(v2['value'].values - v['value'].values)
            density_byelev.loc[elev,bid] = v2['value'].values

    density_summary.sort_index(inplace=True)

    value = copy.deepcopy(density_byelev)
    # lim = np.nanmax(value[snow.plotorder[0]])
    lims = plotlims(snow.basin, snow.plotorder)

    sns.set_style('darkgrid')
    sns.set_context("notebook")

    nf = len(snow.masks)
    cvalue = copy.deepcopy(density)
    colors1 = cmocean.cm.speed(np.linspace(0., 1, 255))
    colors2 = plt.cm.binary(np.linspace(0, 1, 1))
    colors = np.vstack((colors2, colors1))
    mymap = mcolors.LinearSegmentedColormap.from_list('my_colormap', colors)

    sns.set_style('darkgrid')
    sns.set_context("notebook")

    # This is to get the background white
    pmask = snow.masks[snow.plotorder[0]]['mask']
    ixo = pmask == 0
    cvalue[ixo] = np.nan
    mymap.set_bad('white',1.)

    ixf = cvalue == 0
    cvalue[ixf] = -1
    mymap.set_under('lightgrey',1.)

    # qMin,qMax = np.nanpercentile(density[density > 0],[3,95])
    # clims = (qMin,qMax)
    #
    # # hack density to force colormap for imshow()
    # density[density < qMin] = qMin + 1
    # density[density > qMax] = qMax - 1
    # print('density, ', clims)

    plt.close(4)
    fig,(ax,ax1) = plt.subplots(num=4, figsize = snow.figsize,
                                dpi=snow.dpi, nrows = 1, ncols = 2)
    h = ax.imshow(cvalue, vmin=300, vmax=550, cmap = mymap)

    # Basin boundaries
    for name in snow.masks:
        ax.contour(snow.masks[name]['mask'],cmap = 'Greys',linewidths = 1)

    # Do pretty stuff
    h.axes.get_xaxis().set_ticks([])
    h.axes.get_yaxis().set_ticks([])
    divider = make_axes_locatable(ax)

    cax = divider.append_axes("right", size="4%", pad=0.2)
    cbar = plt.colorbar(h, cax = cax)
    cbar.ax.tick_params()
    cbar.set_label('[kg/$m^3$]')

    h.axes.set_title('Density\n%s'%(snow.report_date.date().strftime("%Y-%-m-%-d")))

    ####################
    #       ax1
    ####################
    # sep = 0.05
    # wid = 1/len(snow.plotorder)-sep
    # widths = (-wid-sep, 0, wid+sep)
    # for i,name in enumerate(lims.sumorder):
    #     for iter,edge in enumerate(snow.edges):
    #         if sum(snow.density[name][edge]) > 0:
    #             bp = ax1.boxplot(snow.density[name][edge],positions=[iter + widths[i]],widths=wid)
    #
    #         for element in ['boxes', 'whiskers', 'caps']:
    #             plt.setp(bp[element], color=snow.barcolors[i])
    #         # for patch in bp['boxes']:
    #         #     patch.set(facecolor='b')
    #
    # # ax1.legend([bp["boxes"][0], bp["boxes"][1]], ['A', 'B'], loc='upper right')
    # ax1.set_xticks(np.arange(0,len(snow.edges)))
    # ax1.set_xticklabels([str(x) for x in snow.edges])
    # ax1.set_xlim((-0.5, len(snow.edges) - 0.5))
    #
    # for tick in ax1.get_xticklabels():
    #     tick.set_rotation(30)
    #
    # ax1.set_title('Elevation distribution')
    # ax1.set_ylabel(r'density [kg/$m^3$]')

    # print('density, ', density_summary)
    # print('density, ', density_summary['Lakes Basin'])

    for iters, name in enumerate(lims.sumorder):
        density_summary[name].plot(ax=ax1, color = snow.barcolors[iters],
                                   label = name,
                                   linewidth = 0.75)
        density_summary_8[name].plot(ax=ax1, color = snow.barcolors[iters],
                                     linestyle = ':',
                                     linewidth = 0.5,
                                     label='__nolabel__')
        density_summary_11[name].plot(ax=ax1, color = snow.barcolors[iters],
                                     linestyle = ':',
                                     linewidth = 0.5,
                                     label='__nolabel__')

        density_summary_8.dropna(how='all',inplace=True)
        density_summary_11.dropna(how='all',inplace=True)

        # idx = density_summary_11.index.intersection(density_summary_8.index)
        # density_summary_11 = density_summary_11[idx]
        # idx = density_summary_8.index.intersection(density_summary_11.index)
        # density_summary_8 = density_summary_8[idx]

        ax1.fill_between(density_summary_8.index,
                         density_summary_8[name].values,
                         density_summary_11[name].values)

    ax1.yaxis.set_label_position("right")
    ax1.set_xlim((datetime(snow.wy -1 , 10, 1),snow.end_date))
    ax1.tick_params(axis='y')
    ax1.yaxis.tick_right()
    ax1.set_ylabel('mean density [kg/$m^3$]')
    ax1.legend(loc=2,fontsize=8)

    ax1.yaxis.set_label_position("right")
    ax1.yaxis.tick_right()

    plt.tight_layout()
    fig.subplots_adjust(top=0.88)

    if sum(sum(ixf)) > 1000:
        patches = [mpatches.Patch(color='grey', label='no SWI')]
        ax.legend(handles=patches, bbox_to_anchor=(lims.pbbx, 0.05),
                  loc=2, borderaxespad=0. )

    snow._logger.info(' saving {}density_{}.png'.format(snow.figs_path,snow.name_append))
    plt.savefig('{}density_{}.png'.format(snow.figs_path,snow.name_append))
