
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
from snowav import database
from snowav.database.tables import Basins
from snowav.plotting.plotlims import plotlims as plotlims

def image_change(snow, forecast = None):
    '''
    Change in SWE volume.

    '''

    # Make df from database
    delta_swe_byelev = pd.DataFrame(index = snow.edges, columns = snow.plotorder)

    # Two different calculation schemes, depending on WRF forecast
    if forecast is None:
        run_name = snow.run_name
        outputs = snow.outputs
        ixs = snow.ixs
        ixe = snow.ixe
        start_date = snow.start_date
        end_date = snow.end_date
        name_append = snow.name_append
        title = 'Change in SWE \n {} to {}'.format(
                                snow.report_start.date().strftime("%Y-%-m-%-d"),
                                snow.report_date.date().strftime("%Y-%-m-%-d"))

        # Get change in swe during the specified period
        delta_swe = outputs['swe_z'][ixe] - outputs['swe_z'][ixs]

        for bid in snow.plotorder:
            r = database.database.query(snow,
                                        start_date,
                                        end_date,
                                        run_name,
                                        bid,
                                        'swe_vol')

            for elev in snow.edges:
                v = r[(r['elevation'] == str(elev)) & (r['date_time'] == start_date)]
                v2 = r[(r['elevation'] == str(elev)) & (r['date_time'] == end_date)]
                delta_swe_byelev.loc[elev,bid] = np.nansum(v2['value'].values - v['value'].values)

    else:
        run_name = snow.for_run_name
        outputs = snow.for_outputs
        ixs = 0
        ixe = -1
        start_date = snow.for_start_date
        end_date = snow.for_end_date
        name_append = snow.name_append + '_forecast'
        title = 'Forecast Change in SWE \n {} to {}'.format(
                                start_date.date().strftime("%Y-%-m-%-d"),
                                end_date.date().strftime("%Y-%-m-%-d"))

        # Get change in swe during the specified period
        # print(outputs)
        # print(len(outputs['swe_z']), ixe, ixs)
        delta_swe = outputs['swe_z'][ixs] - snow.outputs['swe_z'][ixe]

        for bid in snow.plotorder:
            # this is for the end of the actual run
            r = database.database.query(snow,
                                         snow.start_date,
                                         snow.end_date,
                                         snow.run_name,
                                         bid,
                                         'swe_vol')

            # this is for the end of the WRF run
            r2 = database.database.query(snow,
                                        start_date,
                                        end_date,
                                        run_name,
                                        bid,
                                        'swe_vol')

            for elev in snow.edges:
                v = r[(r['elevation'] == str(elev)) & (r['date_time'] == snow.end_date)]
                v2 = r2[(r2['elevation'] == str(elev)) & (r2['date_time'] == end_date)]
                delta_swe_byelev.loc[elev,bid] = np.nansum(v2['value'].values - v['value'].values)

    delta_swe = np.multiply(delta_swe,snow.depth_factor)

    qMin,qMax = np.nanpercentile(delta_swe,[0.5,99.5])

    # ix = np.logical_and(delta_swe < qMin, delta_swe >= np.nanmin(np.nanmin(delta_swe)))
    # delta_swe[ix] = qMin + qMin*0.2
    vMin,vMax = np.nanpercentile(delta_swe,[0.5,99.5])

    colorsbad = plt.cm.Set1_r(np.linspace(0., 1, 1))
    colors1 = cmocean.cm.matter_r(np.linspace(0., 1, 127))
    colors2 = plt.cm.Blues(np.linspace(0, 1, 128))
    colors = np.vstack((colorsbad,colors1, colors2))
    mymap = mcolors.LinearSegmentedColormap.from_list('my_colormap', colors)

    ixf = delta_swe == 0
    delta_swe[ixf] = -100000
    pmask = snow.masks[snow.plotorder[0]]['mask']
    ixo = pmask == 0
    delta_swe[ixo] = np.nan
    cmap = copy.copy(mymap)
    cmap.set_bad('white',1.)

    # Get basin-specific lims
    lims = plotlims(snow.basin, snow.plotorder)

    sns.set_style('darkgrid')
    sns.set_context("notebook")

    plt.close(6)
    fig,(ax,ax1) = plt.subplots(num=6, figsize=snow.figsize,
                                dpi=snow.dpi, nrows = 1, ncols = 2)
    h = ax.imshow(delta_swe, interpolation='none',
        cmap = cmap, norm=MidpointNormalize(midpoint=0,
                                            vmin = vMin-0.01,vmax=vMax+0.01))

    # if snow.basin == 'LAKES':
    #     ax.set_xlim(snow.imgx)
    #     ax.set_ylim(snow.imgy)

    # Basin boundaries
    for name in snow.masks:
        ax.contour(snow.masks[name]['mask'],cmap = "Greys",linewidths = 1)

    # Do pretty stuff
    h.axes.get_xaxis().set_ticks([])
    h.axes.get_yaxis().set_ticks([])
    divider = make_axes_locatable(ax)
    cax = divider.append_axes("right", size="5%", pad=0.2)
    cbar = plt.colorbar(h, cax = cax)
    cbar.set_label(r'$\Delta$ SWE [{}]'.format(snow.depthlbl))

    h.axes.set_title(title)

    if snow.dplcs == 0:
        tlbl = '{}: {} {}'.format(snow.plotorder[0],
                             str(int(delta_swe_byelev[snow.plotorder[0]].sum())),
                             snow.vollbl)
    else:
        tlbl = '{}: {} {}'.format(snow.plotorder[0],
                             str(np.round(delta_swe_byelev[snow.plotorder[0]].sum(),
                                          snow.dplcs)),snow.vollbl)

    for iters,name in enumerate(lims.sumorder):

        if snow.dplcs == 0:
            lbl = '{}: {} {}'.format(name,
                                str(int(delta_swe_byelev[name].sum())),
                                snow.vollbl)
        else:
            lbl = '{}: {} {}'.format(name,
                                str(np.round(delta_swe_byelev[name].sum(),
                                snow.dplcs)),snow.vollbl)

        if iters == 0:
            ax1.bar(range(0,len(snow.edges)),delta_swe_byelev[name],
                    color = snow.barcolors[iters],
                    edgecolor = 'k',label = lbl)

        else:
            ax1.bar(range(0,len(snow.edges)),delta_swe_byelev[name],
                    bottom = pd.DataFrame(delta_swe_byelev[lims.sumorder[0:iters]]).sum(axis = 1).values,
                    color = snow.barcolors[iters], edgecolor = 'k',label = lbl)


    ax1.xaxis.set_ticks(range(0,len(snow.edges)))
    plt.tight_layout()
    ax1.set_xlim((snow.xlims[0]-0.5,snow.xlims[1]+0.5))

    edges_lbl = []
    for i in range(0,len(snow.edges)):
        edges_lbl.append(str(int(snow.edges[int(i)])))

    ax1.set_xticklabels(str(i) for i in edges_lbl)
    for tick in ax1.get_xticklabels():
        tick.set_rotation(30)

    ax1.set_xlim((snow.xlims[0]-0.5,snow.xlims[1]+0.5))

    ylims = ax1.get_ylim()
    max = ylims[1] + abs(ylims[1] - ylims[0])
    min = ylims[0] - abs(ylims[1] - ylims[0])*0.1
    ax1.set_ylim((min, max))

    ax1.set_ylabel('{} - per elevation band'.format(snow.vollbl))
    ax1.set_xlabel('elevation [{}]'.format(snow.elevlbl))
    ax1.axes.set_title('Change in SWE')

    ax1.yaxis.set_label_position("right")
    ax1.yaxis.tick_right()

    # snow-free
    # patches = [mpatches.Patch(color='grey', label='snow free')]
    # if sum(sum(ixf)) > 1000:
    #     ax.legend(handles=patches, bbox_to_anchor=(lims.pbbx, 0.05),
    #               loc=2, borderaxespad=0.)

    # basin total and legend
    if len(snow.plotorder) > 1:
        ax1.legend(loc=(lims.legx,lims.legy))

    ax1.text(lims.btx,lims.bty,tlbl,horizontalalignment='center',
             transform=ax1.transAxes,fontsize = 10)

    plt.tight_layout()
    fig.subplots_adjust(top=0.88)

    snow._logger.info('saving {}swe_change_{}.png'.format(snow.figs_path,name_append))
    plt.savefig('{}swe_change_{}.png'.format(snow.figs_path,name_append))
