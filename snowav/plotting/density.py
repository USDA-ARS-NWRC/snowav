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

    # Make df from database
    density_byelev = pd.DataFrame(index = snow.edges, columns = snow.plotorder)
    density_summary = pd.DataFrame(columns = snow.plotorder)
    density_summary_8 = pd.DataFrame(columns = snow.plotorder)
    # density_summary_11 = pd.DataFrame(columns = snow.plotorder)
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
        # vd_11 = dr[(dr['elevation'] == '11000')]

        try:
            for iter,d in enumerate(vd['date_time'].values):
                density_summary.loc[d,bid] = vd['value'].values[iter]
                if (vd_8['value'].values[iter] > 0):
                    density_summary_8.loc[d,bid] = vd_8['value'].values[iter]
                #     density_summary_11.loc[d,bid] = vd_11['value'].values[iter]

            for elev in snow.edges:
                v = r[(r['elevation'] == str(elev)) & (r['date_time'] == snow.start_date)]
                v2 = r[(r['elevation'] == str(elev)) & (r['date_time'] == snow.end_date)]

                density_byelev.loc[elev,bid] = v2['value'].values
        except:
            snow._logger.info(' skipping density figure, empty values')
            return

    density_summary.sort_index(inplace=True)
    value = copy.deepcopy(density_byelev)
    lims = plotlims(snow.basin, snow.plotorder)

    sns.set_style('darkgrid')
    sns.set_context("notebook")

    nf = len(snow.masks)
    cvalue = copy.deepcopy(density)

    mymap = plt.cm.get_cmap('BuGn', 5) 
    
    # colors1 = cmocean.cm.speed(np.linspace(0., 1, 255))
    # colors2 = plt.cm.binary(np.linspace(0, 1, 1))
    # colors = np.vstack((colors2, colors1))
    # mymap = mcolors.LinearSegmentedColormap.from_list('my_colormap', colors)

    # This is to get the background white
    pmask = snow.masks[snow.plotorder[0]]['mask']
    ixo = pmask == 0
    cvalue[ixo] = np.nan
    mymap.set_bad('white',1.)

    ixf = cvalue == 0
    cvalue[ixf] = -1
    mymap.set_under('lightgrey',1.)

    plt.close(4)
    fig,(ax,ax1) = plt.subplots(num=4, figsize = snow.figsize,
                                dpi=snow.dpi, nrows = 1, ncols = 2)
    h = ax.imshow(cvalue, vmin = 250, cmap = mymap)

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

    h.axes.set_title('Density\n{}'.format(snow.report_date.date().strftime("%Y-%-m-%-d")))


    sep = 0.05
    wid = 1/len(snow.plotorder)-sep
    # widths = (-wid-sep, 0, wid+sep)
    widths = np.arange((-1 + wid), (1 - wid), wid)

    for i,name in enumerate(lims.sumorder):
        for iter,edge in enumerate(snow.edges):
            if sum(snow.density[name][edge]) > 0:
                bp = ax1.boxplot(snow.density[name][edge],positions=[iter + widths[i]],widths=wid)

                for element in ['boxes', 'whiskers', 'caps']:
                    plt.setp(bp[element], color=snow.barcolors[i])

    ax1.set_xticks(np.arange(0,len(snow.edges)))
    ax1.set_xticklabels([str(x) for x in snow.edges])
    ax1.set_xlim((-0.5, len(snow.edges) - 0.5))

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

    snow._logger.info(' saving {}density_{}.png'.format(snow.figs_path,snow.name_append))
    fig.savefig('{}density_{}.png'.format(snow.figs_path,snow.name_append))


    ###########################################
    #   2nd density figure
    ###########################################

    plt.close(4)
    f = plt.figure(num=4, figsize = snow.figsize, dpi=snow.dpi)
    a = plt.gca()

    for iters,name in enumerate(snow.plotorder):
        if iters == 0:
            lname = '8k ft'
        else:
            lname = '__nolabel__'

        a.plot(density_summary[name],
               color = snow.barcolors[iters],
               label = name)
        a.plot(density_summary_8[name],
               color = snow.barcolors[iters],
               linewidth = 0.75,
               linestyle = ':',
               label = lname)

    if snow.flight_dates is not None:
        for d in snow.flight_dates:
            a.axvline(x=d,linestyle = ':',linewidth = 0.75, color = 'k')

    a.legend(loc=2)
    a.set_ylabel('density [kg/$m^3$]')
    a.set_xlim((snow.start_date, snow.end_date))

    for tick in a.get_xticklabels():
        tick.set_rotation(30)

    f.tight_layout()
    f.savefig('{}density_change_{}.png'.format(snow.figs_path,snow.name_append))
