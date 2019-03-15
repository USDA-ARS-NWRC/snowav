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

def density(snow):
    '''

    '''

    density = snow.outputs['density'][snow.ixe]

    # Make df from database
    density_byelev = pd.DataFrame(index = snow.edges, columns = snow.plotorder)
    # delta_density_byelev = pd.DataFrame(index = snow.edges, columns = snow.plotorder)

    for bid in snow.plotorder:
        r = database.database.query(snow, snow.start_date, snow.end_date,
                                    snow.run_name, bid, 'density')

        for elev in snow.edges:
            v = r[(r['elevation'] == str(elev)) & (r['date_time'] == snow.start_date)]
            v2 = r[(r['elevation'] == str(elev)) & (r['date_time'] == snow.end_date)]

            # delta_density_byelev.loc[elev,bid] = np.nansum(v2['value'].values - v['value'].values)
            density_byelev.loc[elev,bid] = v2['value'].values

    value = copy.deepcopy(density_byelev)
    lim = np.nanmax(value[snow.plotorder[0]])
    ylim = (0,600)
    # color = 'xkcd:windows blue'

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

    plt.close(4)
    fig,(ax,ax1) = plt.subplots(num=4, figsize = snow.figsize,
                                dpi=snow.dpi, nrows = 1, ncols = 2)
    h = ax.imshow(cvalue, clim = (50,550), interpolation='none', cmap = mymap)

    # Basin boundaries
    for name in snow.masks:
        ax.contour(snow.masks[name]['mask'],cmap = 'Greys',linewidths = 1)

    if snow.basin == 'LAKES':
        ax.set_xlim(snow.imgx)
        ax.set_ylim(snow.imgy)

    # Do pretty stuff
    h.axes.get_xaxis().set_ticks([])
    h.axes.get_yaxis().set_ticks([])
    divider = make_axes_locatable(ax)

    cax = divider.append_axes("right", size="4%", pad=0.2)
    cbar = plt.colorbar(h, cax = cax)
    cbar.ax.tick_params()
    cbar.set_label('[kg/$m^3$]')

    h.axes.set_title('Density\n%s'%(snow.report_date.date().strftime("%Y-%-m-%-d")))

    ax1.bar(range(0,len(snow.edges)),value[snow.plotorder[0]],
            color = 'g', edgecolor = 'k')

    plt.rcParams['hatch.linewidth'] = 1
    plt.rcParams['hatch.color'] = 'k'

    ax1.set_xlim((snow.xlims[0],snow.xlims[1]))
    xts = ax1.get_xticks()
    edges_lbl = []
    for i in xts[0:len(xts)-1]:
        edges_lbl.append(str(int(snow.edges[int(i)])))

    ax1.set_xticklabels(str(i) for i in edges_lbl)
    for tick in ax1.get_xticklabels():
        tick.set_rotation(30)

    # ax1.set_xlim((snow.xlims[0]-0.5,snow.xlims[1] + 0.5))

    ax1.set_ylabel('density - per elevation band')
    ax1.set_xlabel('elevation [%s]'%(snow.elevlbl))

    ax1.yaxis.set_label_position("right")
    ax1.yaxis.tick_right()
    ax1.set_ylim((0,600))

    plt.tight_layout()
    fig.subplots_adjust(top=0.88)

    if sum(sum(ixf)) > 1000:
        patches = [mpatches.Patch(color='grey', label='snow free')]
        if snow.basin == 'SJ':
            ax.legend(handles=patches, bbox_to_anchor=(0.3, 0.05), loc=2, borderaxespad=0. )
        else:
            ax.legend(handles=patches, bbox_to_anchor=(0.05, 0.05), loc=2, borderaxespad=0. )

    # plt.tight_layout()
    snow._logger.info('saving {}density_{}.png'.format(snow.figs_path,snow.name_append))
    plt.savefig('{}density_{}.png'.format(snow.figs_path,snow.name_append))
