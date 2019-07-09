
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

def swe_change(snow):
    '''
    Change in SWE as a depth.

    '''

    delta_swe = snow.outputs['swe_z'][snow.ixe] - snow.outputs['swe_z'][snow.ixs]
    delta_swe = np.multiply(delta_swe,snow.depth_factor)

    # Make df from database
    delta_swe_byelev = pd.DataFrame(index = snow.edges, columns = snow.plotorder)

    for bid in snow.plotorder:
        r = database.database.query(snow, snow.start_date, snow.end_date,
                                    snow.run_name, bid, 'swe_z')

        for elev in snow.edges:
            v = r[(r['elevation'] == str(elev)) & (r['date_time'] == snow.start_date)]
            v2 = r[(r['elevation'] == str(elev)) & (r['date_time'] == snow.end_date)]
            delta_swe_byelev.loc[elev,bid] = np.nansum(v2['value'].values - v['value'].values)

    qMin,qMax = np.nanpercentile(delta_swe,[1,99.9])

    # ix = np.logical_and(delta_swe < qMin, delta_swe >= np.nanmin(np.nanmin(delta_swe)))
    # delta_swe[ix] = qMin + qMin*0.2
    vMin,vMax = np.nanpercentile(delta_swe,[1,99.9])

    colorsbad = plt.cm.Set1_r(np.linspace(0., 1, 1))
    colors1 = cmocean.cm.matter_r(np.linspace(0., 1, 127))
    colors2 = plt.cm.Blues(np.linspace(0, 1, 128))
    colors = np.vstack((colorsbad,colors1, colors2))
    mymap = mcolors.LinearSegmentedColormap.from_list('my_colormap', colors)

    ixf = delta_swe == 0
    delta_swe[ixf] = -100000 # set snow-free
    pmask = snow.masks[snow.plotorder[0]]['mask']
    ixo = pmask == 0
    delta_swe[ixo] = np.nan
    cmap = copy.copy(mymap)
    cmap.set_bad('white',1.)

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

    # if snow.basin == 'SJ':
    #     fix1 = np.arange(1275,1377)
    #     fix2 = np.arange(1555,1618)
    #     ax.plot(fix1*0,fix1,'k')
    #     ax.plot(fix2*0,fix2,'k')

    # Do pretty stuff
    h.axes.get_xaxis().set_ticks([])
    h.axes.get_yaxis().set_ticks([])
    divider = make_axes_locatable(ax)
    cax = divider.append_axes("right", size="5%", pad=0.2)
    cbar = plt.colorbar(h, cax = cax)
    cbar.set_label(r'$\Delta$ SWE [%s]'%(snow.depthlbl))

    h.axes.set_title('Change in SWE \n %s to %s'
                     %(snow.start_date.date().strftime("%Y-%-m-%-d"),
                       snow.report_date.date().strftime("%Y-%-m-%-d")))

    # Total basin label
    if len(snow.plotorder) == 1:
        sumorder = snow.plotorder
        swid = 0.45
        wid = np.linspace(-0.3,0.3,len(sumorder))
    elif len(snow.plotorder) <= 4:
        sumorder = snow.plotorder[1::]
        swid = 0.25
        wid = np.linspace(-0.3,0.3,len(sumorder))
    elif len(snow.plotorder) >= 5:
        sumorder = snow.plotorder[1::]
        swid = 0.1
        wid = np.linspace(-0.4,0.4,len(sumorder))

    wid = np.linspace(-0.25,0.25,len(sumorder))

    for iters,name in enumerate(sumorder):
        lbl = name
        ax1.bar(range(0,len(snow.edges))-wid[iters],
                delta_swe_byelev[name],
                color = snow.barcolors[iters], width = swid,
                                                edgecolor = 'k',
                                                label = lbl)

    plt.tight_layout()
    xts = ax1.get_xticks()
    edges_lbl = []
    for i in xts[0:len(xts)-1]:
        edges_lbl.append(str(int(snow.edges[int(i)])))

    ax1.set_xticklabels(str(i) for i in edges_lbl)
    for tick in ax1.get_xticklabels():
        tick.set_rotation(30)

    ax1.set_xlim((snow.xlims[0]-0.5,snow.xlims[1]+0.5))

    ylims = ax1.get_ylim()
    if ylims[0] < 0 and ylims[1] == 0:
        ax1.set_ylim((ylims[0]+(ylims[0]*0.3),ylims[1]+ylims[1]*0.3))
    if ylims[0] < 0 and ylims[1] > 0:
        ax1.set_ylim((ylims[0]+(ylims[0]*0.3),(ylims[1] + ylims[1]*0.9)))
        if (ylims[1] + ylims[1]*0.9) < abs(ylims[0]):
            ax1.set_ylim((ylims[0]+(ylims[0]*0.3),(-(ylims[0]*0.6))))

        if ylims[1] == 0:
            # ax1.set_ylim((ylims[0]+(ylims[0]*0.3),(-ylims[0])*0.5))
            ax1.set_ylim((ylims[0]+(ylims[0]*0.3),(-ylims[0])*0.65))
        if ylims[0] == 0:
            ax1.set_ylim((ylims[0]+(ylims[0]*0.3),ylims[1]+ylims[1]*0.3))

    ax1.set_ylabel(r'$\Delta$ [{}] - by elevation band'.format(snow.depthlbl))
    ax1.set_xlabel('elevation [{}]'.format(snow.elevlbl))
    ax1.axes.set_title('Change in SWE')

    ax1.yaxis.set_label_position("right")
    ax1.tick_params(axis='x')
    ax1.tick_params(axis='y')
    ax1.yaxis.tick_right()

    patches = [mpatches.Patch(color='grey', label='snow free')]
    if snow.basin == 'SJ':
        ax.legend(handles=patches, bbox_to_anchor=(0.3, 0.05),
                  loc=2, borderaxespad=0. )
    elif snow.basin == 'RCEW':
        ax.legend(handles=patches, bbox_to_anchor=(-0.3, 0.05),
                  loc=2, borderaxespad=0. )
    else:
        ax.legend(handles=patches, bbox_to_anchor=(0.05, 0.05),
                  loc=2, borderaxespad=0. )

    if snow.basin != 'LAKES' and snow.basin != 'RCEW':
        # more ifs for number subs...
        if len(snow.plotorder) == 5:
            ax1.legend(loc= (0.01,0.68))
        elif len(snow.plotorder) == 4:
            ax1.legend(loc= (0.01,0.76))

    plt.tight_layout()
    fig.subplots_adjust(top=0.88)

    snow._logger.info('saving {}swe_change_depth{}.png'.format(snow.figs_path,snow.name_append))
    plt.savefig('{}swe_change_depth{}.png'.format(snow.figs_path,snow.name_append))
