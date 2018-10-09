
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

def precip_depth(snow):
    '''
    SWI, precipitation, and rain as depths.

    '''

    # Get all images first, so we can set global colorlims

    # SWI
    ixs = np.where(snow.outputs['dates'] == snow.start_date)[0][0]
    ixe = np.where(snow.outputs['dates'] == snow.end_date)[0][0]
    accum = np.zeros((snow.nrows,snow.ncols))

    for n in range(ixs,ixe):
        accum = accum + snow.outputs['swi_z'][n]

    accum = np.multiply(accum, snow.depth_factor)
    precip = np.multiply(snow.precip_total, snow.depth_factor)
    rain = np.multiply(snow.rain_total, snow.depth_factor)

    if np.percentile(accum,[0,99.9])[1] > np.percentile(precip,[0,99.9])[1]:
        qMin,qMax = np.percentile(accum,[0,99.9])
    else:
        qMin,qMax = np.percentile(precip,[0,99.9])

    clims = (0,qMax)

    # Get accum by elevation
    accum_byelev = pd.DataFrame(index = snow.edges, columns = snow.plotorder)

    for bid in snow.plotorder:
        r = database.database.query(snow, snow.start_date, snow.end_date,
                                    snow.run_name, bid, 'swi_z')

        for elev in snow.edges:
            v = r[r['elevation'] == str(elev)]
            accum_byelev.loc[elev,bid] = np.nansum(v['value'].values)

    # Get precip by elevation
    precip_byelev = pd.DataFrame(index = snow.edges, columns = snow.plotorder)

    for bid in snow.plotorder:
        r = database.database.query(snow, snow.start_date, snow.end_date,
                                    snow.run_name, bid, 'precip_z')

        for elev in snow.edges:
            v = r[(r['elevation'] == str(elev)) & (r['date_time'] == snow.start_date)]
            v2 = r[(r['elevation'] == str(elev)) & (r['date_time'] == snow.end_date)]
            precip_byelev.loc[elev,bid] = np.nansum(v2['value'].values - v['value'].values)

    # Get rain by elevation
    rain_byelev = pd.DataFrame(index = snow.edges, columns = snow.plotorder)

    for bid in snow.plotorder:
        r = database.database.query(snow, snow.start_date, snow.end_date,
                                    snow.run_name, bid, 'rain_z')

        for elev in snow.edges:
            v = r[(r['elevation'] == str(elev)) & (r['date_time'] == snow.start_date)]
            v2 = r[(r['elevation'] == str(elev)) & (r['date_time'] == snow.end_date)]
            rain_byelev.loc[elev,bid] = np.nansum(v2['value'].values - v['value'].values)

    # Get bar plot ylims
    if accum_byelev.values.max() > precip_byelev.values.max():
        yMax = accum_byelev.values.max() + accum_byelev.values.max()*0.2
    else:
        yMax = precip_byelev.values.max() + precip_byelev.values.max()*0.2

    cmap = cmocean.cm.dense
    cmap1 = plt.cm.nipy_spectral_r

    sns.set_style('darkgrid')
    sns.set_context("notebook")

    plt.close(0)
    fig, ax = plt.subplots(num=0, figsize = (10,10),
                                dpi=snow.dpi, nrows = 3, ncols = 2)

    ################################################
    #           SWI                                #
    ################################################

    # White background plasma_r cmocean.cm.thermal
    mymap = copy.deepcopy(cmap)
    pmask = snow.masks[snow.plotorder[0]]['mask']
    ixo = pmask == 0
    accum[ixo] = np.nan
    mymap.set_bad('white',1.)

    # Now set SWI-free to some color
    r = ~np.isnan(accum)
    r[r] &= accum[r] < 0.05
    accum[r] = -1
    mymap.set_under('grey',1.)

    h = ax[0,0].imshow(accum, clim = clims, cmap = mymap)

    # Basin boundaries
    for name in snow.masks:
        ax[0,0].contour(snow.masks[name]['mask'],cmap = 'Greys',linewidths = 1)

    if snow.basin == 'SJ':
        fix1 = np.arange(1275,1377)
        fix2 = np.arange(1555,1618)
        ax[0,0].plot(fix1*0,fix1,'k')
        ax[0,0].plot(fix2*0,fix2,'k')

    if snow.basin == 'LAKES':
        ax[0,0].set_xlim(snow.imgx)
        ax[0,0].set_ylim(snow.imgy)

    # Do pretty stuff
    h.axes.get_xaxis().set_ticks([])
    h.axes.get_yaxis().set_ticks([])
    divider = make_axes_locatable(ax[0,0])
    cax = divider.append_axes("right", size="4%", pad=0.2)
    cbar = plt.colorbar(h, cax = cax)
    cbar.set_label('SWI [%s]'%(snow.depthlbl))
    h.axes.set_title('Accumulated SWI')

    # Total basin label
    sumorder = snow.plotorder[1:]

    if snow.basin == 'LAKES' or snow.basin == 'RCEW':
        sumorder = [snow.plotorder[0]]
        swid = 0.45
    else:
        sumorder = snow.plotorder[1::]
        swid = 0.25

    wid = np.linspace(-0.25,0.25,len(sumorder))

    for iters,name in enumerate(sumorder):
        lbl = name
        ax[0,1].bar(range(0,len(snow.edges))-wid[iters],
                accum_byelev[name],
                color = snow.barcolors[iters], width = swid,
                                                edgecolor = 'k',
                                                label = lbl)

    xts = ax[0,1].get_xticks()
    edges_lbl = []
    for i in xts[0:len(xts)-1]:
        edges_lbl.append(str(int(snow.edges[int(i)])))

    ax[0,1].set_xticklabels(str(i) for i in edges_lbl)
    for tick in ax[0,1].get_xticklabels():
        tick.set_rotation(30)

    ax[0,1].set_ylabel('SWI [%s]'%(snow.depthlbl))
    ax[0,1].set_xlim((snow.xlims[0]-0.5,snow.xlims[1]))
    ax[0,1].yaxis.set_label_position("right")
    ax[0,1].yaxis.tick_right()

    ax[0,1].set_ylim((0,yMax))

    if snow.basin != 'LAKES' and snow.basin != 'RCEW':
        # more ifs for number subs...
        if len(snow.plotorder) == 5:
            ax[0,1].legend(loc= (0.01,0.65))
        elif len(snow.plotorder) == 4:
            ax[0,1].legend(loc= (0.01,0.71))

    # Make SWI-free legend if we need one
    if sum(sum(r)) > 1000:
        patches = [mpatches.Patch(color='grey', label='no SWI')]
        if snow.basin == 'SJ':
            ax[0,0].legend(handles=patches, bbox_to_anchor=(0.3, 0.05),
                      loc=2, borderaxespad=0. )
        else:
            ax[0,0].legend(handles=patches, bbox_to_anchor=(0.05, 0.05),
                      loc=2, borderaxespad=0. )

    ################################################
    #           Precip                             #
    ################################################

    mymap1 = copy.deepcopy(cmap1)
    pmask = snow.masks[snow.plotorder[0]]['mask']
    ixo = pmask == 0
    precip[ixo] = np.nan
    mymap1.set_bad('white',1.)
    r = ~np.isnan(precip)
    r[r] &= precip[r] < 0.05
    precip[r] = -1
    mymap1.set_under('grey',1.)

    h2 = ax[1,0].imshow(precip, interpolation='none', cmap = mymap1, clim = clims)

    if snow.basin == 'LAKES':
        ax[1,0].set_xlim(snow.imgx)
        ax[1,0].set_ylim(snow.imgy)

    # Basin boundaries
    for name in snow.masks:
        ax[1,0].contour(snow.masks[name]['mask'],cmap = "Greys",linewidths = 1)

    if snow.basin == 'SJ':
        fix1 = np.arange(1275,1377)
        fix2 = np.arange(1555,1618)
        ax[1,0].plot(fix1*0,fix1,'k')
        ax[1,0].plot(fix2*0,fix2,'k')

    # Do pretty stuff
    h2.axes.get_xaxis().set_ticks([])
    h2.axes.get_yaxis().set_ticks([])
    divider = make_axes_locatable(ax[1,0])
    cax = divider.append_axes("right", size="5%", pad=0.2)
    cbar = plt.colorbar(h2, cax = cax)
    cbar.set_label(r'precip [%s]'%(snow.depthlbl))

    h2.axes.set_title('Total Precipitation')

    sumorder = snow.plotorder[1:]
    if snow.basin == 'LAKES' or snow.basin == 'RCEW':
        sumorder = [snow.plotorder[0]]
        swid = 0.45
    else:
        sumorder = snow.plotorder[1::]
        swid = 0.25

    wid = np.linspace(-0.25,0.25,len(sumorder))

    for iters,name in enumerate(sumorder):
        lbl = name
        ax[1,1].bar(range(0,len(snow.edges))-wid[iters],
                precip_byelev[name],
                color = snow.barcolors[iters], width = swid,
                                                edgecolor = 'k',
                                                label = lbl)

    xts = ax[1,1].get_xticks()
    edges_lbl = []
    for i in xts[0:len(xts)-1]:
        edges_lbl.append(str(int(snow.edges[int(i)])))

    ax[1,1].set_xticklabels(str(i) for i in edges_lbl)
    for tick in ax[1,1].get_xticklabels():
        tick.set_rotation(30)

    ax[1,1].set_xlim((snow.xlims[0]-0.5,snow.xlims[1]))
    ax[1,1].set_ylim((0,yMax))
    ax[1,1].set_ylabel(r'precip [%s]'%(snow.depthlbl))
    ax[1,1].yaxis.set_label_position("right")
    ax[1,1].tick_params(axis='x')
    ax[1,1].tick_params(axis='y')
    ax[1,1].yaxis.tick_right()

    patches = [mpatches.Patch(color='grey', label='no precip')]

    if sum(sum(r)) > 1000:
        if snow.basin == 'SJ':
            ax[1,0].legend(handles=patches, bbox_to_anchor=(0.3, 0.05),
                      loc=2, borderaxespad=0. )
        elif snow.basin == 'RCEW':
            ax[1,0].legend(handles=patches, bbox_to_anchor=(-0.3, 0.05),
                      loc=2, borderaxespad=0. )
        else:
            ax[1,0].legend(handles=patches, bbox_to_anchor=(0.05, 0.05),
                      loc=2, borderaxespad=0. )

    if snow.basin != 'LAKES' and snow.basin != 'RCEW':
        # more ifs for number subs...
        if len(snow.plotorder) == 5:
            ax[1,1].legend(loc= (0.01,0.65))
        elif len(snow.plotorder) == 4:
            ax[1,1].legend(loc= (0.01,0.71))

    ################################################
    #           rain                                #
    ################################################

    ix = rain_byelev < 0
    rain_byelev[ix] = 0

    mymap = copy.deepcopy(cmap1)
    pmask = snow.masks[snow.plotorder[0]]['mask']
    ixo = pmask == 0
    rain[ixo] = np.nan
    mymap.set_bad('white',1.)
    r = ~np.isnan(rain)
    r[r] &= rain[r] < 0.05
    rain[r] = -1
    mymap.set_under('grey',1.)

    h2 = ax[2,0].imshow(rain, interpolation='none', cmap = mymap, clim = clims)

    if snow.basin == 'LAKES':
        ax[2,0].set_xlim(snow.imgx)
        ax[2,0].set_ylim(snow.imgy)

    # Basin boundaries
    for name in snow.masks:
        ax[2,0].contour(snow.masks[name]['mask'],cmap = "Greys",linewidths = 1)

    if snow.basin == 'SJ':
        fix1 = np.arange(1275,1377)
        fix2 = np.arange(1555,1618)
        ax[2,0].plot(fix1*0,fix1,'k')
        ax[2,0].plot(fix2*0,fix2,'k')

    # Do pretty stuff
    h2.axes.get_xaxis().set_ticks([])
    h2.axes.get_yaxis().set_ticks([])
    divider = make_axes_locatable(ax[2,0])
    cax = divider.append_axes("right", size="5%", pad=0.2)
    cbar = plt.colorbar(h2, cax = cax)
    cbar.set_label(r'rain [%s]'%(snow.depthlbl))

    h2.axes.set_title('Rain')

    sumorder = snow.plotorder[1:]
    if snow.basin == 'LAKES' or snow.basin == 'RCEW':
        sumorder = [snow.plotorder[0]]
        swid = 0.45
    else:
        sumorder = snow.plotorder[1::]
        swid = 0.25

    wid = np.linspace(-0.25,0.25,len(sumorder))

    for iters,name in enumerate(sumorder):
        lbl = name
        ax[2,1].bar(range(0,len(snow.edges))-wid[iters],
                rain_byelev[name],
                color = snow.barcolors[iters], width = swid,
                                                edgecolor = 'k',
                                                label = lbl)

    xts = ax[2,1].get_xticks()
    edges_lbl = []
    for i in xts[0:len(xts)-1]:
        edges_lbl.append(str(int(snow.edges[int(i)])))

    ax[2,1].set_xticklabels(str(i) for i in edges_lbl)
    for tick in ax[2,1].get_xticklabels():
        tick.set_rotation(30)

    ax[2,1].set_xlim((snow.xlims[0]-0.5,snow.xlims[1]))
    ax[2,1].set_ylim((0,yMax))
    ax[2,1].set_ylabel(r'rain [%s]'%(snow.depthlbl))
    ax[2,1].set_xlabel('elevation [%s]'%(snow.elevlbl))
    ax[2,1].yaxis.set_label_position("right")
    ax[2,1].tick_params(axis='x')
    ax[2,1].tick_params(axis='y')
    ax[2,1].yaxis.tick_right()

    patches = [mpatches.Patch(color='grey', label='no rain')]
    if sum(sum(r)) > 1000:
        if snow.basin == 'SJ':
            ax[2,0].legend(handles=patches, bbox_to_anchor=(0.3, 0.05),
                      loc=2, borderaxespad=0. )
        elif snow.basin == 'RCEW':
            ax[2,0].legend(handles=patches, bbox_to_anchor=(-0.3, 0.05),
                      loc=2, borderaxespad=0. )
        else:
            ax[2,0].legend(handles=patches, bbox_to_anchor=(0.05, 0.05),
                      loc=2, borderaxespad=0. )

    if snow.basin != 'LAKES' and snow.basin != 'RCEW':
        # more ifs for number subs...
        if len(snow.plotorder) == 5:
            ax[2,1].legend(loc= (0.01,0.65))
        elif len(snow.plotorder) == 4:
            ax[2,1].legend(loc= (0.01,0.71))

    plt.suptitle('Depth of SWI, Precipitation, and Rain\n%s to %s'
                         %(snow.start_date.date().strftime("%Y-%-m-%-d"),
                           snow.end_date.date().strftime("%Y-%-m-%-d")))
    plt.tight_layout()
    fig.subplots_adjust(top=0.92)

    snow._logger.info('saving figure to %sprecip_depth%s.png'%(snow.figs_path,snow.name_append))
    plt.savefig('%sprecip_depth%s.png'%(snow.figs_path,snow.name_append))
