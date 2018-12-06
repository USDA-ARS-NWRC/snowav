
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
import colorcet as cc
import pandas as pd
from snowav import database
from snowav.database.tables import Basins

def density(snow):
    '''

    '''

    # Calculate accumulated swi during the specified period
    density = snow.outputs['density'][snow.ixe]

    # Make df from database
    density_byelev = pd.DataFrame(index = snow.edges, columns = snow.plotorder)
    delta_density_byelev = pd.DataFrame(index = snow.edges, columns = snow.plotorder)

    for bid in snow.plotorder:
        r = database.database.query(snow, snow.start_date, snow.end_date,
                                    snow.run_name, bid, 'density')

        for elev in snow.edges:
            v = r[(r['elevation'] == str(elev)) & (r['date_time'] == snow.start_date)]
            v2 = r[(r['elevation'] == str(elev)) & (r['date_time'] == snow.end_date)]

            delta_density_byelev.loc[elev,bid] = np.nansum(v2['value'].values - v['value'].values)
            density_byelev.loc[elev,bid] = v2['value'].values

    value = copy.deepcopy(density_byelev)
    lim = np.nanmax(value[snow.plotorder[0]])
    ylim = (0,600)
    color = 'xkcd:windows blue'

    sns.set_style('darkgrid')
    sns.set_context("notebook")

    nf = len(snow.masks)

    if nf == 1:
        nr = 1
        nc = 1

    if nf > 1 and nf < 5:
        nr = 2
        nc = 2

    elif nf == 5 or nf == 6:
        nr = 2
        nc = 3

    elif nf > 6:
        nr = 3
        nc = 3


    plt.close(3)
    fig, axs  = plt.subplots(num=3, figsize=snow.figsize, dpi=snow.dpi,
                                nrows = nr, ncols = nc)

    axs = np.array(axs)
    axs = axs.ravel()

    if nf == 5:
        fig.delaxes(axs[5])

    for iters,name in enumerate(snow.plotorder):

        axs[iters].bar(range(0,len(snow.edges)),value[name], color = color)

        xts = axs[iters].get_xticks()
        if len(xts) < 6:
            dxt = xts[1] - xts[0]
            xts = np.arange(xts[0],xts[-1] + 1 ,dxt/2)
            axs[iters].set_xticks(xts)

        edges_lbl = []
        for i in xts[0:len(xts)-1]:
            edges_lbl.append(str(int(snow.edges[int(i)])))

        axs[iters].set_xticklabels(str(i) for i in edges_lbl)

        if iters == 0:
            axs[iters].set_ylabel(r'$\rho$ [kg/$m^3$]')

        if iters > nc - 1:
            axs[iters].set_xlabel('elevation [%s]'%(snow.elevlbl))

        # Put yaxis on right
        if iters == nc - 1 or iters == nf -1 :
            axs[iters].yaxis.set_label_position("right")
            axs[iters].yaxis.tick_right()
            axs[iters].set_ylabel(r'$\rho$ [kg/$m^3$]')

        if iters == 1 and nc == 3:
            axs[iters].set_yticklabels([])

        if iters <= nc - 1 and nf != 1:
            axs[iters].set_xticklabels([])

        axs[iters].set_ylim((ylim))
        for tick in axs[iters].get_xticklabels():
            tick.set_rotation(30)

        axs[iters].set_xlim((snow.xlims[0]-0.5,snow.xlims[1]+0.5))

        axs[iters].text(0.5, 0.92, name, horizontalalignment = 'center',
                        transform=axs[iters].transAxes,
                        fontsize = 10)

    # for n in range(0,len(axs)):
    #     axs[n].set_xticks(xts)
    #     # axs[n].set_xlim(snow.xlims)
    #     axs[iters].set_xlim((snow.xlims[0]-0.5,snow.xlims[1]+0.5))

    fig.tight_layout()
    fig.subplots_adjust(top=0.92,wspace = 0.1)
    fig.suptitle(r'Density, %s'%(snow.report_date.date().strftime("%Y-%-m-%-d")) )

    snow._logger.info('saving figure to %sdensity_subs_%s.png'%(snow.figs_path,snow.name_append))
    plt.savefig('%sdensity_subs_%s.png'%(snow.figs_path,snow.name_append))

    ###############################
    # 2nd density fig
    ###############################

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

    if snow.basin == 'SJ':
        fix1 = np.arange(1275,1377)
        fix2 = np.arange(1555,1618)
        ax.plot(fix1*0,fix1,'k')
        ax.plot(fix2*0,fix2,'k')

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

    ax1.set_xlim((snow.xlims[0]-0.5,snow.xlims[1] + 0.5))

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

    plt.tight_layout()
    snow._logger.info('saving figure to %sdensity_%s.png'%(snow.figs_path,snow.name_append))
    plt.savefig('%sdensity_%s.png'%(snow.figs_path,snow.name_append))


    ###############################
    # 2nd B
    ###############################

    dvalue = copy.deepcopy(density) - snow.outputs['density'][snow.ixs]

    zf = snow.outputs['swe_z'][snow.ixe] == 0
    dvalue[zf] = np.nan

    qMin,qMax = np.percentile(dvalue,[0.1,99.9])
    clims = (qMin,qMax)

    colors1 = cmocean.cm.speed(np.linspace(0., 1, 255))
    colors2 = plt.cm.binary(np.linspace(0, 1, 1))
    colors = np.vstack((colors2, colors1))
    mymap = mcolors.LinearSegmentedColormap.from_list('my_colormap', colors)

    sns.set_style('darkgrid')
    sns.set_context("notebook")

    # This is to get the background white
    pmask = snow.masks[snow.plotorder[0]]['mask']
    ixo = pmask == 0
    dvalue[ixo] = np.nan
    mymap.set_bad('white',1.)

    ixf = dvalue == 0
    dvalue[ixf] = -1
    mymap.set_under('lightgrey',1.)

    plt.close(4)
    fig,(ax,ax1) = plt.subplots(num=4, figsize = snow.figsize,
                                dpi=snow.dpi, nrows = 1, ncols = 2)
    h = ax.imshow(dvalue, clim = (0,100),interpolation='none', cmap = mymap)

    # Basin boundaries
    for name in snow.masks:
        ax.contour(snow.masks[name]['mask'],cmap = 'Greys',linewidths = 1)

    if snow.basin == 'SJ':
        fix1 = np.arange(1275,1377)
        fix2 = np.arange(1555,1618)
        ax.plot(fix1*0,fix1,'k')
        ax.plot(fix2*0,fix2,'k')

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

    h.axes.set_title('Change in Density\n{} to {}'.format(snow.start_date.date().strftime("%Y-%-m-%-d"),
                                snow.report_date.date().strftime("%Y-%-m-%-d")))

    h1 = ax1.bar(range(0,len(snow.edges)),delta_density_byelev[snow.plotorder[0]].values,color = 'g', edgecolor = 'k')
    plt.tight_layout()
    ax1.set_xlim((snow.xlims[0]-0.5,snow.xlims[1] + 0.5))

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

    ax1.set_xlim((snow.xlims[0]-0.5,snow.xlims[1] + 0.5))

    ax1.set_ylabel(r'$\Delta$ density [kg/$m^3$]')
    ax1.set_xlabel('elevation [%s]'%(snow.elevlbl))

    ax1.set_title('Change in Density\n{} to {}'.format(snow.start_date.date().strftime("%Y-%-m-%-d"),
                                snow.report_date.date().strftime("%Y-%-m-%-d")))

    ax1.yaxis.set_label_position("right")
    ax1.yaxis.tick_right()

    yl = ax1.get_ylim()
    ax1.set_ylim((yl[0]-yl[0]*0.2,yl[1]+yl[1]*0.2))

    # plt.tight_layout()
    fig.subplots_adjust(top=0.88)

    if sum(sum(ixf)) > 1000:
        patches = [mpatches.Patch(color='grey', label='snow free')]
        if snow.basin == 'SJ':
            ax.legend(handles=patches, bbox_to_anchor=(0.3, 0.05), loc=2, borderaxespad=0. )
        else:
            ax.legend(handles=patches, bbox_to_anchor=(0.05, 0.05), loc=2, borderaxespad=0. )

    plt.tight_layout()
    snow._logger.info('saving figure to %sdensity_b_%s.png'%(snow.figs_path,snow.name_append))
    plt.savefig('%sdensity_b_%s.png'%(snow.figs_path,snow.name_append))

    ###############################
    # 3rd density fig
    ###############################

    nsub = 100
    depth  = copy.deepcopy(snow.outputs['depth'][snow.ixe])
    density = cvalue
    swe = copy.deepcopy(snow.outputs['swe_z'][snow.ixe])

    depths = np.reshape(depth,(1,len(depth[:,0])*len(depth[0,:])))
    densitys = np.reshape(density,(1,len(depth[:,0])*len(density[0,:])))
    swesats = np.reshape(swe,(1,len(swe[:,0])*len(swe[0,:])))
    densitys = densitys[0,:]
    depths = depths[0,:]
    swesats = swesats[0,:]
    depthsub = depths[::nsub]
    densitysub = densitys[::nsub]
    swesub = swesats[::nsub]

    z = snow.dem
    zs = np.reshape(z,(1,len(z[:,0])*len(z[0,:])))
    zs = zs[0,:]
    zsub = zs[::100]

    cm = cc.m_bgy

    plt.close(20)
    fig = plt.figure(num=20, figsize = (6,4), dpi=snow.dpi)
    ax = plt.gca()
    h = ax.scatter(swesub,densitysub,c=zsub,vmin=5000,vmax=12500 ,cmap=cm,s=5)

    ax.set_ylabel(r'$\rho$ [kg/$m^3$]')
    ax.set_xlabel('SWE [in]')

    divider = make_axes_locatable(ax)
    cax = divider.append_axes("right", size="4%", pad=0.2)
    cbar = plt.colorbar(h, cax = cax)
    cbar.ax.tick_params()
    cbar.set_label('elevation [ft]')
    plt.tight_layout()

    snow._logger.info('saving figure to %sdensity_swe_%s.png'%(snow.figs_path,snow.name_append))
    plt.savefig('%sdensity_swe_%s.png'%(snow.figs_path,snow.name_append))
