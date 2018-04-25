from snowav.methods.MidpointNormalize import MidpointNormalize
import numpy as np
from matplotlib import pyplot as plt
import matplotlib.colors as mcolors
from mpl_toolkits.axes_grid1 import make_axes_locatable
import seaborn as sns
import copy
import cmocean
import matplotlib.patches as mpatches

def image_change(snow): 
    '''
    This plots snow.delta_state
    
    Should add in more functionality for plotting differences between various images/runs
    
    '''  
    
    # if len(args) != 0:
    #     snow.name_append = snow.name_append + args[0]
          
    # Make copy so that we can add nans for the plots, but not mess up the original
    delta_state = copy.deepcopy(snow.delta_state)
    qMin,qMax = np.percentile(delta_state,[1,99.5])
    
    ix = np.logical_and(delta_state < qMin, delta_state >= np.nanmin(np.nanmin(delta_state)))
    delta_state[ix] = qMin + qMin*0.2
    vMin,vMax = np.percentile(delta_state,[1,99])
    # clims = (qMin,qMax )
       
    # Override if absolute limits are provide in the config
    if hasattr(snow,'ch_clminabs') and hasattr(snow,'ch_clmaxabs'):
        clims       = (snow.ch_clminabs,snow.ch_clmaxabs)
 
    colorsbad = plt.cm.Accent_r(np.linspace(0., 1, 1))           
    colors1 = cmocean.cm.matter_r(np.linspace(0., 1, 127))
    colors2 = plt.cm.Blues(np.linspace(0, 1, 128))
    colors = np.vstack((colorsbad,colors1, colors2))
    mymap = mcolors.LinearSegmentedColormap.from_list('my_colormap', colors)

    ixf = delta_state == 0
    delta_state[ixf] = -100000 # set snow-free  
    pmask = snow.masks[snow.total_lbl]['mask']
    ixo = pmask == 0
    delta_state[ixo] = np.nan
    cmap = copy.copy(mymap)
    cmap.set_bad('white',1.)   
         
    sns.set_style('darkgrid')
    sns.set_context("notebook")
   
    plt.close(6) 
    fig,(ax,ax1) = plt.subplots(num=6, figsize=snow.figsize, 
                                dpi=snow.dpi, nrows = 1, ncols = 2)
    h = ax.imshow(delta_state, interpolation='none', 
        cmap = cmap, norm=MidpointNormalize(midpoint=0,
                                            vmin = vMin-0.01,vmax=vMax+0.01))

    if snow.basin == 'LAKES':    
        ax.set_xlim(snow.imgx)
        ax.set_ylim(snow.imgy)

    # Basin boundaries
    for name in snow.masks:
        ax.contour(snow.masks[name]['mask'],cmap = "Greys",linewidths = 1)

    if snow.basin == 'SJ':
        fix1 = np.arange(1275,1377)
        fix2 = np.arange(1555,1618)
        ax.plot(fix1*0,fix1,'k')
        ax.plot(fix2*0,fix2,'k')    
    
    # Do pretty stuff
    h.axes.get_xaxis().set_ticks([])
    h.axes.get_yaxis().set_ticks([])
    divider = make_axes_locatable(ax)
    cax = divider.append_axes("right", size="5%", pad=0.2)
    cbar = plt.colorbar(h, cax = cax)     
    # cbar.ax.tick_params() 
    
    if snow.units == 'KAF': 
        cbar.set_label(r'$\Delta$ SWE [in]')           
    if snow.units == 'SI':
        cbar.set_label(r'$\Delta$ SWE [mm]')   
        
    h.axes.set_title('Change in SWE \n %s to %s'
                     %(snow.dateFrom.date().strftime("%Y-%-m-%-d"),
                       snow.dateTo.date().strftime("%Y-%-m-%-d")))    
  
    # Plot the bar in order
    sumorder  = snow.plotorder[1:]  
    if snow.basin == 'LAKES' or snow.basin == 'RCEW':
        sumorder = [snow.plotorder[0]]  
    if snow.dplcs == 0:
        tlbl = '%s = %s %s'%(snow.plotorder[0],
                             str(int(snow.delta_state_byelev[snow.plotorder[0]].sum())),
                             snow.vollbl)
    else: 
        tlbl = '%s = %s %s'%(snow.plotorder[0],
                             str(np.round(snow.delta_state_byelev[snow.plotorder[0]].sum(),
                                          snow.dplcs)),snow.vollbl)            
    
    for iters,name in enumerate(sumorder):  
            
        if snow.dplcs == 0:
            lbl = '%s = %s %s'%(name,
                                str(int(snow.delta_state_byelev[name].sum())),
                                snow.vollbl)
        else: 
            lbl = '%s = %s %s'%(name,
                                str(np.round(snow.delta_state_byelev[name].sum(),
                                snow.dplcs)),snow.vollbl)                                 

        if iters == 0:
            ax1.bar(range(0,len(snow.edges)),snow.delta_state_byelev[name],
                    color = snow.barcolors[iters],
                    edgecolor = 'k',label = lbl)
        elif iters == 1:   
            ax1.bar(range(0,len(snow.edges)),snow.delta_state_byelev[name], 
                    bottom = snow.delta_state_byelev[sumorder[iters-1]], 
                    color = snow.barcolors[iters], edgecolor = 'k',label = lbl)
        elif iters == 2:   
            ax1.bar(range(0,len(snow.edges)),snow.delta_state_byelev[name], 
                    bottom = snow.delta_state_byelev[sumorder[iters-1]]
                    + snow.delta_state_byelev[sumorder[iters-2]], 
                    color = snow.barcolors[iters], edgecolor = 'k',label = lbl)
        elif iters == 3:   
            ax1.bar(range(0,len(snow.edges)),snow.delta_state_byelev[name], 
                    bottom = snow.delta_state_byelev[sumorder[iters-1]]
                    + snow.delta_state_byelev[sumorder[iters-2]]
                    + snow.delta_state_byelev[sumorder[iters-3]], 
                    color = snow.barcolors[iters], edgecolor = 'k',label = lbl)
           
    ax1.set_xlim((snow.xlims[0]-0.5,snow.xlims[1]))
    plt.tight_layout()                            
    xts = ax1.get_xticks()
    edges_lbl = []
    for i in xts[0:len(xts)-1]:
        edges_lbl.append(str(int(snow.edges[int(i)])))

    ax1.set_xticklabels(str(i) for i in edges_lbl)
    for tick in ax1.get_xticklabels():
        tick.set_rotation(30) 
         
    if hasattr(snow,"ch_ylims"):
        ax1.set_ylim(snow.ch_ylims)
    else:
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
                     
    if snow.units == 'KAF':
        ax1.set_ylabel('KAF - per elevation band')
        ax1.set_xlabel('elevation [ft]')
        ax1.axes.set_title('Change in SWE')
    if snow.units == 'SI':
        ax1.set_ylabel(r'$km^3$')
        ax1.set_xlabel('elevation [m]')
        ax1.axes.set_title(r'Change in SWE')
           
    ax1.yaxis.set_label_position("right")
    ax1.tick_params(axis='x')
    ax1.tick_params(axis='y')
    ax1.yaxis.tick_right()
    
    patches = [mpatches.Patch(color='grey', label='snow free')]
    if snow.basin == 'SJ':
        ax.legend(handles=patches, bbox_to_anchor=(0.3, 0.05),
                  loc=2, borderaxespad=0. )
    elif snow.basin == 'RCEW':
        ax.legend(handles=patches, bbox_to_anchor=(-0.1, 0.05),
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
        
    if snow.basin == 'BRB':
        ax1.text(0.26,0.96,tlbl,horizontalalignment='center',
                 transform=ax1.transAxes,fontsize = 10)
    
    if snow.basin == 'TUOL' or snow.basin == 'SJ':
        ax1.text(0.3,0.94,tlbl,horizontalalignment='center',
                 transform=ax1.transAxes,fontsize = 10)
  
    plt.tight_layout() 
    fig.subplots_adjust(top=0.88)
    
    print('saving figure to %sswe_change%s.png'%(snow.figs_path,snow.name_append))
    plt.savefig('%sswe_change%s.png'%(snow.figs_path,snow.name_append))    
    