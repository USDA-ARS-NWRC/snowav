
import numpy as np
from matplotlib import pyplot as plt
import matplotlib.colors as mcolors
from mpl_toolkits.axes_grid1 import make_axes_locatable
import seaborn as sns
import copy
import cmocean
import matplotlib.patches as mpatches


def water_balance(snow):
    
    precip = snow.precip
    precip_byelev = snow.precip_byelev
    
    sns.set_style('darkgrid')
    sns.set_context("notebook")

    plt.close(12)
    fig,(ax,ax1) = plt.subplots(num=12, 
                                figsize = snow.figsize, 
                                dpi=snow.dpi, 
                                nrows = 1, ncols = 2)      
    h = ax.imshow(precip)
    
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
    cbar.set_label('[%s]'%(snow.depthlbl))
    
    h.axes.set_title('Accumulated Precipitation \n %s to %s'
                     %(snow.dateFrom.date().strftime("%Y-%-m-%-d"),
                       snow.dateTo.date().strftime("%Y-%-m-%-d")))  
    
    # Total basin label
    sumorder = snow.plotorder[1:]
    if snow.basin == 'LAKES' or snow.basin == 'RCEW':
        sumorder = [snow.plotorder[0]]
        
    if snow.dplcs == 0:
        tlbl = '%s = %s %s'%(snow.plotorder[0],
                             str(int(precip_byelev[snow.plotorder[0]].sum())),
                             snow.vollbl)
    else: 
        tlbl = '%s = %s %s'%(snow.plotorder[0],
                             str(np.round(precip_byelev[snow.plotorder[0]].sum(),
                            snow.dplcs)),snow.vollbl)
          
    # Plot the bars
    for iters,name in enumerate(sumorder):
        if snow.dplcs == 0:           
            lbl = '%s = %s %s'%(name,str(int(precip_byelev[name].sum())),
                                snow.vollbl)
        else:
            lbl = '%s = %s %s'%(name,str(np.round(precip_byelev[name].sum(),
                                snow.dplcs)),snow.vollbl)
    
        if iters == 0:
            ax1.bar(range(0,len(snow.edges)),precip_byelev[name], 
                    color = snow.barcolors[iters], 
                    edgecolor = 'k',label = lbl)
        elif iters == 1:   
            ax1.bar(range(0,len(snow.edges)),precip_byelev[name], 
                    bottom = precip_byelev[sumorder[iters-1]], 
                    color = snow.barcolors[iters], edgecolor = 'k',label = lbl)
          
        elif iters == 2:   
            ax1.bar(range(0,len(snow.edges)),precip_byelev[name], 
                    bottom = (precip_byelev[sumorder[iters-1]] + precip_byelev[sumorder[iters-2]]), 
                    color = snow.barcolors[iters], edgecolor = 'k',label = lbl)
    
        elif iters == 3:   
            ax1.bar(range(0,len(snow.edges)),precip_byelev[name], 
                    bottom = (precip_byelev[sumorder[iters-1]] + precip_byelev[sumorder[iters-2]] + precip_byelev[sumorder[iters-3]]), 
                    color = snow.barcolors[iters], edgecolor = 'k',label = lbl)
            