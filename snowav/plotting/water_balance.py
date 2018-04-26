
import numpy as np
from matplotlib import pyplot as plt
import matplotlib.colors as mcolors
from mpl_toolkits.axes_grid1 import make_axes_locatable
import seaborn as sns
import copy
import cmocean
import matplotlib.patches as mpatches


def water_balance(snow):
    
    cmap = plt.cm.PuBu
    
    precip = snow.precip
    
    sns.set_style('darkgrid')
    sns.set_context("notebook")

    plt.close(12)
    fig,(ax,ax1) = plt.subplots(num=12, 
                                figsize = snow.figsize, 
                                dpi=snow.dpi, 
                                nrows = 1, ncols = 2)      
    h = ax.imshow(precip,cmap = cmap)
    
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
    
    h.axes.set_title('Total Precipitation \n %s to %s'
                     %(snow.dateFrom.date().strftime("%Y-%-m-%-d"),
                       snow.dateTo.date().strftime("%Y-%-m-%-d")))  
          
    # Plot the bars
    # for iters,name in enumerate(snow.plotorder):
    iters = 0
    name = snow.plotorder[iters]
    
    ax1.plot(snow.precip_summary[name],label = 'precipitation') 
    
    ax1.plot(snow.accum_summary[name] 
             + snow.state_summary[name] 
             - snow.evap_summary[name],label = 'SWI + SWE - evap')
       
    ax1.legend()
    
    for tick in ax1.get_xticklabels():
        tick.set_rotation(30)    
        
    plt.tight_layout()  
    fig.subplots_adjust(top=0.88)  
       
    print('saving figure to %swater_balance%s.png'%(snow.figs_path,snow.name_append))   
    plt.savefig('%swater_balance%s.png'%(snow.figs_path,snow.name_append)) 
    
                  