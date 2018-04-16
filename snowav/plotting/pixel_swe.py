
import numpy as np
from matplotlib import pyplot as plt
import matplotlib.colors as mcolors
from mpl_toolkits.axes_grid1 import make_axes_locatable
import seaborn as sns
import copy
import cmocean
import matplotlib.patches as mpatches


def pixel_swe(snow):
    
    
    plt.close(2) 
    fig,(ax,ax1) = plt.subplots(num=2, figsize=snow.figsize,
                                dpi=snow.dpi, nrows = 1, ncols = 2)

    sumorder = snow.plotorder[1:]
    if snow.basin == 'LAKES' or snow.basin == 'RCEW':
        sumorder = [snow.plotorder[0]]  
        swid = 0.45
    else:
        sumorder = snow.plotorder[1::]
        swid = 0.25
    
    wid = np.linspace(-0.3,0.3,len(sumorder))
    
    for iters,name in enumerate(sumorder): 
        # iters = 0
        # name = sumorder[iters]                             
        ax.bar(range(0,len(snow.edges))-wid[iters],
                snow.depth_mdep_byelev[name], 
                color = snow.barcolors[iters], width = swid, edgecolor = 'k',label = name) 
                    
        ax1.bar(range(0,len(snow.edges))-wid[iters],
                snow.state_mswe_byelev[name], 
                color = snow.barcolors[iters], width = swid, edgecolor = 'k',label = name)  
        
    ax.set_xlim((0,len(snow.edges)))    
    ax1.set_xlim((0,len(snow.edges)))        
    
    plt.tight_layout()                            
    xts         = ax1.get_xticks()
    edges_lbl   = []
    for i in xts[0:len(xts)-1]:
        edges_lbl.append(str(int(snow.edges[int(i)])))

    ax.set_xticklabels(str(i) for i in edges_lbl)
    ax1.set_xticklabels(str(i) for i in edges_lbl)
    for tick,tick1 in zip(ax.get_xticklabels(),ax1.get_xticklabels()):
        tick.set_rotation(30)
        tick1.set_rotation(30)  
         
    ylims = ax1.get_ylim()
    ax1.set_ylim((ylims[0],ylims[1]+ylims[1]*0.2))  
    
    ylims = ax.get_ylim()
    ax.set_ylim((ylims[0],ylims[1]+ylims[1]*0.2))               

    ax.set_ylabel('mean depth [%s]'%(snow.depthlbl))
    ax.set_xlabel('elevation [ft]')
    ax1.set_ylabel('mean SWE [%s]'%(snow.depthlbl))
    ax1.set_xlabel('elevation [ft]')
    ax.legend(loc='upper left')
    ax.grid('on')
    ax1.grid('on')
    
    ax1.yaxis.set_label_position("right")
    ax1.yaxis.tick_right()       
    
    ax.set_title('Mean Depth, %s'%(snow.dateTo.date().strftime("%Y-%-m-%-d")))
    ax1.set_title('Mean SWE, %s'%(snow.dateTo.date().strftime("%Y-%-m-%-d")))
    
    plt.tight_layout()
    print('saving figure to %smean_swe_depth%s.png'%(snow.figs_path,snow.name_append))
    plt.savefig('%smean_swe_depth%s.png'%(snow.figs_path,snow.name_append))   
    