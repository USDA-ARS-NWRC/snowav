
import numpy as np
from matplotlib import pyplot as plt
import matplotlib.colors as mcolors
from mpl_toolkits.axes_grid1 import make_axes_locatable
import seaborn as sns
import copy
import cmocean
import matplotlib.patches as mpatches


def accumulated(snow):
    '''
    Issues:
    - matplotlib, multiple colormaps, and qMin/qMax do not play nice
    
    '''

    if snow.acc_flag == 'False':
        print('Config option, not creating accumulated() fig')
        return
         
    # Only report accum by the subsection between 
    accum = copy.deepcopy(snow.accum_sub) 
    accum_byelev = copy.deepcopy(snow.accum_byelev_sub) 
    
    qMin,qMax = np.percentile(accum,[0,99.8])
    clims = (0,qMax)
    colors1 = cmocean.cm.dense(np.linspace(0., 1, 255))
    colors2 = plt.cm.binary(np.linspace(0, 1, 1))
    colors = np.vstack((colors2, colors1))
    mymap = mcolors.LinearSegmentedColormap.from_list('my_colormap', colors)
     
    sns.set_style('darkgrid')
    sns.set_context("notebook")
    
    # White background
    pmask = snow.masks[snow.total_lbl]['mask']
    ixo = pmask == 0
    accum[ixo] = np.nan
    mymap.set_bad('white',1.) 
    
    # Now set SWI-free to some color

    ixf = accum == 0
    accum[ixf] = -1
    mymap.set_under('grey',1.) 
    
    plt.close(0)
    fig,(ax,ax1) = plt.subplots(num=0, 
                                figsize = snow.figsize, 
                                dpi=snow.dpi, 
                                nrows = 1, ncols = 2)      
    h = ax.imshow(accum, clim = clims, cmap = mymap)
    
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
    # cbar.ax.tick_params() 
    cbar.set_label('[%s]'%(snow.depthlbl))
    
    h.axes.set_title('Accumulated SWI \n %s to %s'
                     %(snow.dateFrom.date().strftime("%Y-%-m-%-d"),
                       snow.dateTo.date().strftime("%Y-%-m-%-d")))  
    
    # Total basin label
    sumorder = snow.plotorder[1:]
    if snow.basin == 'LAKES' or snow.basin == 'RCEW':
        sumorder = [snow.plotorder[0]]
        
    if snow.dplcs == 0:
        tlbl = '%s = %s %s'%(snow.plotorder[0],
                             str(int(accum_byelev[snow.plotorder[0]].sum())),
                             snow.vollbl)
    else: 
        tlbl = '%s = %s %s'%(snow.plotorder[0],
                             str(np.round(accum_byelev[snow.plotorder[0]].sum(),
                            snow.dplcs)),snow.vollbl)
          
    # Plot the bars
    for iters,name in enumerate(sumorder):
        if snow.dplcs == 0:           
            lbl = '%s = %s %s'%(name,str(int(accum_byelev[name].sum())),
                                snow.vollbl)
        else:
            lbl = '%s = %s %s'%(name,str(np.round(accum_byelev[name].sum(),
                                snow.dplcs)),snow.vollbl)
    
        if iters == 0:
            ax1.bar(range(0,len(snow.edges)),accum_byelev[name], 
                    color = snow.barcolors[iters], 
                    edgecolor = 'k',label = lbl)
        elif iters == 1:   
            ax1.bar(range(0,len(snow.edges)),accum_byelev[name], 
                    bottom = accum_byelev[sumorder[iters-1]], 
                    color = snow.barcolors[iters], edgecolor = 'k',label = lbl)
          
        elif iters == 2:   
            ax1.bar(range(0,len(snow.edges)),accum_byelev[name], 
                    bottom = (accum_byelev[sumorder[iters-1]] + accum_byelev[sumorder[iters-2]]), 
                    color = snow.barcolors[iters], edgecolor = 'k',label = lbl)
    
        elif iters == 3:   
            ax1.bar(range(0,len(snow.edges)),accum_byelev[name], 
                    bottom = (accum_byelev[sumorder[iters-1]] + accum_byelev[sumorder[iters-2]] + accum_byelev[sumorder[iters-3]]), 
                    color = snow.barcolors[iters], edgecolor = 'k',label = lbl)
            
           
        plt.rcParams['hatch.linewidth'] = 1
        plt.rcParams['hatch.color'] = 'k'                
        ax1.set_xlim((snow.xlims[0]-0.5,snow.xlims[1]))
    
    plt.tight_layout()
    xts         = ax1.get_xticks()
    edges_lbl   = []
    for i in xts[0:len(xts)-1]:
        edges_lbl.append(str(int(snow.edges[int(i)])))
    
    ax1.set_xticklabels(str(i) for i in edges_lbl)
    for tick in ax1.get_xticklabels():
        tick.set_rotation(30)        
    
    ax1.set_ylabel('%s - per elevation band'%(snow.vollbl))
    ax1.set_xlabel('elevation [%s]'%(snow.elevlbl))
    
    ax1.yaxis.set_label_position("right")
    ax1.yaxis.tick_right()
    
    if hasattr(snow,'acc_ylims'):
        ax1.set_ylim((snow.acc_ylims))
    else:
        ylims = ax1.get_ylim()
        ax1.set_ylim((0,ylims[1] + ylims[1]*0.2))
    
    plt.tight_layout()
    fig.subplots_adjust(top=0.88)
    
    if snow.basin != 'LAKES' and snow.basin != 'RCEW':
        # more ifs for number subs...
        if len(snow.plotorder) == 5:
            ax1.legend(loc= (0.01,0.68))
        elif len(snow.plotorder) == 4:
            ax1.legend(loc= (0.01,0.74))
        
    if snow.basin == 'BRB':
        ax1.text(0.26,0.94,tlbl,horizontalalignment='center',
                 transform=ax1.transAxes,fontsize = 10)
    else:
        ax1.text(0.3,0.94,tlbl,horizontalalignment='center',
                 transform=ax1.transAxes,fontsize = 10)
    
    # Make SWI-free legend if we need one
    if sum(sum(ixf)) > 1000:
        patches = [mpatches.Patch(color='grey', label='no SWI')]
        if snow.basin == 'SJ':
            ax.legend(handles=patches, bbox_to_anchor=(0.3, 0.05), 
                      loc=2, borderaxespad=0. )
        else:
            ax.legend(handles=patches, bbox_to_anchor=(0.05, 0.05), 
                      loc=2, borderaxespad=0. )
         
    # Save if true
    if snow.acc_flag:
        print('saving figure to %sswi%s.png'%(snow.figs_path,snow.name_append))
        plt.savefig('%sswi%s.png'%(snow.figs_path,snow.name_append))  
    
     