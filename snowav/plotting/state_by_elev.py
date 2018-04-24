
import numpy as np
from matplotlib import pyplot as plt
import matplotlib.colors as mcolors
from mpl_toolkits.axes_grid1 import make_axes_locatable
import seaborn as sns
import copy
import cmocean
import matplotlib.patches as mpatches  
    
def state_by_elev(snow): 
    '''
    Plots SWE by elevation, delineated by melt/nonmelt
    
    '''
        
    lim = np.max(snow.melt[snow.total_lbl]) + np.max(snow.nonmelt[snow.total_lbl])
    ylim = np.max(lim) + np.max(lim)*0.3 
    colors = ['xkcd:rose red','xkcd:cool blue']
    fs = list(snow.figsize)
    fs[0] = fs[0]*0.9
    fs = tuple(fs)
    
    sns.set_style('darkgrid')
    sns.set_context("notebook")
    
    nf = len(snow.masks)
    
    if nf > 1 and nf < 5:
        nr = 2
        nc = 2
    elif nf < 7:
        nr = 2
        nc = 3
        
    if nf > 1:
        plt.close(7)
        fig,ax  = plt.subplots(num=7, figsize=fs, dpi=snow.dpi,
                               nrows = nr, ncols = nc)
        axs = ax.ravel()
        if nf == 5:
            fig.delaxes(axs[5])
               
        for iters,name in enumerate(snow.plotorder):
            # iters = 0
            # name = snow.plotorder[iters]
            axs[iters].bar(range(0,len(snow.edges)),snow.melt[name], 
                           color = colors[0], bottom = snow.nonmelt[name])
            axs[iters].bar(range(0,len(snow.edges)),snow.nonmelt[name], 
                           color = colors[1], label = 'unavail ')
            
            axs[iters].set_xlim(snow.xlims)                          
            xts = axs[iters].get_xticks()
            if len(xts) < 6:
                dxt = xts[1] - xts[0]
                xts = np.arange(xts[0],xts[-1] + 1 ,dxt/2)
                axs[iters].set_xticks(xts)                 
    
            edges_lbl = []
            for i in xts[0:len(xts)-1]:
                edges_lbl.append(str(int(snow.edges[int(i)])))
    
            axs[iters].set_xticklabels(str(i) for i in edges_lbl)
            
            if iters > nc - 1:
                if snow.units == 'KAF':
                    axs[iters].set_xlabel('elevation [ft]')
                if snow.units == 'SI':
                    axs[iters].set_xlabel('elevation [m]')                    
                  
            # Put yaxis on right 
            if iters == nc - 1 or iters == nf -1 :
                axs[iters].yaxis.set_label_position("right")
                axs[iters].yaxis.tick_right()
            
            if iters == 1 and nc == 3:
                axs[iters].set_yticklabels([])
                
            if iters <= nc - 1:
                axs[iters].set_xticklabels([])
            
            axs[iters].tick_params(axis='x')
            axs[iters].tick_params(axis='y')
            
            # Get basin total storage in strings for label
            kaf = str(np.int(sum(snow.melt[name])
                             + sum(snow.nonmelt[name]))) 
            
            if snow.dplcs == 0:
                kaf = str(np.int(sum(snow.melt[name])
                                 + sum(snow.nonmelt[name]))) 
            else: 
                kaf = str(np.round(sum(snow.melt[name])
                                   + sum(snow.nonmelt[name]),snow.dplcs))             
    
            if snow.units == 'KAF':          
                axs[iters].text(0.5,0.92,'%s - %s KAF'
                                %(snow.masks[name]['label'],kaf),
                                horizontalalignment='center',
                                transform=axs[iters].transAxes, fontsize = 10)
            if snow.units == 'SI':          
                axs[iters].text(0.5,0.92,'%s - %s $km^3$'
                                %(snow.masks[name]['label'],kaf),
                                horizontalalignment='center',
                                transform=axs[iters].transAxes)                
            
            if iters == 1 and nc == 3:
                axs[iters].set_yticklabels([])
            else:
                axs[iters].set_ylabel(snow.units)
            
            lbl = []
            for n in (0,1):
                if n == 0:
                    
                    if snow.dplcs == 0:
                        kafa = str(np.int(sum(snow.melt[name]))) 
                    else: 
                        kafa = str(np.round(sum(snow.melt[name]),snow.dplcs)) 
    
                    tmpa = (r'avail = %s')%(kafa)                         
                    lbl.append(tmpa)
    
                if snow.dplcs == 0:
                    kafna = str(np.int(sum(snow.nonmelt[name]))) 
                else: 
                    kafna = str(np.round(sum(snow.nonmelt[name]),snow.dplcs))            
    
                tmpna = ('unavail = %s')%(kafna)
                lbl.append(tmpna)                     
    
            axs[iters].legend(lbl, loc = (0.025, 0.65),fontsize = 9)
            
            if hasattr(snow,'elv_ylims'):
                axs[iters].set_ylim((snow.elv_ylims))
            else:
                axs[iters].set_ylim((0,ylim))          
            
            for tick in axs[iters].get_xticklabels():
                tick.set_rotation(30) 
        
        fig.tight_layout()
        for n in range(0,len(axs)):
            axs[n].set_xticks(xts)
            axs[n].set_xlim(snow.xlims)       
    
    else:
        name = snow.plotorder[0]
        plt.close(1)
        fig,ax  = plt.subplots(num=1, figsize=fs, dpi=snow.dpi)
               
        ax.bar(range(0,len(snow.edges)),snow.melt[name],
               color = colors[0], bottom = snow.nonmelt[name])
        ax.bar(range(0,len(snow.edges)),snow.nonmelt[name],
               color = colors[1], label = 'unavail ')
        
        ax.set_xlim(snow.xlims)                          
        xts = ax.get_xticks()
        edges_lbl = []
        for i in xts[0:len(xts)-1]:
            edges_lbl.append(str(int(snow.edges[int(i)])))
    
        ax.set_xticklabels(str(i) for i in edges_lbl)
    
        if snow.units == 'KAF':
            ax.set_xlabel('elevation [ft]')
        if snow.units == 'SI':
            ax.set_xlabel('elevation [m]')                    
        
        ax.tick_params(axis='x')
        ax.tick_params(axis='y')
        
        # Get basin total storage in strings for label
        kaf = str(np.int(sum(snow.melt[name]) + sum(snow.nonmelt[name]))) 
        
        if snow.dplcs == 0:
            kaf = str(np.int(sum(snow.melt[name]) 
                             + sum(snow.nonmelt[name]))) 
        else: 
            kaf = str(np.round(sum(snow.melt[name]) 
                               + sum(snow.nonmelt[name]),snow.dplcs))             
    
        if snow.units == 'KAF':          
            ax.text(.25,0.95,'%s - %s KAF'
                    %(snow.masks[name]['label'],kaf),
                    horizontalalignment='center',transform=ax.transAxes)
            ax.set_ylabel('KAF')
        if snow.units == 'SI':          
            ax.text(-0.88,2.1,'%s - %s $km^3$'
                    %(snow.masks[name]['label'],kaf),
                    horizontalalignment='center',transform=ax.transAxes)
            ax.set_ylabel(r'$km^3$')                
        
        lbl = []
        for n in (0,1):
            if n == 0:
                if snow.dplcs == 0:
                    kafa = str(np.int(sum(snow.melt[name]))) 
                else: 
                    kafa = str(np.round(sum(snow.melt[name]),snow.dplcs)) 
                
                tmpa = ('avail = %s')%(kafa)
                                      
                lbl.append(tmpa)
            
            # kafna = str(np.int(sum(snow.nonmelt[name]))) 
            if snow.dplcs == 0:
                kafna = str(np.int(sum(snow.nonmelt[name]))) 
            else: 
                kafna = str(np.round(sum(snow.nonmelt[name]),snow.dplcs))            
    
            tmpna = ('unavail = %s')%(kafna)                     
            lbl.append(tmpna)                     
    
        # axs[iters].legend(lbl,loc='upper left') 
        ax.legend(lbl, loc = (0.025, 0.8),fontsize = 9)
        ax.set_ylim((0,ylim)) 
        for tick in ax.get_xticklabels():
            tick.set_rotation(30)             
    
        fig.tight_layout()
        ax.set_xticks(xts)
        ax.set_xlim(snow.xlims)
    
    fig.subplots_adjust(top=0.92,wspace = 0.1)
    fig.suptitle('SWE, %s'%snow.dateTo.date().strftime("%Y-%-m-%-d"))
          
    if snow.elv_flag:
        print('saving figure to %sswe_elev%s.png'%(snow.figs_path,snow.name_append))   
        plt.savefig('%sswe_elev%s.png'%(snow.figs_path,snow.name_append)) 
         
