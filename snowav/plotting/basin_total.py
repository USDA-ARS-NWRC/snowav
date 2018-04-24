
import numpy as np
from matplotlib import pyplot as plt
import matplotlib.colors as mcolors
from mpl_toolkits.axes_grid1 import make_axes_locatable
import seaborn as sns
import copy
import cmocean
import matplotlib.patches as mpatches
from datetime import datetime
import pandas as pd
from matplotlib.dates import DateFormatter
  
    
    
def basin_total(snow):       
         
    sns.set_style('darkgrid')
    sns.set_context("notebook")
    
    plt.close(8)
    fig,(ax,ax1) = plt.subplots(num=8, figsize=snow.figsize, 
                                dpi=snow.dpi, nrows = 1, ncols = 2)
    
    snow.barcolors.insert(0,'black')
    
    if snow.basin == 'LAKES' or snow.basin == 'RCEW':
        plotorder = [snow.plotorder[0]]
    else:
        plotorder = snow.plotorder
    
    for iters,name in enumerate(plotorder):
        # name = snow.plotorder[0]
        # iters = 0
        snow.state_summary[name].plot(ax=ax, color = snow.barcolors[iters])             
        ax1.plot(snow.accum_summary[name], 
                 color = snow.barcolors[iters], label='_nolegend_')
    
    ax1.yaxis.set_label_position("right")
    ax1.set_xlim((datetime(snow.wy -1 , 10, 1),snow.dateTo))
    ax.set_xlim((datetime(snow.wy - 1, 10, 1),snow.dateTo))
    ax1.tick_params(axis='y')
    ax1.yaxis.tick_right()
    ax.legend(loc='upper left')
    
    # Put on the same yaxis
    swey = ax.get_ylim()
    swiy = ax1.get_ylim()
    
    if swey[1] < swiy[1]:
        ax1.set_ylim((-0.1,swiy[1]))
        ax.set_ylim((-0.1,swiy[1]))
                    
    if swey[1] >= swiy[1]:
        ax1.set_ylim((-0.1,swey[1]))
        ax.set_ylim((-0.1,swey[1]))
    
    for tick,tick1 in zip(ax.get_xticklabels(),ax1.get_xticklabels()):
        tick.set_rotation(30) 
        tick1.set_rotation(30) 
         
    ax1.set_ylabel(r'[%s]'%(snow.vollbl))
    ax.axes.set_title('Basin SWE')
    ax1.axes.set_title('Accumulated Basin SWI')
    ax.set_ylabel(r'[%s]'%(snow.vollbl)) 
    
    # plt.tight_layout()   
    del snow.barcolors[0] 
    
    print('saving figure to %sbasin_total%s.png'%(snow.figs_path,
                                                  snow.name_append))   
    plt.savefig('%sbasin_total%s.png'%(snow.figs_path,snow.name_append))       
    
    ########################################
    #         Second figure
    ########################################
    
    # This needs to be improved...
    
    if snow.basin == 'BRB':
        accum_summary = snow.accum_summary
        main = 'Boise River Basin'
        multiswe = pd.DataFrame.from_csv(snow.summary_swe) 
        multiswi = pd.DataFrame.from_csv(snow.summary_swi)
        
        multiswe.wy17.iloc[304:] = 0
        multiswi.wy17 = np.cumsum(multiswi.wy17)
        multiswi.wy17.iloc[304:] = multiswi.wy17[303]
        
        state_summary = snow.state_summary.asfreq('D')
        
        # Put in this year
        multiswe.wy18.iloc[:len(state_summary[main])] = (
                                            state_summary[main].values) 
        multiswi.wy18.iloc[:len(accum_summary[main])] = (
                                            accum_summary[main].values)          
        
        if snow.units == 'SI':
            multiswe.wy17 = np.multiply(multiswe.wy17,0.00123348)
            multiswi.wy17 = np.multiply(multiswi.wy17,0.00123348)
            multiswe.wy16 = np.multiply(multiswe.wy16,0.00123348)
            multiswi.wy16 = np.multiply(multiswi.wy16,0.00123348)
            multiswe.wy15 = np.multiply(multiswe.wy15,0.00123348)
            multiswi.wy15 = np.multiply(multiswi.wy15,0.00123348)                
        
        plt.close(8)
        fig,(ax,ax1)    = plt.subplots(num=8, figsize=snow.figsize,
                                       dpi=snow.dpi, nrows = 1, ncols = 2)
    
        ax.plot(multiswe['wy15'], color = 'g',label = 'wy2015')
        ax.plot(multiswe['wy16'], color = 'r',label = 'wy2016')
        ax.plot(multiswe['wy17'], color = 'k',label = 'wy2017')
        ax.plot(multiswe['wy18'], color = 'b', label = 'wy2018')
    
        ax1.plot(multiswi['wy15'], color = 'g',label = 'wy2015')
        ax1.plot(multiswi['wy16'], color = 'r',label = 'wy2016')
        ax1.plot(multiswi['wy17'], color = 'k',label = 'wy2017')
        ax1.plot(multiswi['wy18'], color = 'b', label = 'wy2018')
        
        formatter = DateFormatter('%b')
        ax.xaxis.set_major_formatter(formatter)
        ax1.xaxis.set_major_formatter(formatter)
    
        ax1.yaxis.set_label_position("right")
        ax1.set_xlim((datetime(2017, 10, 1),datetime(2018, 8, 1)))
        ax.set_xlim((datetime(2017, 10, 1),datetime(2018, 8, 1)))
        ax1.tick_params(axis='y')
        ax1.yaxis.tick_right()
        ax.legend(loc='upper left')
        
        for tick,tick1 in zip(ax.get_xticklabels(),ax1.get_xticklabels()):
            tick.set_rotation(30) 
            tick1.set_rotation(30) 
             
        if snow.units == 'KAF':
            ax.set_ylabel('[KAF]')  
            ax1.set_ylabel('[KAF]')
            ax.axes.set_title('Water Year SWE')
            ax1.axes.set_title('Accumulated Basin SWI')
        
        if snow.units == 'SI':
            ax.set_ylabel(r'[$km^3$]') 
            ax1.set_ylabel(r'[$km^3$]')            
            ax.axes.set_title('Basin SWE [$km^3$]')
            ax1.axes.set_title('Accumulated Basin SWI [$km^3$]')
        
        ax.set_ylim((0,ax1.get_ylim()[1]))
        plt.tight_layout()      
        
        print('saving figure to %sbasin_total_multiyr%s.png'%(snow.figs_path,snow.name_append))   
        plt.savefig('%sbasin_total_multiyr%s.png'%(snow.figs_path,snow.name_append))   
        
          