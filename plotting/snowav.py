
import numpy as np
from matplotlib import pyplot as plt
# import matplotlib
from MidpointNormalize import MidpointNormalize
import matplotlib.colors as mcolors
from mpl_toolkits.axes_grid1 import make_axes_locatable
import seaborn as sns
from smrf import ipw
from shutil import copyfile
import os
import copy
import pandas as pd
from jinja2 import FileSystemLoader
from latex.jinja2 import make_env
from latex import build_pdf
import auto_report.str_replacer as sr
import ConfigParser as cfp
from datetime import datetime
import wyhr_to_datetime as wy
from operator import itemgetter
import utm
import mysql.connector
import netCDF4 as nc
import wyhr_to_datetime as wy
import cmocean

class snowav(object):

    def __init__(self,config_file):
        print('Reading the config file...')
        cfg         = cfp.ConfigParser()
        cfg.read(config_file)

        ####################################################
        #             Basin section                        #
        ####################################################
        self.basin          = cfg.get('Basin','basin')
        self.pixel          = int(cfg.get('Basin','pixel'))
        self.save_path      = cfg.get('Basin','save_path')
        self.name_append    = cfg.get('Basin','name_append')
        self.nrows          = int(cfg.get('Basin','nrows'))
        self.ncols          = int(cfg.get('Basin','ncols'))
        self.wy             = int(cfg.get('Basin','wy'))
        
        # Determine units and appropriate conversion
        self.units          = cfg.get('Basin','units')
        
        # KAF and inches
        if self.units == 'KAF':
            self.conversion_factor  = (self.pixel**2)*0.000000810713194*0.001
            self.depth_factor       = 0.03937
        
        # M m^3 and meters
        if self.units == 'SI':
            self.conversion_factor  = (self.pixel**2)*0.000000810713194*1233.48/1e6 
            self.depth_factor       = 0.001            
        
        # Would like to determine these automatically       
        # self.dateFrom    = cfg.get('Time','dateFrom')        
        # self.dateTo      = cfg.get('Time','dateTo')
        
        ####################################################
        #           Outputs section                        #
        ####################################################        
        
        # If the config file specifies csnowFile and psnowFile use those
        # and only read outputs between them, otherwise read in the full 
        # output directory
        self.snowband       = int(cfg.get('Outputs','snowband'))
        self.emband         = int(cfg.get('Outputs','emband'))          
        self.run_dir        = cfg.get('Outputs','run_dir')
        self.run_files      = sorted(os.listdir(self.run_dir))
        
        if cfg.has_option('Outputs','csnowFile') and cfg.has_option('Outputs','psnowFile'):
            self.psnowFile  = cfg.get('Outputs','psnowFile')
            self.csnowFile  = cfg.get('Outputs','csnowFile')
            self.cemFile    = cfg.get('Outputs','cemFile')
            
            spts            = self.csnowFile.split('.')
            sthr            = int(spts[-1])
            epts            = self.psnowFile.split('.')
            enhr            = int(epts[-1])
            run_files_filt  = []
            
            # This isn't very clean...
            for name in self.run_files:
                file_hr = 0
                if len(name) == 9 and name[0] is not 'e':
                    file_hr = int(name[5:9])
                    
                elif len(name) == 7:
                    file_hr = int(name[3:7])
                    
                if file_hr >= int(sthr) and file_hr <= int(enhr):
                    run_files_filt.append(name) 

            self.em_files   = [value for value in run_files_filt if 'em' in value]
            self.snow_files = copy.deepcopy(self.em_files)
            self.snow_files = [r.replace('em', 'snow') for r in self.snow_files]

        else:
            self.snow_files = [value for value in self.run_files if 'snow' in value]
            self.em_files   = [value for value in self.run_files if 'em' in value]
            self.psnowFile  = self.snow_files[0] 
            self.csnowFile  = self.snow_files[len(self.snow_files)-1] 
            self.cemFile    = self.em_files[len(self.em_files)-1] 
      
        ####################################################
        #           Accumulated                            #
        ####################################################   
        self.acc_clmin      = cfg.get('Accumulated','clmin')
        self.acc_clmax      = cfg.get('Accumulated','clmax')
        
        if cfg.has_option('Accumulated','ymin') and cfg.has_option('Accumulated','ymax'):
            self.acc_ylims  = (int(cfg.get('Accumulated','ymin')),int(cfg.get('Accumulated','ymax')))        

        ####################################################
        #           Elevation                              #
        ####################################################             
        if cfg.has_option('Elevation','ymin') and cfg.has_option('Elevation','ymax'):
            self.el_ylims   = (int(cfg.get('Elevation','ymin')),int(cfg.get('Elevation','ymax')))                      
        
        ####################################################
        #           Changes                                #
        ####################################################           
        self.ch_clmin       = cfg.get('Changes','clmin')
        self.ch_clmax       = cfg.get('Changes','clmax')
             
        if cfg.has_option('Changes','clminabs') and cfg.has_option('Changes','clmaxabs'):   
            self.ch_clminabs   = int(cfg.get('Changes','clminabs'))
            self.ch_clmaxabs   = int(cfg.get('Changes','clmaxabs'))
        if cfg.has_option('Changes','ymin') and cfg.has_option('Changes','ymax'):
            self.ch_ylims = (int(cfg.get('Changes','ymin')),int(cfg.get('Changes','ymax')))            
        
        ####################################################
        #           Results                                #
        ####################################################           
        self.r_clmin        = cfg.get('Results','clmin')
        self.r_clmax        = cfg.get('Results','clmax')     
        
        ####################################################
        #           DEM                                    #
        ####################################################         
        self.demPath        = cfg.get('DEM','demPath')
        self.total          = cfg.get('DEM','total')        
        self.subbasin1      = cfg.get('DEM','subbasin1')
        self.subbasin2      = cfg.get('DEM','subbasin2')
        self.subbasin3      = cfg.get('DEM','subbasin3')
        self.total_lbl      = cfg.get('DEM','total_lbl')
        self.sub1_lbl       = cfg.get('DEM','sub1_lbl')
        self.sub2_lbl       = cfg.get('DEM','sub2_lbl')
        self.sub3_lbl       = cfg.get('DEM','sub3_lbl')        

        ####################################################
        #          Plots                                   #
        ####################################################  
        self.figsize        = (int(cfg.get('Plots','fig_length')),int(cfg.get('Plots','fig_height')))
        self.dpi            = int(cfg.get('Plots','dpi'))
        if self.basin == 'BRB':
            self.barcolors = ['navy','royalblue','darkgrey','lightblue']
        if self.basin == 'TUOL':
            self.barcolors = ['darkgreen','palegreen','darkgrey','mediumseagreen']

        ####################################################
        #          Report                                  #
        ####################################################         
        self.env_path       = cfg.get('Report','env_path')
        self.tex_file       = cfg.get('Report','tex_file')
        self.rep_path       = cfg.get('Report','rep_path')
        self.conf_path      = cfg.get('Report','conf_path')
        self.templ_path     = cfg.get('Report','templ_path')
        self.rep_title      = cfg.get('Report','rep_title')

        
        ####################################################   
        # Now let's make some calculations, etc....
        
        # Set up the masks
        nrows               = self.nrows
        ncols               = self.ncols 

        if self.basin == 'BRB':
            self.masks  = { self.total_lbl: {'border': np.zeros((nrows,ncols)), 'mask': np.genfromtxt(self.total),'label':self.total_lbl,'runoff':[],'SWE':[]},
                            self.sub1_lbl: {'border': np.zeros((nrows,ncols)), 'mask': np.genfromtxt(self.subbasin1),'label':self.sub1_lbl,'runoff':[],'SWE':[]},
                            self.sub2_lbl: {'border': np.zeros((nrows,ncols)), 'mask': np.genfromtxt(self.subbasin2),'label':self.sub2_lbl,'runoff':[],'SWE':[]},
                            self.sub3_lbl: {'border': np.zeros((nrows,ncols)), 'mask': np.genfromtxt(self.subbasin3),'label':self.sub3_lbl,'runoff':[],'SWE':[]} 
                            } 
            
        elif self.basin == 'TUOL':
            self.masks  = { self.total_lbl: {'border': np.zeros((nrows,ncols)), 'mask': np.genfromtxt(self.total,skip_header=6),'label':self.total_lbl,'runoff':[],'SWE':[]},
                            self.sub1_lbl: {'border': np.zeros((nrows,ncols)), 'mask': np.genfromtxt(self.subbasin1,skip_header=6),'label':self.sub1_lbl,'runoff':[],'SWE':[]},
                            self.sub2_lbl: {'border': np.zeros((nrows,ncols)), 'mask': np.genfromtxt(self.subbasin2,skip_header=6),'label':self.sub2_lbl,'runoff':[],'SWE':[]},
                            self.sub3_lbl: {'border': np.zeros((nrows,ncols)), 'mask': np.genfromtxt(self.subbasin3,skip_header=6),'label':self.sub3_lbl,'runoff':[],'SWE':[]} 
                            }         
            self.barcolors = ['navy','royalblue','cornflowerblue','powderblue']

        # Get the DEM and set up the elevation bands
        self.dem            = np.genfromtxt(self.demPath,skip_header=6)
        if self.units == 'KAF':
            self.dem        = self.dem*3.28
            self.step       = 250
            self.edges      = np.arange(2750,13000+self.step,self.step)
            self.ixd        = np.digitize(self.dem,self.edges)        
        
        if self.units == 'SI':
            self.dem        = self.dem          
            self.step       = 100
            self.edges      = np.arange(800,3600+self.step,self.step)
            self.ixd        = np.digitize(self.dem,self.edges)             
        
        
        ####################################################   
        # These might be plot-specific...        
       
#         # Previous snow file, the starting point for change in SWE, etc
#         pdata           = ipw.IPW('%s%s'%(self.run_dir,self.psnowFile))
#         self.pswe       = pdata.bands[2].data
#         
#         # Current snow file, the ending point for change in SWE, etc
#         cdata           = ipw.IPW('%s%s'%(self.run_dir,self.csnowFile))
#         cemdata         = ipw.IPW('%s%s'%(self.run_dir,self.cemFile))
#         self.cswe       = cdata.bands[2].data
#         self.cc         = cemdata.bands[9].data         
# 
#         # Previous snow file, the starting point for change in SWE, etc
#         pdata           = ipw.IPW(self.psnowFile)
#         self.pswe       = pdata.bands[2].data
#     
#         # Current snow file, the ending point for change in SWE, etc
#         cdata           = ipw.IPW(self.csnowFile)
#         self.cemdata    = ipw.IPW(self.cemFile)
#         self.cswe       = cdata.bands[2].data
#         self.cc         = self.cemdata.bands[9].data   
        
                       
        # Finally, let's make a copy of the config file in the same place
        # that the figs will be saved. This copy of the config file will
        # get written with summary SWE/SWI info per sub basin
        pathf           = os.path.split(config_file)
        extf            = os.path.splitext(pathf[1])
        path_shr        = os.path.split(self.psnowFile)
        path_ehr        = os.path.split(self.csnowFile)
        ext_shr         = os.path.splitext(path_shr[1])
        ext_ehr         = os.path.splitext(path_ehr[1])
        
        self.figs_path  = self.save_path + '%s_%s/'%(ext_shr[1][1:5],ext_ehr[1][1:5])
        
        # Check run dates, and make a new folder if necessary
        if not os.path.exists(self.figs_path):
            os.makedirs(self.figs_path)

        self.config_copy = self.figs_path + extf[0] + self.name_append + '_%s_%s'%(ext_shr[1][1:5],ext_ehr[1][1:5]) + extf[1]
        
        # If it doesn't already exist, make it
        if not os.path.isfile(self.config_copy): 
            copyfile(config_file,self.config_copy)
            
        print('Finished reading config file')                   

    def calc(self):
        
        # Need to calc self.dateFrom and self.dateTo
        
        cclimit     = -5*1000*1000; #  based on an average of 60 W/m2 from my TL paper
        
        # Initialize containers for:
        #     accum: accumated values (swi,evap)
        #     state: (swe, cold content, depth, etc)
        accum           = np.zeros((self.nrows,self.ncols))
        state           = np.zeros((self.nrows,self.ncols))
        pstate          = np.zeros((self.nrows,self.ncols))
        cold            = np.zeros((self.nrows,self.ncols))
        accum_byelev    = pd.DataFrame(index = self.edges, columns = self.masks.keys()) 
        state_byelev    = pd.DataFrame(index = self.edges, columns = self.masks.keys())
        
        melt            = pd.DataFrame(index = self.edges, columns = self.masks.keys())
        nonmelt         = pd.DataFrame(index = self.edges, columns = self.masks.keys())
        
        # Loop through the iSnobal output files
        # Currently we need to load in each em.XXXX file to 'accumulate',
        # but only the first and last snow.XXXX file
        for iters,(em_name,snow_name) in enumerate(zip(self.em_files,self.snow_files)):
            
            em_file     = ipw.IPW('%s%s'%(self.run_dir,em_name))
            band        = em_file.bands[self.emband].data
            accum       = accum + band
            
            # When it is the first snow file, save that
            if snow_name == self.psnowFile.split('/')[-1]:
                snow_file   = ipw.IPW('%s%s'%(self.run_dir,snow_name))
                pstate      = snow_file.bands[self.snowband].data
            
            # When it hits the current snow file, load it in
            if snow_name == self.csnowFile.split('/')[-1]:
                snow_file   = ipw.IPW('%s%s'%(self.run_dir,snow_name))
                state       = snow_file.bands[self.snowband].data
                self.cold   = em_file.bands[9].data        

        # Time range based on ipw file headers
        self.dateFrom = wy.wyhr_to_datetime(self.wy,int(self.psnowFile.split('.')[-1]))
        self.dateTo = wy.wyhr_to_datetime(self.wy,int(self.csnowFile.split('.')[-1]))
        
        # Mask by subbasin and elevation band
        for iters,mask_name in enumerate(self.masks):
            # mask_name = 'Boise River Basin'
            accum_mask      = np.multiply(accum,self.masks[mask_name]['mask']) 
            elevbin         = np.multiply(self.ixd,self.masks[mask_name]['mask'])
            
            # Do it by elevation band
            for n in range(0,len(self.edges)):
                # n = 0
                
                ind         = elevbin == n
                state_bin   = state[ind]
                
                # Cold content
                ccb         = cold[ind]
                cind        = ccb > cclimit
                
                accum_byelev.loc[self.edges[n],mask_name]  = np.nansum(accum_mask[ind])
                melt.loc[self.edges[n],mask_name]        = np.nansum(state_bin[cind])
                nonmelt.loc[self.edges[n],mask_name]     = np.nansum(state_bin[~cind])  
        
            self.masks[mask_name]['SWE'] = (melt[mask_name].sum() + nonmelt[mask_name].sum())*self.conversion_factor    
        
        # Do some calculations
        delta_state     = state - pstate
        
        # Convert to desired units
        self.accum           = np.multiply(accum,self.depth_factor)
        self.state           = np.multiply(state,self.depth_factor)
        self.pstate          = np.multiply(pstate,self.conversion_factor)
        self.delta_state     = np.multiply(delta_state,self.conversion_factor)
        self.accum_byelev    = np.multiply(accum_byelev,self.conversion_factor)
        self.state_byelev    = np.multiply(state_byelev,self.conversion_factor)
        self.melt            = np.multiply(melt,self.conversion_factor)
        self.nonmelt         = np.multiply(melt,self.conversion_factor)
        
        # Write the config file
    
       
    # and then these should just be the plots themselves...
    def accumulated(self):

        qMin,qMax       = np.percentile(self.accum,[0,self.acc_clmax])
        clims           = (0,qMax)
        colors1         = cmocean.cm.deep(np.linspace(0., 1, 255))
        colors2         = plt.cm.Set2_r(np.linspace(0, 1, 1))
        colors          = np.vstack((colors2, colors1))
        mymap           = mcolors.LinearSegmentedColormap.from_list('my_colormap', colors)

        sns.set_style('darkgrid')
        sns.set_context("notebook")
        
        # This is to get the background white
        pmask           = self.masks[self.total_lbl]['mask']
        ixo             = pmask == 0
        self.accum[ixo] = np.nan
        ixz             = self.accum == 0
        self.accum[ixz] = 999
        mymap.set_bad('white',1.) 
        
        plt.close(0)
        fig,(ax,ax1)    = plt.subplots(num=0, figsize = self.figsize, dpi=self.dpi, nrows = 1, ncols = 2)      
        h               = ax.imshow(self.accum, clim = clims, interpolation='none', cmap = mymap)
        
        # Basin boundaries
        for name in self.masks:
            ax.contour(self.masks[name]['mask'],cmap = 'Greys',linewidths = 1)
                   
        # Do pretty stuff
        h.axes.get_xaxis().set_ticks([])
        h.axes.get_yaxis().set_ticks([])
       
        divider = make_axes_locatable(ax)
        cax     = divider.append_axes("right", size="5%", pad=0.2)
        cbar    = plt.colorbar(h, cax = cax)
        cbar.ax.tick_params() 
        
        if self.units == 'KAF':
            cbar.set_label('SWI [in]')
            h.axes.set_title('Snow Water Input Accumulated [in] \n %s to %s'%(self.dateFrom.date(),self.dateTo.date()))
        if self.units == 'SI':
            cbar.set_label('SWI [m]')
            h.axes.set_title('Snow Water Input Accumulated [m] \n %s to %s'%(self.dateFrom.date(),self.dateTo.date()))            
        
        # Sort by ascending melt amounts for the bar plots
        ordered_runoff     = np.zeros(len(self.masks))
        bname               = []
        for iters,name in enumerate(self.accum_byelev.columns.tolist()):
            # print(iters,name)
            ordered_runoff[iters] = self.accum_byelev[name].sum()
            bname.append(name)
        
        # Order by area
        # This is pretty ugly right now too, works in most cases, but could use work...
        x           = np.argsort(abs(ordered_runoff), axis=0)
        ordered_runoff_area = ordered_runoff[x]
        bname   = [bname[i] for i in x]
        
        # Sort ascending if the total value is negative, otherwise sort descending
        if self.accum_byelev[self.total_lbl].sum() < 0:
            sort_melt = sorted(enumerate(ordered_runoff_area), key=itemgetter(1))
        else:
            sort_melt = sorted(enumerate(ordered_runoff_area), key=itemgetter(1),reverse = True)
        
        # Plot the bar in order
        for iters,name in enumerate(sort_melt):
            a = sort_melt[iters][0] # where a is the index from lowest to highest...
            if self.units == 'KAF':
                b = ax1.bar(range(0,len(self.edges)),self.accum_byelev[bname[a]], label = '%s = %s KAF'%(bname[a],str(int(sort_melt[iters][1]))))
                for n in range(0,len(b)):
                    b[n].set_color(self.barcolors[iters])                
                
            if self.units == 'SI':
                b = ax1.bar(range(0,len(self.edges)),self.accum_byelev[bname[a]], label = r'%s = %s $M m^3$'%(bname[a],str(int(sort_melt[iters][1]))))
                for n in range(0,len(b)):
                    b[n].set_color(self.barcolors[iters])
        
        edges_lbl = self.edges[1::4]
        ax1.set_xticklabels(str(i) for i in edges_lbl)
        ax1.set_xlim(0,len(self.edges))
        if self.units == 'KAF':
            ax1.set_ylabel('KAF - per elevation band')
            ax1.set_xlabel('elevation [ft]')
        if self.units == 'SI':
            ax1.set_ylabel(r'M $m^3$ - per elevation band')
            ax1.set_xlabel('elevation [m]') 
                       
        ax1.yaxis.set_label_position("right")
        ax1.yaxis.tick_right()
        ylims = ax1.get_ylim()
        
        # If ylims were specified
        if hasattr(self,"asylims"):
            ax1.set_ylim(self.asylims)
        else:
            ax1.set_ylim((0,ylims[1]+ylims[1]*0.3))
        ax1.legend()
        
        plt.tight_layout()
        print('saving figure to %sresults%s'%(self.figs_path,self.name_append))
        # plt.savefig('%sswi%s.png'%(self.figs_path,self.save_app))   
        
    def current_image(self):
        
        # Color limits
        qMin,qMax   = np.nanpercentile(self.state,[self.r_clmin,self.r_clmax])
        clims       = (qMin,qMax)
        clims2      = (-10,0) # currently this is only appropriate for cold content
           
        # Prepare no-snow and outside of the basin for the colormaps   
        ixz                 = self.state == 0
        self.state[ixz]     = 999
        self.cold[ixz]      = -999 
        pmask               = self.masks[self.total_lbl]['mask']
        ixo                 = pmask == 0   
        
        # Colormap for self.state
        colorsbad       = plt.cm.Set2_r(np.linspace(0., 1, 1))
        colors1         = plt.cm.viridis_r(np.linspace(0., 1, 255))
        colors          = np.vstack((colors1,colorsbad))
        mymap           = mcolors.LinearSegmentedColormap.from_list('my_colormap', colors) 
        self.state[ixo] = np.nan        
        mymap.set_bad('white',1.)   
        
        # Colormap for cold content
        self.cold[ixo]  = np.nan
        colorsbad       = plt.cm.Set2_r(np.linspace(0., 1, 1))
        colors1         = plt.cm.plasma(np.linspace(0., 1, 255))
        colors          = np.vstack((colorsbad, colors1))
        mymap1          = mcolors.LinearSegmentedColormap.from_list('my_colormap', colors) 
        mymap1.set_bad('white',1.)                  
        
        sns.set_style('dark')
        sns.set_context("notebook")
        
        plt.close(3)
        fig,(ax,ax1)  = plt.subplots(num=3, figsize=self.figsize, dpi=self.dpi, nrows = 1, ncols = 2)
        h             = ax.imshow(self.state, clim = clims, interpolation='none', cmap = mymap)
        h1            = ax1.imshow(self.cold, clim = clims2, interpolation='none', cmap = mymap1)
        
        # Basin boundaries
        for name in self.masks:
            print(name)
            ax.contour(self.masks[name]['mask'],cmap = "Greys",linewidths = 1)
            ax1.contour(self.masks[name]['mask'],cmap = "Greys",linewidths = 1)
        
        # Do pretty stuff for the left plot
        h.axes.get_xaxis().set_ticks([])
        h.axes.get_yaxis().set_ticks([])
        divider = make_axes_locatable(ax)
        cax     = divider.append_axes("right", size="5%", pad=0.2)
        cbar    = plt.colorbar(h, cax = cax)
        cbar.ax.tick_params()  
        
        # Do pretty stuff for the right plot
        h1.axes.get_xaxis().set_ticks([])
        h1.axes.get_yaxis().set_ticks([])
        h1.axes.set_title('Cold Content [MJ/$m^3$] \n %s'%(self.dateTo.date()))
        divider = make_axes_locatable(ax1)
        cax     = divider.append_axes("right", size="5%", pad=0.2)
        cbar    = plt.colorbar(h1, cax = cax)
        cbar.set_label('Cold Content [MJ/$m^3$]')
        cbar.ax.tick_params() 
        
        # Labels that depend on units
        if self.units == 'KAF':
            h.axes.set_title('SWE [in] \n %s'%(self.dateTo.date()))
            cbar.set_label(r'SWE [in]')
        if self.units == 'SI':
            h.axes.set_title('SWE [m] \n %s'%(self.dateTo.date()))
            cbar.set_label(r'SWE [m]')
        
        plt.tight_layout() 
        print('saving figure to %sresults%s'%(self.figs_path,self.name_append))
        # plt.savefig('%sresults%s.png'%(self.figs_path,self.save_app))           
        
        
        
