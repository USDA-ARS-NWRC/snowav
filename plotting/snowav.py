
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
import tqdm

class snowav(object):

    def __init__(self,config_file):
        '''
        Notes: 
        - pixel, nrows, ncols are assigned automatically for BRB, TUOL - other things could be still...
        - consider adding defaults and overwrite options for all plot labels (KAF per elevation band, 
            elevation [m], etc.), would eliminate a lot of if/else statements in plotting, and allow
            for plotting in cm too...
        
        '''
        
        try:
            print('Reading the config file...')
            cfg         = cfp.ConfigParser()
            cfg.read(config_file)
    
            ####################################################
            #             Basin section                        #
            ####################################################
            self.basin          = cfg.get('Basin','basin')
            self.save_path      = cfg.get('Basin','save_path')
            self.name_append    = cfg.get('Basin','name_append')
            self.wy             = int(cfg.get('Basin','wy'))
            
            # Determine units and appropriate conversion
            self.units          = cfg.get('Basin','units') 
            
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
                
                spts            = self.psnowFile.split('.')
                sthr            = int(spts[-1])
                epts            = self.csnowFile.split('.')
                enhr            = int(epts[-1])
                run_files_filt  = []

                # This could be made more robust
                try:
                    for name in self.run_files:
                        if not 'error.out' in name:
                            file_hr     = int(name.split('.')[-1])
  
                            if file_hr >= int(sthr) and file_hr <= int(enhr):
                                run_files_filt.append(name) 
                except:
                    print('error parsing self.run_files...')
                            
                self.em_files   = [value for value in run_files_filt if 'em' in value]
                self.snow_files = [value for value in run_files_filt if 'snow' in value]
    
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
            self.templ_path     = cfg.get('Report','templ_path')
            
            # These will later get appended with self.dateTo once it is calculated
            if self.basin == 'BRB':
                self.report_name    = ('BoiseRiverBasin_SnowpackSummary_.pdf')
                self.rep_title      = 'Boise River Basin Snowpack Summary'
            if self.basin == 'TUOL':
                self.report_name    = ('TuolumneRiverBasin_SnowpackSummary_.pdf') 
                self.rep_title      = 'Tuolumne River Basin Snowpack Summary'
                
            # These are strings for the report
            if self.units == 'KAF':
                self.reportunits = 'KAF'
                
            if self.units == 'SI':
                self.reportunits = 'M m^3'
                           
            ####################################################   
            # Assign some basin-specific defaults and create mask dicts
            if self.basin == 'BRB':
                self.pixel  = 100
                self.nrows  = 1500
                self.ncols  = 1500   
                
                self.masks  = { self.total_lbl: {'border': np.zeros((self.nrows,self.ncols)), 'mask': np.genfromtxt(self.total),'label':self.total_lbl},
                                self.sub1_lbl: {'border': np.zeros((self.nrows,self.ncols)), 'mask': np.genfromtxt(self.subbasin1),'label':self.sub1_lbl},
                                self.sub2_lbl: {'border': np.zeros((self.nrows,self.ncols)), 'mask': np.genfromtxt(self.subbasin2),'label':self.sub2_lbl},
                                self.sub3_lbl: {'border': np.zeros((self.nrows,self.ncols)), 'mask': np.genfromtxt(self.subbasin3),'label':self.sub3_lbl} 
                                } 
                
            elif self.basin == 'TUOL':
                self.pixel  = 50
                self.nrows  = 1339
                self.ncols  = 1374               
                
                self.masks  = { self.total_lbl: {'border': np.zeros((self.nrows,self.ncols)), 'mask': np.genfromtxt(self.total,skip_header=6),'label':self.total_lbl},
                                self.sub1_lbl: {'border': np.zeros((self.nrows,self.ncols)), 'mask': np.genfromtxt(self.subbasin1,skip_header=6),'label':self.sub1_lbl},
                                self.sub2_lbl: {'border': np.zeros((self.nrows,self.ncols)), 'mask': np.genfromtxt(self.subbasin2,skip_header=6),'label':self.sub2_lbl},
                                self.sub3_lbl: {'border': np.zeros((self.nrows,self.ncols)), 'mask': np.genfromtxt(self.subbasin3,skip_header=6),'label':self.sub3_lbl} 
                                }  
              
            # Get the DEM 
            self.dem                    = np.genfromtxt(self.demPath,skip_header=6)
            
            # Do unit-specific things
            if self.units == 'KAF':
                self.conversion_factor  = (self.pixel**2)*0.000000810713194*0.001       # storage in KAF
                self.depth_factor       = 0.03937                                       # depth in inches
                self.dem                = self.dem*3.28
                self.step               = 250
                self.edges              = np.arange(2750,13000+self.step,self.step)
                self.ixd                = np.digitize(self.dem,self.edges)        
            
            if self.units == 'SI':
                self.conversion_factor  = (self.pixel**2)*0.000000810713194*1233.48/1e6 # storage in M m^3
                self.depth_factor       = 0.001                                         # depth in meters
                self.dem                = self.dem          
                self.step               = 100
                self.edges              = np.arange(800,3600+self.step,self.step)
                self.ixd                = np.digitize(self.dem,self.edges)             
                                 
            # Make a copy of the config file in the same place that the figs will be saved
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
    
            # Only need to store this name if we decide to write more to the copied config file
            self.config_copy = self.figs_path + extf[0] + self.name_append + '_%s_%s'%(ext_shr[1][1:5],ext_ehr[1][1:5]) + extf[1]
            
            # If it doesn't already exist, make it
            if not os.path.isfile(self.config_copy): 
                copyfile(config_file,self.config_copy)
                
            print('finished reading config file')    
        
        except:
            print('error reading config file')

    def process(self):
        '''
        This function calculates everything we will need for the plots.
        
        Does not currently save to csv...
        
        accum:        is used for accumulated values in the 'em' files. This will 
                        typically be swi, but could also be evaporation if the proper band 
                        is specified in the config file
        state:        captures the current model state in the 'snow' files. This will 
                        typically be SWE, but could also be depth, etc.
        cold:         is cold content only
        delta_state:  is the difference between psnowFile and csnowFile
        _byelev:      are dataframes of totals for the period, by subbasin and elevation
        melt:         dataframe of totals for the period, defined by cold content
        nonmelt:      dataframe of totals for the period, defined by cold content
        '''
  
        cclimit         = -5*1000*1000  #  based on an average of 60 W/m2 from TL paper

        accum           = np.zeros((self.nrows,self.ncols))
        state           = np.zeros((self.nrows,self.ncols))
        pstate          = np.zeros((self.nrows,self.ncols))
        cold            = np.zeros((self.nrows,self.ncols))
        accum_byelev    = pd.DataFrame(index = self.edges, columns = self.masks.keys()) 
        state_byelev    = pd.DataFrame(index = self.edges, columns = self.masks.keys())
        delta_state_byelev    = pd.DataFrame(index = self.edges, columns = self.masks.keys())
        melt            = pd.DataFrame(index = self.edges, columns = self.masks.keys())
        nonmelt         = pd.DataFrame(index = self.edges, columns = self.masks.keys())
        
        # Loop through output files
        # Currently we need to load in each em.XXXX file to 'accumulate',
        # but only the first and last snow.XXXX file for changes
        for iters,(em_name,snow_name) in enumerate(zip(self.em_files,self.snow_files)):
            em_file     = ipw.IPW('%s%s'%(self.run_dir,em_name))
            band        = em_file.bands[self.emband].data
            accum       = accum + band
            
            # When it is the first snow file, save that
            if snow_name == self.psnowFile.split('/')[-1]:
                # snow_file   = ipw.IPW('%s%s'%(self.run_dir,snow_name))
                
                # If we force this back to read in self.psnowFile, rather than
                # '%s%s'%(self.run_dir,snow_name), then we can compare images from
                # two different run folders if it was specified that way in the config
                # file
                snow_file   = ipw.IPW(self.psnowFile)
                pstate      = snow_file.bands[self.snowband].data
            
            # When it hits the current snow file, load it in
            if snow_name == self.csnowFile.split('/')[-1]:
                snow_file   = ipw.IPW('%s%s'%(self.run_dir,snow_name))
                state       = snow_file.bands[self.snowband].data
                self.cold   = em_file.bands[9].data        

        # Time range based on ipw file headers
        self.dateFrom       = wy.wyhr_to_datetime(self.wy,int(self.psnowFile.split('.')[-1]))
        self.dateTo         = wy.wyhr_to_datetime(self.wy,int(self.csnowFile.split('.')[-1]))
        
        # Append date to report name
        parts               = self.report_name.split('.')
        self.report_name    = parts[0] + self.dateTo.date().strftime("%Y%m%d") + '.' + parts[1]
        
        # Difference in state (SWE)
        delta_state         = state - pstate
        
        # Mask by subbasin and elevation band
        for mask_name in self.masks:
            accum_mask              = np.multiply(accum,self.masks[mask_name]['mask'])
            state_mask              = np.multiply(state,self.masks[mask_name]['mask'])
            delta_state_byelev_mask = np.multiply(delta_state,self.masks[mask_name]['mask'])
            state_byelev_mask       = np.multiply(state,self.masks[mask_name]['mask'])
            elevbin                 = np.multiply(self.ixd,self.masks[mask_name]['mask'])
            
            # Do it by elevation band
            for n in range(0,len(self.edges)):
                ind         = elevbin == n
                state_bin   = state_mask[ind]
                
                # Cold content
                ccb         = cold[ind]
                cind        = ccb > cclimit
                
                accum_byelev.loc[self.edges[n],mask_name]       = np.nansum(accum_mask[ind])
                state_byelev.loc[self.edges[n],mask_name]       = np.nansum(state_byelev_mask[ind])
                delta_state_byelev.loc[self.edges[n],mask_name] = np.nansum(delta_state_byelev_mask[ind])
                melt.loc[self.edges[n],mask_name]               = np.nansum(state_bin[cind])
                nonmelt.loc[self.edges[n],mask_name]            = np.nansum(state_bin[~cind])  
        
            self.masks[mask_name]['SWE'] = (melt[mask_name].sum() + nonmelt[mask_name].sum())*self.conversion_factor        
        
        # Convert to desired units
        self.accum                  = np.multiply(accum,self.depth_factor)
        self.state                  = np.multiply(state,self.depth_factor)
        self.pstate                 = np.multiply(pstate,self.conversion_factor)
        self.delta_state            = np.multiply(delta_state,self.conversion_factor)
        self.delta_state_byelev     = np.multiply(delta_state_byelev,self.conversion_factor)
        self.accum_byelev           = np.multiply(accum_byelev,self.conversion_factor)
        self.state_byelev           = np.multiply(state_byelev,self.conversion_factor)
        self.melt                   = np.multiply(melt,self.conversion_factor)
        self.nonmelt                = np.multiply(nonmelt,self.conversion_factor)
        self.cold                   = np.multiply(self.cold,0.000001) # [MJ]
        
        # Consider writing summaries to the config file...?
    
       
    def accumulated(self):
        '''
        This function plots self.accum (typically SWI)
        
        '''
        # Make copy so that we can add nans for the plots, but not mess up the original
        accum           = copy.deepcopy(self.accum)

        qMin,qMax       = np.percentile(accum,[0,self.acc_clmax])
        clims           = (0,qMax)
        colors1         = cmocean.cm.dense(np.linspace(0., 1, 255))
        colors2         = plt.cm.binary(np.linspace(0, 1, 1))
        colors          = np.vstack((colors2, colors1))
        mymap           = mcolors.LinearSegmentedColormap.from_list('my_colormap', colors)

        sns.set_style('darkgrid')
        sns.set_context("notebook")
        
        # This is to get the background white
        pmask           = self.masks[self.total_lbl]['mask']
        ixo             = pmask == 0
        accum[ixo]      = np.nan
        # ixz             = accum == 0
        # accum[ixz]      = 999
        mymap.set_bad('white',1.) 
        
        plt.close(0)
        fig,(ax,ax1)    = plt.subplots(num=0, figsize = self.figsize, dpi=self.dpi, nrows = 1, ncols = 2)      
        h               = ax.imshow(accum, clim = clims, interpolation='none', cmap = mymap)
        
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
            h.axes.set_title('Snow Water Input Accumulated [in] \n %s to %s'%(self.dateFrom.date().strftime("%Y-%-m-%-d"),self.dateTo.date().strftime("%Y-%-m-%-d")))
        if self.units == 'SI':
            cbar.set_label('SWI [m]')
            h.axes.set_title('Snow Water Input Accumulated [m] \n %s to %s'%(self.dateFrom.date().strftime("%Y-%-m-%-d"),self.dateTo.date().strftime("%Y-%-m-%-d")))            
        
        # Sort by ascending melt amounts for the bar plots
        # This is ugly and could be improved...
        ordered_runoff          = np.zeros(len(self.masks))
        bname                   = []
        for iters,name in enumerate(self.accum_byelev.columns.tolist()):
            ordered_runoff[iters] = self.accum_byelev[name].sum()
            bname.append(name)
            
        x                       = np.argsort(abs(ordered_runoff), axis=0)
        ordered_runoff_area     = ordered_runoff[x]
        bname                   = [bname[i] for i in x]
        
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
        for tick in ax1.get_xticklabels():
            tick.set_rotation(30)        
        
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
        fig.subplots_adjust(top=0.88)
        print('saving figure to %sswi%s.png'%(self.figs_path,self.name_append))
        plt.savefig('%sswi%s.png'%(self.figs_path,self.name_append))   
        
    def current_image(self):
        '''
        This plots self.state (typically SWE) and self.cold
        
        '''
        # Make a copy so we can edit for plots
        state           = copy.deepcopy(self.state)
        cold            = copy.deepcopy(self.cold)
        
        # Color limits
        qMin,qMax       = np.nanpercentile(state,[self.r_clmin,self.r_clmax])
        clims           = (qMin,qMax)
        clims2          = (-10,0) # currently this is only appropriate for cold content
        
        # Areas outside basin
        pmask           = self.masks[self.total_lbl]['mask']
        ixo             = pmask == 0   
        
        # Colormap for self.state
        colorsbad       = plt.cm.binary(np.linspace(0., 1, 1))
        colors1         = cmocean.cm.haline_r(np.linspace(0., 1, 255))
        colors          = np.vstack((colors1,colorsbad))
        mymap           = mcolors.LinearSegmentedColormap.from_list('my_colormap', colors) 
        state[ixo]      = np.nan        
        mymap.set_bad('white',1.)   
        
        # Colormap for cold content
        colorsbad       = plt.cm.binary(np.linspace(0., 1, 1))
        colors1         = cmocean.cm.thermal(np.linspace(0., 1, 255))
        colors          = np.vstack((colorsbad, colors1))
        mymap1          = mcolors.LinearSegmentedColormap.from_list('my_colormap', colors) 
        cold[ixo]  = np.nan
        mymap1.set_bad('white',1.) 
        
        # Prepare no-snow and outside of the basin for the colormaps   
        ixz             = state == 0
        state[ixz]      = 999
        cold[ixz]       = -999                  
        
        sns.set_style('dark')
        sns.set_context("notebook")
        
        plt.close(3)
        fig,(ax,ax1)  = plt.subplots(num=3, figsize=self.figsize, dpi=self.dpi, nrows = 1, ncols = 2)
        h             = ax.imshow(state, clim = clims, interpolation='none', cmap = mymap)
        h1            = ax1.imshow(cold, clim=clims2, interpolation='none', cmap = mymap1)
        
        # Basin boundaries
        for name in self.masks:
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
        h1.axes.set_title('Cold Content [MJ/$m^3$] \n %s'%(self.dateTo.date().strftime("%Y-%-m-%-d")))
        divider = make_axes_locatable(ax1)
        cax     = divider.append_axes("right", size="5%", pad=0.2)
        cbar1    = plt.colorbar(h1, cax = cax)
        cbar1.set_label('Cold Content [MJ/$m^3$]')
        cbar1.ax.tick_params() 
        
        # Labels that depend on units
        if self.units == 'KAF':
            h.axes.set_title('SWE [in] \n %s'%(self.dateTo.date().strftime("%Y-%-m-%-d")))
            cbar.set_label(r'SWE [in]')
        if self.units == 'SI':
            h.axes.set_title('SWE [m] \n %s'%(self.dateTo.date().strftime("%Y-%-m-%-d")))
            cbar.set_label(r'SWE [m]') 
        
        plt.tight_layout() 
        print('saving figure to %sresults%s.png'%(self.figs_path,self.name_append))
        plt.savefig('%sresults%s.png'%(self.figs_path,self.name_append))  
        
    def image_change(self): 
        '''
        This plots self.delta_state
        
        Should add in more functionality for plotting differences between various images/runs
        
        '''        
        # Make copy so that we can add nans for the plots, but not mess up the original
        delta_state     = copy.deepcopy(self.delta_state)
        
        qMin,qMax       = np.percentile(delta_state,[self.ch_clmin,self.ch_clmax])
        clims           = (qMin,qMax) 
        
        # Override if absolute limits are provide in the config
        if hasattr(self,'ch_clminabs') and hasattr(self,'ch_clmaxabs'):
            clims       = (self.ch_clminabs,self.ch_clmaxabs)
     
        # if qMin is 0, need to change things 
        if qMax > 0:
            colorsbad   = plt.cm.Set2_r(np.linspace(0., 1, 1))           
            colors1     = cmocean.cm.matter_r(np.linspace(0., 1, 127))
            colors2     = plt.cm.Blues(np.linspace(0, 1, 128))
            colors      = np.vstack((colors1, colors2,colorsbad))
            mymap       = mcolors.LinearSegmentedColormap.from_list('my_colormap', colors)

        else:
            colors1     = cmocean.cm.matter_r(np.linspace(0., 1, 255))
            colors2     = plt.cm.Set2_r(np.linspace(0, 1, 1))
            colors      = np.vstack((colors1, colors2))
            mymap       = mcolors.LinearSegmentedColormap.from_list('my_colormap', colors)
    
        ixf                 = delta_state == 0
        delta_state[ixf]    = 999 # set snow-free  
        pmask               = self.masks[self.total_lbl]['mask']
        ixo                 = pmask == 0
        delta_state[ixo]    = np.nan
        cmap                = copy.copy(mymap)
        cmap.set_bad('white',1.)   
        cmap.set_over('lightslategrey') 
        
        # Sort by ascending melt amounts for the bar plots
        ordered_melt        = np.zeros(len(self.masks))
        bname               = []
        for iters,name in enumerate(self.delta_state_byelev):
            ordered_melt[iters] = self.delta_state_byelev[name].sum()
            bname.append(name)
        
        x               = np.argsort(abs(ordered_melt), axis=0)
        ordered_melt_area = ordered_melt[x]
        bname           = [bname[i] for i in x]
        
        # Sort ascending if the total value is negative, otherwise sort descending
        if sum(self.delta_state_byelev[self.total_lbl]) < 0:
            sort_melt = sorted(enumerate(ordered_melt_area), key=itemgetter(1))
        else:
            sort_melt = sorted(enumerate(ordered_melt_area), key=itemgetter(1),reverse = True)
        
        sns.set_style('darkgrid')
        sns.set_context("notebook")
       
        plt.close(2) 
        fig,(ax,ax1)    = plt.subplots(num=2, figsize=self.figsize, dpi=self.dpi, nrows = 1, ncols = 2)
        h               = ax.imshow(delta_state, clim = clims, interpolation='none', cmap = cmap, norm=MidpointNormalize(midpoint= 0 ,vmin = qMin,vmax = qMax))

        # Basin boundaries
        for name in self.masks:
            ax.contour(self.masks[name]['mask'],cmap = "Greys",linewidths = 1)
        
        # Do pretty stuff
        h.axes.get_xaxis().set_ticks([])
        h.axes.get_yaxis().set_ticks([])
        divider = make_axes_locatable(ax)
        cax     = divider.append_axes("right", size="5%", pad=0.2)
        cbar = plt.colorbar(h, cax = cax)
        cbar.ax.tick_params() 
        
        if self.units == 'KAF':
            h.axes.set_title('Change in SWE [in] \n %s to %s'%(self.dateFrom.date().strftime("%Y-%-m-%-d"),self.dateTo.date().strftime("%Y-%-m-%-d")))
            cbar.set_label(r'$\Delta$ SWE [in]')           
        if self.units == 'SI':
            h.axes.set_title('Change in SWE [m] \n %s to %s'%(self.dateFrom.date().strftime("%Y-%-m-%-d"),self.dateTo.date().strftime("%Y-%-m-%-d")))
            cbar.set_label(r'$\Delta$ SWE [m]')   
          
        # Plot the bar in order
        for iters,name in enumerate(sort_melt):
            a = sort_melt[iters][0] # where a is the index from lowest to highest...
            if self.units == 'KAF':
                b = ax1.bar(range(0,len(self.edges)),self.delta_state_byelev[bname[a]], label = '%s = %s KAF'%(bname[a],str(int(sort_melt[iters][1]))))
                for n in range(0,len(b)):
                    b[n].set_color(self.barcolors[iters])
            if self.units == 'SI':
                b = ax1.bar(range(0,len(self.edges)),self.delta_state_byelev[bname[a]], label = r'%s = %s M $m^3$'%(bname[a],str(int(sort_melt[iters][1]))))    
                for n in range(0,len(b)):
                    b[n].set_color(self.barcolors[iters])
                            
        edges_lbl = self.edges[1::4]
        ax1.set_xticklabels(str(i) for i in edges_lbl)
        ax1.set_xlim(0,len(self.edges))
        for tick in ax1.get_xticklabels():
            tick.set_rotation(30) 
             
        if hasattr(self,"ch_ylims"):
            ax1.set_ylim(self.ch_ylims)
        else:
            ylims = ax1.get_ylim()
            if ylims[0] < 0:
                ax1.set_ylim((ylims[0]+(ylims[0]*0.3),ylims[1]+ylims[1]*0.3))
            if ylims[1] == 0:
                ax1.set_ylim((ylims[0]+(ylims[0]*0.3),ylims[1]-ylims[0]*0.1))
           
        if self.units == 'KAF':
            ax1.set_ylabel('KAF - per elevation band')
            ax1.set_xlabel('elevation [ft]')
            ax1.axes.set_title('Change in SWE [KAF]')
        if self.units == 'SI':
            ax1.set_ylabel(r'M $m^3$')
            ax1.set_xlabel('elevation [m]')
            ax1.axes.set_title(r'Change in SWE [M $m^3$]')
               
        ax1.yaxis.set_label_position("right")
        ax1.tick_params(axis='x')
        ax1.tick_params(axis='y')
        ax1.yaxis.tick_right()
        ax1.legend(loc='best')   
        plt.tight_layout() 
        
        print('saving figure to %sswe_change%s.png'%(self.figs_path,self.name_append))
        plt.savefig('%sswe_change%s.png'%(self.figs_path,self.name_append))    
        
    def state_by_elev(self): 
        '''
        Plots SWE by elevation, delineated by melt/nonmelt
        
        '''
            
        lim     = np.max(self.melt[self.total_lbl]) + np.max(self.nonmelt[self.total_lbl])
        ylim    = np.max(lim) + np.max(lim)*0.2 
        colors  = ['xkcd:light mustard','xkcd:dark violet']
        fs      = list(self.figsize)
        fs[0]   = fs[0]*0.8
        fs      = tuple(fs)

        sns.set_style('darkgrid')
        sns.set_context("notebook")
       
        plt.close(1)
        fig,ax  = plt.subplots(num=1, figsize=fs, dpi=self.dpi, nrows = 2, ncols = 2)
        axs     = ax.ravel()
        
        # This is somewhat of a hack - we won't want alpha order, and my sort by size thing is weird sometimes...
        # Probably ought to pull this from somewhere instead...
        if self.basin == 'BRB':
            mask_order = ['Boise River Basin','Featherville','Twin Springs','Mores Creek']    
        if self.basin == 'TUOL':
            mask_order = ['Extended Tuolumne','Tuolumne','Cherry','Eleanor']
        
        for iters,name in enumerate(mask_order):
            sns.set_style('darkgrid')
            axs[iters].bar(range(0,len(self.edges)),self.melt[name], color = colors[0], bottom = self.nonmelt[name])
            axs[iters].bar(range(0,len(self.edges)),self.nonmelt[name], color = colors[1], label = 'nonmelt ')
            edges_lbl = self.edges[1::4]
            axs[iters].set_xticklabels(str(i) for i in edges_lbl)
            axs[iters].set_xlim((0,len(self.edges)))
            
            if iters > 1:
                if self.units == 'KAF':
                    axs[iters].set_xlabel('elevation [ft]')
                if self.units == 'SI':
                    axs[iters].set_xlabel('elevation [m]')                    
                
            axs[iters].tick_params('x',labelsize=8)    
            axs[iters].yaxis.set_label_position("right")
            axs[iters].tick_params(axis='x')
            axs[iters].tick_params(axis='y')
            axs[iters].yaxis.tick_right()
            
            # Get basin total storage in strings for label
            kaf    = str(np.int(sum(self.melt[name]) + sum(self.nonmelt[name]))) 

            if self.units == 'KAF':          
                axs[iters].axes.set_title('%s - %s KAF'%(self.masks[name]['label'],kaf))
                axs[iters].set_ylabel('KAF')
            if self.units == 'SI':          
                axs[iters].axes.set_title(r'%s - %s M $m^3$'%(self.masks[name]['label'],kaf))
                axs[iters].set_ylabel(r'M $m^3$')                
            
            lbl = []
            for n in (0,1):
                if n == 0:
                    kafa = str(np.int(sum(self.melt[name]))) 
                    if self.units == 'KAF':   
                        tmpa = ('melt = %s KAF')%(kafa)
                    if self.units == 'SI':   
                        tmpa = (r'melt = %s M $m^3$')%(kafa)  
                                          
                    lbl.append(tmpa)
                kafna = str(np.int(sum(self.nonmelt[name]))) 

                if self.units == 'KAF':
                    tmpna = ('nonmelt = %s KAF')%(kafna)
                if self.units == 'SI':
                    tmpna = (r'nonmelt = %s M $m^3$')%(kafna)                    
                lbl.append(tmpna)                     

            axs[iters].legend(lbl,loc='upper left') 
            axs[iters].set_ylim((0,ylim)) 
            for tick in axs[iters].get_xticklabels():
                tick.set_rotation(30) 
        
        fig.tight_layout()
        fig.subplots_adjust(top=0.88)
        if self.units == 'KAF': 
            fig.suptitle('SWE [KAF], %s'%self.dateTo.date().strftime("%Y-%-m-%-d"))
        if self.units == 'SI': 
            fig.suptitle(r'SWE [M $m^3$], %s'%self.dateTo.date().strftime("%Y-%-m-%-d"))
        
        
        print('saving figure to %sswe_elev%s.png'%(self.figs_path,self.name_append))   
        plt.savefig('%sswe_elev%s.png'%(self.figs_path,self.name_append))  
        
        
             
        
