
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
import mysql.connector

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
            
            # Going to change this to be the full path so that we can also go between two directories
            self.snow_files     = [self.run_dir + s for s in self.snow_files] 
            self.em_files       = [self.run_dir + s for s in self.em_files]  
            # self.psnowFile      = self.run_dir + self.psnowFile
            # self.csnowFile      = self.run_dir + self.csnowFile
            # self.cemFile        = self.run_dir + self.cemFile
            
            # If a previous run dir is listed, add those paths
            if cfg.has_option('Outputs','run_dir_p'):
                run_dir_p       = cfg.get('Outputs','run_dir_p')
                run_files_p     = [run_dir_p + s for s in sorted(os.listdir(run_dir_p))]
                snow_files_p    = [value for value in run_files_p if 'snow' in value]
                em_files_p      = [value for value in run_files_p if 'em' in value]
                self.snow_files = snow_files_p + self.snow_files
                self.em_files   = em_files_p + self.em_files
            
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
    
            ####################################################
            #          Plots                                   #
            ####################################################  
            self.figsize        = (int(cfg.get('Plots','fig_length')),int(cfg.get('Plots','fig_height')))
            self.dpi            = int(cfg.get('Plots','dpi'))
            if self.basin == 'BRB':
                # self.barcolors = ['navy','royalblue','darkgrey','lightblue']
                self.barcolors = ['darkgreen','palegreen','xkcd:dusty green','xkcd:hot green']  
            if self.basin == 'TUOL':
                # self.barcolors = ['darkgreen','palegreen','darkgrey','mediumseagreen']
                self.barcolors = ['darkgreen','palegreen','xkcd:dusty green','xkcd:hot green']  
            if self.basin == 'SJ':
                # self.barcolors = ['red','orangered','darkgrey','lightcoral'] 
                self.barcolors = ['darkgreen','palegreen','xkcd:dusty green','xkcd:vibrant green']     
    
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
            if self.basin == 'SJ':
                self.report_name    = ('SanJoaquinRiverBasin_SnowpackSummary_.pdf') 
                self.rep_title      = 'San Joaquin River Basin Snowpack Summary'
                
            # These are strings for the report
            if self.units == 'KAF':
                self.reportunits = 'KAF'
                
            if self.units == 'SI':
                self.reportunits = 'M m^3'
                           
            ####################################################   
            # Assign some basin-specific defaults and create mask dicts
            if self.basin == 'BRB':
                self.pixel          = 100
                self.nrows          = 1500
                self.ncols          = 1500 
                self.subbasin1      = cfg.get('DEM','subbasin1')
                self.subbasin2      = cfg.get('DEM','subbasin2')
                self.subbasin3      = cfg.get('DEM','subbasin3')
                self.total_lbl      = cfg.get('DEM','total_lbl')
                self.sub1_lbl       = cfg.get('DEM','sub1_lbl')
                self.sub2_lbl       = cfg.get('DEM','sub2_lbl')
                self.sub3_lbl       = cfg.get('DEM','sub3_lbl')  
                self.plotorder      = [self.total_lbl, self.sub2_lbl, self.sub1_lbl, self.sub3_lbl]                 
                
                self.masks  = { self.total_lbl: {'border': np.zeros((self.nrows,self.ncols)), 'mask': np.genfromtxt(self.total),'label':self.total_lbl},
                                self.sub1_lbl: {'border': np.zeros((self.nrows,self.ncols)), 'mask': np.genfromtxt(self.subbasin1),'label':self.sub1_lbl},
                                self.sub2_lbl: {'border': np.zeros((self.nrows,self.ncols)), 'mask': np.genfromtxt(self.subbasin2),'label':self.sub2_lbl},
                                self.sub3_lbl: {'border': np.zeros((self.nrows,self.ncols)), 'mask': np.genfromtxt(self.subbasin3),'label':self.sub3_lbl} 
                                } 
                
            elif self.basin == 'TUOL':
                self.pixel          = 50
                self.nrows          = 1339
                self.ncols          = 1374 
                self.subbasin1      = cfg.get('DEM','subbasin1')
                self.subbasin2      = cfg.get('DEM','subbasin2')
                self.subbasin3      = cfg.get('DEM','subbasin3')
                self.total_lbl      = cfg.get('DEM','total_lbl')
                self.sub1_lbl       = cfg.get('DEM','sub1_lbl')
                self.sub2_lbl       = cfg.get('DEM','sub2_lbl')
                self.sub3_lbl       = cfg.get('DEM','sub3_lbl') 
                self.plotorder      = [self.total_lbl, self.sub1_lbl, self.sub2_lbl, self.sub3_lbl]                               
                
                self.masks  = { self.total_lbl: {'border': np.zeros((self.nrows,self.ncols)), 'mask': np.genfromtxt(self.total,skip_header=6),'label':self.total_lbl},
                                self.sub1_lbl: {'border': np.zeros((self.nrows,self.ncols)), 'mask': np.genfromtxt(self.subbasin1,skip_header=6),'label':self.sub1_lbl},
                                self.sub2_lbl: {'border': np.zeros((self.nrows,self.ncols)), 'mask': np.genfromtxt(self.subbasin2,skip_header=6),'label':self.sub2_lbl},
                                self.sub3_lbl: {'border': np.zeros((self.nrows,self.ncols)), 'mask': np.genfromtxt(self.subbasin3,skip_header=6),'label':self.sub3_lbl} 
                                }  

            elif self.basin == 'SJ':
                self.pixel          = 50
                self.nrows          = 1657 
                self.ncols          = 1875 
                self.subbasin1      = cfg.get('DEM','subbasin1')
                self.subbasin2      = cfg.get('DEM','subbasin2')
                self.subbasin3      = cfg.get('DEM','subbasin3')
                self.total_lbl      = cfg.get('DEM','total_lbl') 
                self.sub1_lbl       = cfg.get('DEM','sub1_lbl')
                self.sub2_lbl       = cfg.get('DEM','sub2_lbl')
                self.sub3_lbl       = cfg.get('DEM','sub3_lbl')
                self.plotorder      = [self.total_lbl, self.sub3_lbl, self.sub2_lbl, self.sub1_lbl]                             
                
                self.masks  = { self.total_lbl: {'border': np.zeros((self.nrows,self.ncols)), 'mask': np.genfromtxt(self.total,skip_header=6),'label':self.total_lbl},
                                self.sub1_lbl: {'border': np.zeros((self.nrows,self.ncols)), 'mask': np.genfromtxt(self.subbasin1,skip_header=6),'label':self.sub1_lbl},
                                self.sub2_lbl: {'border': np.zeros((self.nrows,self.ncols)), 'mask': np.genfromtxt(self.subbasin2,skip_header=6),'label':self.sub2_lbl},
                                self.sub3_lbl: {'border': np.zeros((self.nrows,self.ncols)), 'mask': np.genfromtxt(self.subbasin3,skip_header=6),'label':self.sub3_lbl}
                                }  
              
            # Get the DEM 
            self.dem                    = np.genfromtxt(self.demPath,skip_header=6)
            
            # Bins
            if self.basin == 'BRB':
                if self.units == 'KAF':
                    emin        = 2750      # [ft]
                    emax        = 10500     # [ft]  
                    self.step   = 250
                if self.units == 'SI':
                    emin        = 800       # [m]
                    emax        = 3200      # [m]  
                    self.step   = 100                                 
            if self.basin == 'TUOL':
                if self.units == 'KAF':
                    emin        = 3000      # [ft]
                    emax        = 12000     # [ft]  
                    self.step   = 250
                if self.units == 'SI':
                    emin        = 800       # [m]
                    emax        = 3600      # [m]  
                    self.step   = 100  
            if self.basin == 'SJ':
                if self.units == 'KAF':
                    emin        = 1000      # [ft]
                    emax        = 13000     # [ft]  
                    self.step   = 500
                if self.units == 'SI':
                    emin        = 300       # [m]
                    emax        = 4000      # [m]  
                    self.step   = 200 
            
            # Do unit-specific things
            if self.units == 'KAF':
                self.conversion_factor  = (self.pixel**2)*0.000000810713194*0.001       # storage in KAF
                self.depth_factor       = 0.03937                                       # depth in inches
                self.dem                = self.dem*3.28
                self.edges              = np.arange(emin,emax+self.step,self.step)
                self.ixd                = np.digitize(self.dem,self.edges)        
            
            if self.units == 'SI':
                self.conversion_factor  = (self.pixel**2)*0.000000810713194*1233.48/1e6 # storage in M m^3
                self.depth_factor       = 0.001                                         # depth in meters
                self.dem                = self.dem          
                self.edges              = np.arange(emin,emax+self.step,self.step)
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
                
            # A few more basin-specific things
            if self.basin == 'BRB':
                self.xlims          = (0,31) # this goes from 2750 to 10500ft
            if self.basin == 'TUOL':
                self.xlims          = (0,len(self.edges))                
            if self.basin == 'SJ':
                self.xlims          = (0,len(self.edges))
            
            print('done.')    
        
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
  
        print('Processing iSnobal outputs...')
        
        cclimit             = -5*1000*1000  #  based on an average of 60 W/m2 from TL paper

        accum               = np.zeros((self.nrows,self.ncols))
        snowmelt            = np.zeros((self.nrows,self.ncols))
        state               = np.zeros((self.nrows,self.ncols))
        pstate              = np.zeros((self.nrows,self.ncols))
        state_byday         = np.zeros((self.nrows,self.ncols,len(self.snow_files)))
        accum_byelev        = pd.DataFrame(index = self.edges, columns = self.masks.keys()) 
        state_byelev        = pd.DataFrame(index = self.edges, columns = self.masks.keys())
        delta_state_byelev  = pd.DataFrame(index = self.edges, columns = self.masks.keys())
        melt                = pd.DataFrame(index = self.edges, columns = self.masks.keys())
        nonmelt             = pd.DataFrame(index = self.edges, columns = self.masks.keys())
        snowmelt_byelev     = pd.DataFrame(index = self.edges, columns = self.masks.keys())
        state_summary       = pd.DataFrame(columns = self.masks.keys())
        accum_summary       = pd.DataFrame(columns = self.masks.keys())
        
        
        # Loop through output files
        # Currently we need to load in each em.XXXX file to 'accumulate',
        # but only the first and last snow.XXXX file for changes
        for iters,(em_name,snow_name) in enumerate(zip(self.em_files,self.snow_files)):
            
            # Make date for pd index
            date        = wy.wyhr_to_datetime(self.wy,int(snow_name.split('.')[-1]))
            
            em_file     = ipw.IPW(em_name)
            band        = em_file.bands[self.emband].data
            accum       = accum + band
            snowmelt    = snowmelt + em_file.bands[7].data
            
            # load and calculate sub-basin total
            snow_file   = ipw.IPW(snow_name)
            tmpstate    = snow_file.bands[self.snowband].data  
            state_byday[:,:,iters] = tmpstate
            
            # Store daily sub-basin totals
            for mask_name in self.masks:
                accum_summary.loc[date, mask_name]   = np.nansum(np.multiply(accum,self.masks[mask_name]['mask'])) 
                state_summary.loc[date, mask_name]   = np.nansum(np.multiply(tmpstate,self.masks[mask_name]['mask'])) 
            
            # When it is the first snow file, copy
            if snow_name.split('/')[-1] == self.psnowFile.split('/')[-1]:
                pstate          = copy.deepcopy(tmpstate)
            
            # When it hits the current snow file, copy
            if snow_name.split('/')[-1] == self.csnowFile.split('/')[-1]:
                state           = copy.deepcopy(tmpstate)
                self.cold       = em_file.bands[9].data        

        # Time range based on ipw file headers
        self.dateFrom           = wy.wyhr_to_datetime(self.wy,int(self.psnowFile.split('.')[-1]))
        self.dateTo             = wy.wyhr_to_datetime(self.wy,int(self.csnowFile.split('.')[-1]))
        
        # Append date to report name
        parts                   = self.report_name.split('.')
        self.report_name        = parts[0] + self.dateTo.date().strftime("%Y%m%d") + '.' + parts[1]
        
        # Difference in state (SWE)
        delta_state         = state - pstate
        
        # Mask by subbasin and elevation band
        for mask_name in self.masks:
            accum_mask              = np.multiply(accum,self.masks[mask_name]['mask'])
            snowmelt_mask           = np.multiply(snowmelt,self.masks[mask_name]['mask'])
            state_mask              = np.multiply(state,self.masks[mask_name]['mask'])
            delta_state_byelev_mask = np.multiply(delta_state,self.masks[mask_name]['mask'])
            state_byelev_mask       = np.multiply(state,self.masks[mask_name]['mask'])
            elevbin                 = np.multiply(self.ixd,self.masks[mask_name]['mask'])
            
            # Do it by elevation band
            for n in np.arange(0,len(self.edges)):
                ind         = elevbin == n
                state_bin   = state_mask[ind]
                snowmelt_bin = snowmelt[ind]
                
                # Cold content
                ccb         = self.cold[ind]
                cind        = ccb > cclimit
                
                accum_byelev.loc[self.edges[n],mask_name]       = np.nansum(accum_mask[ind])
                snowmelt_byelev.loc[self.edges[n],mask_name]    = np.nansum(snowmelt_mask[ind])
                state_byelev.loc[self.edges[n],mask_name]       = np.nansum(state_byelev_mask[ind])
                delta_state_byelev.loc[self.edges[n],mask_name] = np.nansum(delta_state_byelev_mask[ind])
                melt.loc[self.edges[n],mask_name]               = np.nansum(state_bin[cind])
                nonmelt.loc[self.edges[n],mask_name]            = np.nansum(state_bin[~cind])    
        
            self.masks[mask_name]['SWE'] = (melt[mask_name].sum() + nonmelt[mask_name].sum())*self.conversion_factor        
        
        # Convert to desired units
        self.accum                  = np.multiply(accum,self.depth_factor)                  # sum over all time steps, spatial
        self.state                  = np.multiply(state,self.depth_factor)                  # at last time step, spatial
        self.state_byday            = np.multiply(state_byday,self.depth_factor)
        self.pstate                 = np.multiply(pstate,self.conversion_factor)            # at first time step, spatial
        self.delta_state            = np.multiply(delta_state,self.depth_factor)            # at last time step, spatial
        self.delta_state_byelev     = np.multiply(delta_state_byelev,self.conversion_factor)# at last time step, by elevation
        self.accum_byelev           = np.multiply(accum_byelev,self.conversion_factor)      # at last time step, by elevation
        self.snowmelt_byelev        = np.multiply(snowmelt_byelev,self.conversion_factor)      # at last time step, by elevation
        self.state_byelev           = np.multiply(state_byelev,self.conversion_factor)      # at last time step, by elevation
        self.melt                   = np.multiply(melt,self.conversion_factor)              # at last time step, based on cold content
        self.nonmelt                = np.multiply(nonmelt,self.conversion_factor)           # at last time step, based on cold content
        self.cold                   = np.multiply(self.cold,0.000001)                       # at last time step, spatial [MJ]
        self.state_summary          = np.multiply(state_summary,self.conversion_factor)     # daily by sub-basin
        self.accum_summary          = np.multiply(accum_summary,self.conversion_factor)     # daily by sub-basin
        
        print('done.')
        
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
        
        if self.basin == 'SJ':
            fix1 = np.arange(1275,1377)
            fix2 = np.arange(1555,1618)
            ax.plot(fix1*0,fix1,'k')
            ax.plot(fix2*0,fix2,'k')
                   
        # Do pretty stuff
        h.axes.get_xaxis().set_ticks([])
        h.axes.get_yaxis().set_ticks([])
        divider = make_axes_locatable(ax)
        cax     = divider.append_axes("right", size="5%", pad=0.2)
        cbar    = plt.colorbar(h, cax = cax)
        cbar.ax.tick_params() 
        
        if self.units == 'KAF':
            cbar.set_label('SWI [in]')
            h.axes.set_title('Accumulated SWI [in] \n %s to %s'%(self.dateFrom.date().strftime("%Y-%-m-%-d"),self.dateTo.date().strftime("%Y-%-m-%-d")))
        if self.units == 'SI':
            cbar.set_label('SWI [m]')
            h.axes.set_title('Accumulated SWI [m] \n %s to %s'%(self.dateFrom.date().strftime("%Y-%-m-%-d"),self.dateTo.date().strftime("%Y-%-m-%-d")))            
        
        # Bar plots
        for iters,name in enumerate(self.plotorder):
            snowmelt    = self.snowmelt_byelev[name]
            rain        = self.accum_byelev[name] - snowmelt
            if self.units == 'KAF':
                ax1.bar(range(0,len(self.edges)),snowmelt, color = self.barcolors[iters], edgecolor = 'k', hatch = '/////')
                plt.rcParams['hatch.linewidth'] = 0.5
                ax1.bar(range(0,len(self.edges)),rain, bottom = snowmelt, color = self.barcolors[iters], label = '%s = %s KAF'%(name,str(int(self.accum_byelev[name].sum()))))
                ax1.set_xlim(self.xlims)
                 
            if self.units == 'SI':
                ax1.bar(range(0,len(self.edges)),snowmelt, color = self.barcolors[iters], edgecolor = 'k', hatch = '/////')
                plt.rcParams['hatch.linewidth'] = 0.5
                ax1.bar(range(0,len(self.edges)),rain, bottom = snowmelt - rain, color = self.barcolors[iters], label = r'%s = %s $M m^3$'%(name,str(int(self.accum_byelev[name].sum()))))
                ax1.set_xlim(self.xlims)
        
        # Just for hatching legend entry
        ax1.bar(range(0,1),0.1, color = self.barcolors[iters], alpha = 0, hatch = '/////', label = 'snowmelt') 
        
        plt.tight_layout()
        xts         = ax1.get_xticks()
        edges_lbl   = []
        for i in xts[0:len(xts)-1]:
            edges_lbl.append(str(int(self.edges[int(i)])))
        
        ax1.set_xticklabels(str(i) for i in edges_lbl)
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
            
        ax1.legend(loc='upper left')

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
        mymap.set_under('darkslategrey',-1)     
        
        # Colormap for cold content
        colorsbad       = plt.cm.binary(np.linspace(0., 1, 1))
        colors1         = plt.cm.ocean_r(np.linspace(0., 1, 128))
        colors2         = plt.cm.YlOrRd(np.linspace(0., 1, 127))
        # colors2         = cmocean.cm.thermal(np.linspace(0., 1, 127))
        colors          = np.vstack((colorsbad, colors1, colors2))
        mymap1          = mcolors.LinearSegmentedColormap.from_list('my_colormap', colors) 
        cold[ixo]       = np.nan
        mymap1.set_bad('white')
        mymap1.set_over('darkslategrey',1) 
        
        # Prepare no-snow and outside of the basin for the colormaps   
        ixz             = state == 0
        state[ixz]      = -1
        cold[ixz]       = 1                  
        
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
            
        if self.basin == 'SJ':
            fix1 = np.arange(1275,1377)
            fix2 = np.arange(1555,1618)
            ax.plot(fix1*0,fix1,'k')
            ax.plot(fix2*0,fix2,'k')
            ax1.plot(fix1*0,fix1,'k')
            ax1.plot(fix2*0,fix2,'k')
        
        # Do pretty stuff for the left plot
        h.axes.get_xaxis().set_ticks([])
        h.axes.get_yaxis().set_ticks([])
        divider = make_axes_locatable(ax)
        cax     = divider.append_axes("right", size="5%", pad=0.2)
        cbar    = plt.colorbar(h, cax = cax, extend='both')
        oldlabels = cbar.ax.get_yticklabels()
        oldlabels[0] = 'snow\nfree'
        cbar.ax.set_yticklabels(oldlabels)
        
        cbar.ax.tick_params()  
        
        # Do pretty stuff for the right plot
        h1.axes.get_xaxis().set_ticks([])
        h1.axes.get_yaxis().set_ticks([])
        h1.axes.set_title('Cold Content [MJ/$m^3$] \n %s'%(self.dateTo.date().strftime("%Y-%-m-%-d")))
        divider = make_axes_locatable(ax1)
        cax     = divider.append_axes("right", size="5%", pad=0.2)
        cbar1    = plt.colorbar(h1, cax = cax, extend='max')
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
            colors      = np.vstack((colorsbad,colors1, colors2))
            mymap       = mcolors.LinearSegmentedColormap.from_list('my_colormap', colors)

        else:
            colors1     = cmocean.cm.matter_r(np.linspace(0., 1, 255))
            colors2     = plt.cm.Set2_r(np.linspace(0, 1, 1))
            colors      = np.vstack((colors1, colors2))
            mymap       = mcolors.LinearSegmentedColormap.from_list('my_colormap', colors)
    
        ixf                 = delta_state == 0
        delta_state[ixf]    = -1 # set snow-free  
        pmask               = self.masks[self.total_lbl]['mask']
        ixo                 = pmask == 0
        delta_state[ixo]    = np.nan
        cmap                = copy.copy(mymap)
        cmap.set_bad('white',1.)   
        cmap.set_over('darkslategrey')         
        
        sns.set_style('dark')
        sns.set_context("notebook")
       
        plt.close(2) 
        fig,(ax,ax1)    = plt.subplots(num=2, figsize=self.figsize, dpi=self.dpi, nrows = 1, ncols = 2)
        h               = ax.imshow(delta_state, clim = clims, interpolation='none', cmap = cmap, norm=MidpointNormalize(midpoint= 0 ,vmin = qMin,vmax = qMax))

        # Basin boundaries
        for name in self.masks:
            ax.contour(self.masks[name]['mask'],cmap = "Greys",linewidths = 1)
 
        if self.basin == 'SJ':
            fix1 = np.arange(1275,1377)
            fix2 = np.arange(1555,1618)
            ax.plot(fix1*0,fix1,'k')
            ax.plot(fix2*0,fix2,'k')    
        
        # Do pretty stuff
        h.axes.get_xaxis().set_ticks([])
        h.axes.get_yaxis().set_ticks([])
        divider = make_axes_locatable(ax)
        cax     = divider.append_axes("right", size="5%", pad=0.2)
        cbar    = plt.colorbar(h, cax = cax,extend='both')
        pos     = cbar.ax.get_position()
        
        oldlabels = cbar.ax.get_yticklabels()
        oldlabels[0] = 'snow\nfree'
        cbar.ax.set_yticklabels(oldlabels)
        
        cbar.ax.tick_params() 
        # cax.grid(False)
        
        # axb = cbar.ax.twinx()
        # axb.set_ylim(clims)
        
        if self.units == 'KAF':
            h.axes.set_title('Change in SWE [in] \n %s to %s'%(self.dateFrom.date().strftime("%Y-%-m-%-d"),self.dateTo.date().strftime("%Y-%-m-%-d")))
            cbar.set_label(r'$\Delta$ SWE [in]')           
        if self.units == 'SI':
            h.axes.set_title('Change in SWE [m] \n %s to %s'%(self.dateFrom.date().strftime("%Y-%-m-%-d"),self.dateTo.date().strftime("%Y-%-m-%-d")))
            cbar.set_label(r'$\Delta$ SWE [m]')   
      
        # Plot the bar in order
        for iters,name in enumerate(self.plotorder):
            if self.units == 'KAF':
                b = ax1.bar(range(0,len(self.edges)),self.delta_state_byelev[name], label = '%s = %s KAF'%(name,str(int(self.delta_state_byelev[name].sum()))))
                for n in range(0,len(b)):
                    b[n].set_color(self.barcolors[iters])
            if self.units == 'SI':
                b = ax1.bar(range(0,len(self.edges)),self.delta_state_byelev[name], label = r'%s = %s M $m^3$'%(name,str(int(self.delta_state_byelev[name].sum()))))    
                for n in range(0,len(b)):
                    b[n].set_color(self.barcolors[iters])
        
        ax1.set_xlim(self.xlims)
        plt.tight_layout()                            
        xts         = ax1.get_xticks()
        edges_lbl   = []
        for i in xts[0:len(xts)-1]:
            edges_lbl.append(str(int(self.edges[int(i)])))
  
        ax1.set_xticklabels(str(i) for i in edges_lbl)
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
            if ylims[0] == 0:
                ax1.set_ylim((ylims[0]+(ylims[0]*0.3),ylims[1]+ylims[1]*0.3))               
                         
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
        ax1.legend(loc='upper left')  
        ax1.grid(True)      
        plt.tight_layout() 
        fig.subplots_adjust(top=0.88)
        
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
               
        for iters,name in enumerate(self.plotorder):
            axs[iters].bar(range(0,len(self.edges)),self.melt[name], color = colors[0], bottom = self.nonmelt[name])
            axs[iters].bar(range(0,len(self.edges)),self.nonmelt[name], color = colors[1], label = 'nonmelt ')
            
            axs[iters].set_xlim(self.xlims)                          
            xts         = axs[iters].get_xticks()
            edges_lbl   = []
            for i in xts[0:len(xts)-1]:
                edges_lbl.append(str(int(self.edges[int(i)])))

            axs[iters].set_xticklabels(str(i) for i in edges_lbl)
            
            if iters > 1:
                if self.units == 'KAF':
                    axs[iters].set_xlabel('elevation [ft]')
                if self.units == 'SI':
                    axs[iters].set_xlabel('elevation [m]')                    
                  
            # Put yaxis on right 
            if iters == 1 or iters == 3:
                axs[iters].yaxis.set_label_position("right")
                axs[iters].yaxis.tick_right()
                
            if iters == 0 or iters == 1:
                axs[iters].set_xticklabels([])
            
            axs[iters].tick_params(axis='x')
            axs[iters].tick_params(axis='y')
            
            # Get basin total storage in strings for label
            kaf    = str(np.int(sum(self.melt[name]) + sum(self.nonmelt[name]))) 

            if self.units == 'KAF':          
                # axs[iters].axes.set_title('%s - %s KAF'%(self.masks[name]['label'],kaf))
                axs[iters].text(0.5,0.92,'%s - %s KAF'%(self.masks[name]['label'],kaf),horizontalalignment='center',transform=axs[iters].transAxes)
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

            # axs[iters].legend(lbl,loc='upper left') 
            axs[iters].legend(lbl,bbox_to_anchor=(0, 0, 0.6, 0.92),fontsize = 9)
            axs[iters].set_ylim((0,ylim)) 
            for tick in axs[iters].get_xticklabels():
                tick.set_rotation(30) 
        
        fig.tight_layout()
        for n in range(0,len(axs)):
            axs[n].set_xticks(xts)
            axs[n].set_xlim(self.xlims)
  
        fig.subplots_adjust(top=0.92)
        if self.units == 'KAF': 
            fig.suptitle('SWE [KAF], %s'%self.dateTo.date().strftime("%Y-%-m-%-d"))
        if self.units == 'SI': 
            fig.suptitle(r'SWE [M $m^3$], %s'%self.dateTo.date().strftime("%Y-%-m-%-d"))
              
        print('saving figure to %sswe_elev%s.png'%(self.figs_path,self.name_append))   
        plt.savefig('%sswe_elev%s.png'%(self.figs_path,self.name_append))  
        
    def basin_total(self): 
        
              
        sns.set_style('darkgrid')
        sns.set_context("notebook")

        plt.close(4)
        fig,(ax,ax1)    = plt.subplots(num=4, figsize=self.figsize, dpi=self.dpi, nrows = 1, ncols = 2)
        axb             = ax.twinx()
        axb.grid()

        for iters,name in enumerate(self.plotorder):
            self.state_summary[name].plot(ax=ax, color = self.barcolors[iters])             
            axb.plot(self.state_summary[name] - self.state_summary[name].iloc[0], color = self.barcolors[iters], linestyle=':')
            ax1.plot(self.accum_summary[name], color = self.barcolors[iters])

        ax1.yaxis.set_label_position("right")
        ax1.set_xlim((self.dateFrom,self.dateTo))
        ax1.tick_params(axis='y')
        ax1.yaxis.tick_right()
        ax.legend(loc='upper left')

        # Put on the same yaxis
        ax.set_ylim(ax1.get_ylim())
        axb.set_ylim(ax1.get_ylim())
        
        for tick,tick1 in zip(ax.get_xticklabels(),ax1.get_xticklabels()):
            tick.set_rotation(30) 
            tick1.set_rotation(30) 
             
        if self.units == 'KAF':
            ax.set_ylabel('storage [KAF]') 
            axb.set_ylabel('change during period [KAF]')  
            ax1.set_ylabel('SWI [KAF]')
            ax.axes.set_title('Total Basin SWE [KAF]')
            ax1.axes.set_title('Accumulated Basin SWI [KAF]')
        
        if self.units == 'SI':
            ax.set_ylabel(r'storage [M $m^3$]') 
            axb.set_ylabel(r'change during period [M $m^3$]')  
            ax1.set_ylabel('SWI [KAF]')            
            ax.axes.set_title('Total Basin SWE [M $m^3$]')
            ax1.axes.set_title('Accumulated Basin SWI [M $m^3$]')
        
        plt.tight_layout()      
        
        print('saving figure to %sbasin_total%s.png'%(self.figs_path,self.name_append))   
        plt.savefig('%sbasin_total%s.png'%(self.figs_path,self.name_append))        
             
    def stn_validate(self):
        stns        = ['ATAI1','BASI1','CCDI1','DHDI1','JKPI1','TRMI1']
        # stns        = ['CCDI1','COZI1','GGSI1','GLNI1','MRKI1','PRAI1']
        lbls        = ['Atlanta Summit','Banner Summit','Camas Creek','Dollarhide','Jackson Peak','Trinity Mountain']

        client      = 'BRB_2017'
        
        # get metadata from the data base from snotel sites
        qry         = ('SELECT tbl_metadata.* FROM tbl_metadata INNER JOIN tbl_stations ON tbl_metadata.primary_id=tbl_stations.station_id'
                       ' WHERE tbl_stations.client="'"%s"'" HAVING network_name = "'"SNOTEL"'";'%client)
        cnx         = mysql.connector.connect(user='markrobertson', password='whatdystm?1',host='10.200.28.137',database='weather_db')
        meta_sno    = pd.read_sql(qry, cnx)
        meta_sno    = meta_sno.loc[meta_sno['source'] == 'NRCS']
        meta_sno.index = meta_sno['secondary_id']
      
        swe_meas    = pd.DataFrame(index = pd.date_range(datetime(2017,10,1), self.dateTo, freq='D'),columns = stns)  
        swe_mod     = pd.DataFrame(index = pd.date_range(datetime(2017,10,1), self.dateTo, freq='D'),columns = stns)  
        tbl         = 'tbl_level1'
        var         = 'snow_water_equiv'
        st_time     = '2017-10-1 00:00:00'
        end_time    = self.dateTo.date().strftime("%Y-%-m-%-d")

        # Get Snotel station results
        for iters,stn in enumerate(stns): 
            cnx     = mysql.connector.connect(user='markrobertson', password='whatdystm?1',host='10.200.28.137', port='32768',database='weather_db')
            var_qry = ('SELECT weather_db.%s.date_time, weather_db.%s.%s ' % (tbl,tbl,var) +
                        'FROM weather_db.%s ' % tbl +
                        "WHERE weather_db.%s.date_time between '" % tbl + st_time+ "' and '"+end_time+"'"
                        "AND weather_db.%s.station_id IN ('" % tbl + stn + "');")
        
            data    = pd.read_sql(var_qry, cnx, index_col='date_time')
            dind    = pd.date_range(st_time,end_time,freq='D')
            swe_meas[stn]   = data.reindex(dind)
        
        # Now get pixel results
        ncpath  = self.run_dir.split('output')
        ncf     = nc.Dataset(ncpath[0] + 'snow.nc', 'r')    # open netcdf file
        nctvec  = ncf.variables['time'][:]
        vswe    = ncf.variables['specific_mass']            # get variable
        ncxvec  = ncf.variables['x'][:]                     # get x vec
        ncyvec  = ncf.variables['y'][:]                     # get y vec  
        
        for stn in stns:
            ll      = utm.from_latlon(meta_sno.ix[stn,'latitude'],meta_sno.ix[stn,'longitude']) # get utm coords from metadata
            xind    = np.where(abs(ncxvec-ll[0]) == min(abs(ncxvec-ll[0])))[0]  # get closest pixel index to the station
            yind    = np.where(abs(ncyvec-ll[1]) == min(abs(ncyvec-ll[1])))[0]  # get closest pixel index to the station
            swe     = pd.Series(vswe[:,yind,xind].flatten(),index=nctvec)  # pull out closest model pixel data
            swe_mod.loc[0:len(swe.values),stn] = swe.values
             
        # Plot
        maxswe = swe_meas.max()
        
        sns.set_style('darkgrid')
        sns.set_context('notebook')
        
        plt.close(6)
        fig, axs = plt.subplots(num = 6,figsize = (10,10),nrows = 3,ncols = 2)   
        axs = axs.flatten() 
        
        for iters,stn in enumerate(stns):
            axs[iters].plot(swe_meas[stn],'k')
            axs[iters].plot(swe_mod[stn],'b')
            
            axs[iters].set_title(lbls[iters])
            axs[iters].set_ylim((-0.1,max(maxswe) + max(maxswe)*0.05))
            axs[iters].set_xlim((self.dateFrom,self.dateTo))
            
            if iters == 1 or iters == 3 or iters == 5:
                axs[iters].yaxis.tick_right()
            
            if iters == 4 or iters == 5:
                print(iters)
                for tick in axs[iters].get_xticklabels():
                    tick.set_rotation(30) 
            else:
                axs[iters].set_xticklabels('') 

        axs[0].legend(['Snotel','model'],loc='upper left')
        axs[0].set_ylabel('SWE [mm]')
        
        plt.suptitle('Validation at Snotel Sites')
        plt.tight_layout()
        plt.subplots_adjust(top=0.92)
       
        print('saving figure to %svalidation%s.png'%(self.figs_path,self.name_append))   
        plt.savefig('%svalidation%s.png'%(self.figs_path,self.name_append))                  
