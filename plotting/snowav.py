
import numpy as np
from matplotlib import pyplot as plt
from MidpointNormalize import MidpointNormalize
import matplotlib.colors as mcolors
from mpl_toolkits.axes_grid1 import make_axes_locatable
import seaborn as sns
from smrf import ipw
from shutil import copyfile
import os
import copy
import pandas as pd
import ConfigParser as cfp
from datetime import datetime
import utm
import netCDF4 as nc
import wyhr_to_datetime as wy
import cmocean
import mysql.connector
from matplotlib.dates import DateFormatter
import matplotlib.patches as mpatches

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
            self.dplcs          = int(cfg.get('Outputs','decimals'))         
            
            ##### do all of this at the end?
            self.run_files      = sorted(os.listdir(self.run_dir))
            
            if 'error.out' in self.run_files:
                self.run_files.remove('error.out')
            if '._error.out' in self.run_files:
                self.run_files.remove('._error.out')
            
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
                for name in self.run_files:
                    # if not 'error.out' or '._error.out' in name:
                    file_hr     = int(name.split('.')[-1])

                    if file_hr >= int(sthr) and file_hr <= int(enhr):
                        run_files_filt.append(name) 
         
                self.em_files   = [value for value in run_files_filt if 'em' in value]
                self.snow_files = [value for value in run_files_filt if not 'em' in value]
                
                # Going to change this to be the full path so that we can also go between two directories
                self.snow_files     = [self.run_dir + s for s in self.snow_files] 
                self.em_files       = [self.run_dir + s for s in self.em_files]       
            else:
                self.em_files       = [value for value in self.run_files if 'em' in value]
                self.snow_files     = [value for value in self.run_files if not 'em' in value]
                self.snow_files     = [self.run_dir + s for s in self.snow_files] 
                self.em_files       = [self.run_dir + s for s in self.em_files]  

                # Now define p and s
                self.psnowFile  = self.snow_files[0] 
                self.csnowFile  = self.snow_files[len(self.snow_files)-1] 
                self.cemFile    = self.em_files[len(self.em_files)-1]     
            
            if cfg.has_option('Outputs','csnowFile_flt') and cfg.has_option('Outputs','psnowFile_flt'):
                self.psnowFile_flt  = cfg.get('Outputs','psnowFile_flt')
                self.csnowFile_flt  = cfg.get('Outputs','csnowFile_flt')
                self.cemFile_flt    = cfg.get('Outputs','csnowFile_flt').replace('snow','em')
                
            ####################################################
            #           Runs                                   #
            ####################################################
        
            # Directories in 'Runs' are added underneath run_dir
            if cfg.has_section('Runs'): 
                runs = list(cfg.items('Runs'))   
                
                psfiles = []
                pefiles = []
                for rdir in runs:
                    run_files     = [rdir[1] + s for s in sorted(os.listdir(rdir[1]))]
                    if 'error.out' in run_files:
                        run_files.remove('error.out')
                    if '._error.out' in self.run_files:
                        self.run_files.remove('._error.out')                        
                    snow_files    = [value for value in run_files if not 'em' in value]
                    em_files      = [value for value in run_files if 'em' in value]
                    psfiles         = psfiles + snow_files
                    pefiles         = pefiles + em_files
                    
                self.snow_files = psfiles + self.snow_files
                self.em_files   = pefiles + self.em_files 
                
                if not cfg.has_option('Outputs','csnowFile') and cfg.has_option('Outputs','psnowFile'):
                
                    # Now define p and s
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
            self.barcolors      = ['xkcd:true green','palegreen','xkcd:dusty green','xkcd:vibrant green','red']     
    
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
            if self.basin == 'LAKES':
                self.report_name    = ('LakesBasin_SnowpackSummary_.pdf') 
                self.rep_title      = 'Lakes Basin Snowpack Summary'
              
            # These are strings for the report
            if self.units == 'KAF':
                self.reportunits = 'KAF'
                
            if self.units == 'SI':
                self.reportunits = 'km^3'
                           
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
                self.suborder       = [self.sub1_lbl,self.sub2_lbl,self.sub3_lbl]                 
                
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
                self.suborder       = [self.sub1_lbl,self.sub2_lbl,self.sub3_lbl]                              
                
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
                self.subbasin4      = cfg.get('DEM','subbasin4')
                self.total_lbl      = cfg.get('DEM','total_lbl') 
                self.sub1_lbl       = cfg.get('DEM','sub1_lbl')
                self.sub2_lbl       = cfg.get('DEM','sub2_lbl')
                self.sub3_lbl       = cfg.get('DEM','sub3_lbl')
                self.sub4_lbl       = cfg.get('DEM','sub4_lbl')
                self.plotorder      = [self.total_lbl, self.sub3_lbl, self.sub2_lbl, self.sub1_lbl,self.sub4_lbl] 
                self.suborder       = [self.sub1_lbl,self.sub2_lbl,self.sub3_lbl,self.sub4_lbl]                            
                
                self.masks  = { self.total_lbl: {'border': np.zeros((self.nrows,self.ncols)), 'mask': np.genfromtxt(self.total,skip_header=6),'label':self.total_lbl},
                                self.sub1_lbl: {'border': np.zeros((self.nrows,self.ncols)), 'mask': np.genfromtxt(self.subbasin1,skip_header=6),'label':self.sub1_lbl},
                                self.sub2_lbl: {'border': np.zeros((self.nrows,self.ncols)), 'mask': np.genfromtxt(self.subbasin2,skip_header=6),'label':self.sub2_lbl},
                                self.sub3_lbl: {'border': np.zeros((self.nrows,self.ncols)), 'mask': np.genfromtxt(self.subbasin3,skip_header=6),'label':self.sub3_lbl},
                                self.sub4_lbl: {'border': np.zeros((self.nrows,self.ncols)), 'mask': np.genfromtxt(self.subbasin4),'label':self.sub4_lbl}
                                }  

            elif self.basin == 'LAKES':
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
                
                self.masks  = { self.total_lbl: {'border': np.zeros((self.nrows,self.ncols)), 'mask': np.genfromtxt(self.total),'label':self.total_lbl},
                                self.sub1_lbl: {'border': np.zeros((self.nrows,self.ncols)), 'mask': np.genfromtxt(self.subbasin1),'label':self.sub1_lbl},
                                self.sub2_lbl: {'border': np.zeros((self.nrows,self.ncols)), 'mask': np.genfromtxt(self.subbasin2),'label':self.sub2_lbl},
                                self.sub3_lbl: {'border': np.zeros((self.nrows,self.ncols)), 'mask': np.genfromtxt(self.subbasin3),'label':self.sub3_lbl}
                                }  
            
            # Get the DEM 
            self.dem                    = np.genfromtxt(self.demPath,skip_header=6)
            
            # Bins
            if self.basin == 'BRB':
                if self.units == 'KAF':
                    emin        = 2500      # [ft]
                    emax        = 10500     # [ft]  
                    self.step   = 1000
                if self.units == 'SI':
                    emin        = 800       # [m]
                    emax        = 3200      # [m]  
                    self.step   = 250                                 
            if self.basin == 'TUOL':
                if self.units == 'KAF':
                    emin        = 3000      # [ft]
                    emax        = 12000     # [ft]  
                    self.step   = 1000
                if self.units == 'SI':
                    emin        = 800       # [m]
                    emax        = 3600      # [m]  
                    self.step   = 500  
            if self.basin == 'SJ': 
                if self.units == 'KAF':
                    emin        = 1000      # [ft]
                    emax        = 13000     # [ft]  
                    self.step   = 1000
                if self.units == 'SI':
                    emin        = 300       # [m]
                    emax        = 4000      # [m]  
                    self.step   = 500 
            if self.basin == 'LAKES':
                if self.units == 'KAF':
                    emin        = 8000      # [ft]
                    emax        = 12500     # [ft]  
                    self.step   = 500
                if self.units == 'SI':
                    emin        = 2400       # [m]
                    emax        = 3800      # [m]  
                    self.step   = 200                         
            
            # Do unit-specific things
            if self.units == 'KAF':
                self.conversion_factor  = (self.pixel**2)*0.000000810713194*0.001       # storage in KAF
                self.depth_factor       = 0.03937                                       # depth in inches
                self.dem                = self.dem*3.28
                self.edges              = np.arange(emin,emax+self.step,self.step)
                self.ixd                = np.digitize(self.dem,self.edges)     
                self.depthlbl           = 'in'
                self.vollbl             = 'KAF' 
                self.elevlbl            = 'ft'   
            
            if self.units == 'SI':
                self.conversion_factor  = (self.pixel**2)*0.000000810713194*1233.48/1e9 # storage in km^3
                # self.depth_factor       = 0.001                                       # depth in meters
                self.depth_factor       = 1                                             # depth in mm
                self.dem                = self.dem          
                self.edges              = np.arange(emin,emax+self.step,self.step)
                self.ixd                = np.digitize(self.dem,self.edges)     
                self.depthlbl           = 'mm'
                self.vollbl             = '$km^3$' 
                self.elevlbl            = 'm'     
                                 
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
                self.xlims          = (0,len(self.edges)) # this goes from 2750 to 10500ft
            if self.basin == 'TUOL':
                self.xlims          = (0,len(self.edges))                
            if self.basin == 'SJ':
                self.xlims          = (0,len(self.edges))
            if self.basin == 'LAKES':
                self.xlims          = (0,len(self.edges))
                self.imgx           = (1200,1375)
                self.imgy           = (425,225)
                       
            print('done.')    
        
        except:
            print('error reading config file')

    def process(self,*args):
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
        
        # If we add pre and current snow files for flights, force those
        # and add a few new things to calculate
        if len(args) != 0:
            self.psnowFile  = args[0]
            self.csnowFile  = args[1]
            self.cemFile    = args[1].replace('snow','em')
            self.snow_files = [self.psnowFile, self.csnowFile]
            self.em_files   = [self.psnowFile.replace('snow','em'), self.csnowFile.replace('snow','em')]
       
        cclimit             = -5*1000*1000  #  based on an average of 60 W/m2 from TL paper

        accum               = np.zeros((self.nrows,self.ncols))
        accum_sub           = np.zeros((self.nrows,self.ncols))
        snowmelt            = np.zeros((self.nrows,self.ncols))
        rain_bg             = np.zeros((self.nrows,self.ncols))
        precip              = np.zeros((self.nrows,self.ncols))
        precip_sub          = np.zeros((self.nrows,self.ncols))  
        snowmelt_sub        = np.zeros((self.nrows,self.ncols))
        state               = np.zeros((self.nrows,self.ncols))
        pstate              = np.zeros((self.nrows,self.ncols))
        state_byday         = np.zeros((self.nrows,self.ncols,len(self.snow_files)))
        accum_byelev        = pd.DataFrame(index = self.edges, columns = self.masks.keys())
        rain_bg_byelev      = pd.DataFrame(index = self.edges, columns = self.masks.keys())
        precip_byelev       = pd.DataFrame(index = self.edges, columns = self.masks.keys())
        precip_byelev_sub   = pd.DataFrame(index = self.edges, columns = self.masks.keys())
        accum_byelev_sub    = pd.DataFrame(index = self.edges, columns = self.masks.keys())  
        state_byelev        = pd.DataFrame(index = self.edges, columns = self.masks.keys())
        delta_state_byelev  = pd.DataFrame(index = self.edges, columns = self.masks.keys())
        melt                = pd.DataFrame(index = self.edges, columns = self.masks.keys())
        nonmelt             = pd.DataFrame(index = self.edges, columns = self.masks.keys())
        snowmelt_byelev     = pd.DataFrame(index = self.edges, columns = self.masks.keys())
        snowmelt_byelev_sub = pd.DataFrame(index = self.edges, columns = self.masks.keys())
        state_summary       = pd.DataFrame(columns = self.masks.keys())
        accum_summary       = pd.DataFrame(columns = self.masks.keys())
        precip_summary      = pd.DataFrame(columns = self.masks.keys())
               
        # Loop through output files
        # Currently we need to load in each em.XXXX file to 'accumulate',
        # but only the first and last snow.XXXX file for changes
        accum_sub_flag = False
        for iters,(em_name,snow_name) in enumerate(zip(self.em_files,self.snow_files)):
            # iters = 0
            # iters = iters + 1
            # em_name = self.em_files[iters]
            # snow_name = self.snow_files[iters]
            print(snow_name.split('.')[-1])
                     
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
            
            # Snow covered area and rain that fell on bare ground
            # sca         = tmpstate > 0
            
            # Get rain from input data
            sf          = em_name.replace('runs','data')
            sf          = sf.replace('run','data')
            sf          = sf.replace('output','ppt_4b')
            ppt_path    = sf.split('em')[0]
            
            out_hr      = snow_name.split('.')[-1]
            hrs         = range(int(out_hr) - 23,int(out_hr) + 1)
            
            pFlag       = False
            ppt_files   = []
            for hr in hrs:
                if os.path.isfile(ppt_path +'ppt.4b_'+ str(hr)):
                    pFlag = True
                    ppt_files = ppt_files + [ppt_path +'ppt.4b_'+ str(hr)]
            
            rain_hrly   = np.zeros((self.nrows,self.ncols))
            precip_hrly = np.zeros((self.nrows,self.ncols))
            
            # Load 'em in
            if pFlag:
                for pfile in ppt_files:
                    ppt             = ipw.IPW(pfile)
                    pre             = ppt.bands[0].data
                    percent_snow    = ppt.bands[1].data
                    rain_hrly       = rain_hrly + np.multiply(pre,(1-percent_snow))
                    precip_hrly     = precip_hrly + pre

            rain_bg     = rain_bg + rain_hrly 
            precip      = precip + precip_hrly    
                
            # Currently this is grabbing the second to last
            # Being used for flight difference
            # Definitely a better way to do this...
            if iters == (len(self.snow_files) - 1):
                # Mean pixel depth [mm] on psnowFile
                self.pre_pm     =  np.nansum(tmpstate*self.masks[self.total_lbl]['mask'])/self.masks[self.total_lbl]['mask'].sum()  
                self.pre_swe    =  np.nansum(tmpstate*self.masks[self.total_lbl]['mask'])*self.conversion_factor          
            
            # Store daily sub-basin totals
            for mask_name in self.masks:
                accum_summary.loc[date, mask_name]   = np.nansum(np.multiply(accum,self.masks[mask_name]['mask'])) 
                state_summary.loc[date, mask_name]   = np.nansum(np.multiply(tmpstate,self.masks[mask_name]['mask'])) 
                precip_summary.loc[date, mask_name]  = np.nansum(np.multiply(precip,self.masks[mask_name]['mask'])) 

            # When it is the first snow file, copy
            if snow_name.split('/')[-1] == self.psnowFile.split('/')[-1]:
                print('psnowfile is %s'%(snow_name))
                accum_sub_flag  = True
                pstate          = copy.deepcopy(tmpstate)
                
            # This only add between psnowFile and csnowFile
            if accum_sub_flag:
                accum_sub       = accum_sub + copy.deepcopy(band)
                snowmelt_sub    = snowmelt_sub + copy.deepcopy(em_file.bands[7].data) 
                precip_sub      = precip_sub + precip_hrly  
            
            # When it hits the current snow file, copy
            if snow_name.split('/')[-1] == self.csnowFile.split('/')[-1]:
                print('csnowfile is %s'%(snow_name))
                
                # Turn off, but last one will still have been added
                accum_sub_flag  = False
             
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
            rain_bg_mask            = np.multiply(rain_bg,self.masks[mask_name]['mask'])
            precip_mask             = np.multiply(precip,self.masks[mask_name]['mask'])
            precip_sub_mask         = np.multiply(precip_sub,self.masks[mask_name]['mask'])
            accum_mask_sub          = np.multiply(accum_sub,self.masks[mask_name]['mask'])
            snowmelt_mask           = np.multiply(snowmelt,self.masks[mask_name]['mask'])
            snowmelt_mask_sub       = np.multiply(snowmelt_sub,self.masks[mask_name]['mask'])
            state_mask              = np.multiply(state,self.masks[mask_name]['mask'])
            delta_state_byelev_mask = np.multiply(delta_state,self.masks[mask_name]['mask'])
            state_byelev_mask       = np.multiply(state,self.masks[mask_name]['mask'])
            elevbin                 = np.multiply(self.ixd,self.masks[mask_name]['mask'])
            
            # Do it by elevation band
            for n in np.arange(0,len(self.edges)):
                ind         = elevbin == n
                state_bin   = state_mask[ind]
                
                # Cold content
                ccb         = self.cold[ind]
                cind        = ccb > cclimit
                
                accum_byelev.loc[self.edges[n],mask_name]       = np.nansum(accum_mask[ind])
                rain_bg_byelev.loc[self.edges[n],mask_name]     = np.nansum(rain_bg_mask[ind])
                precip_byelev.loc[self.edges[n],mask_name]      = np.nansum(precip_mask[ind])
                precip_byelev_sub.loc[self.edges[n],mask_name]  = np.nansum(precip_sub_mask[ind])
                accum_byelev_sub.loc[self.edges[n],mask_name]   = np.nansum(accum_mask_sub[ind])
                snowmelt_byelev.loc[self.edges[n],mask_name]    = np.nansum(snowmelt_mask[ind])
                snowmelt_byelev_sub.loc[self.edges[n],mask_name] = np.nansum(snowmelt_mask_sub[ind])
                state_byelev.loc[self.edges[n],mask_name]       = np.nansum(state_byelev_mask[ind])
                delta_state_byelev.loc[self.edges[n],mask_name] = np.nansum(delta_state_byelev_mask[ind])
                melt.loc[self.edges[n],mask_name]               = np.nansum(state_bin[cind])
                nonmelt.loc[self.edges[n],mask_name]            = np.nansum(state_bin[~cind])    
        
            self.masks[mask_name]['SWE'] = (melt[mask_name].sum() + nonmelt[mask_name].sum())*self.conversion_factor        
        
        # Convert to desired units
        self.precip                 = np.multiply(precip,self.depth_factor) 
        self.accum                  = np.multiply(accum,self.depth_factor)                      # sum over all time steps, spatial
        self.accum_sub              = np.multiply(accum_sub,self.depth_factor)                  # sum over all time steps, spatial
        self.state                  = np.multiply(state,self.depth_factor)                      # at last time step, spatial
        self.state_byday            = np.multiply(state_byday,self.depth_factor)
        self.pstate                 = np.multiply(pstate,self.conversion_factor)                # at first time step, spatial
        self.delta_state            = np.multiply(delta_state,self.depth_factor)                # at last time step, spatial
        self.delta_state_byelev     = np.multiply(delta_state_byelev,self.conversion_factor)    # at last time step, by elevation
        self.accum_byelev           = np.multiply(accum_byelev,self.conversion_factor)          # at last time step, by elevation
        self.rain_bg_byelev         = np.multiply(rain_bg_byelev,self.conversion_factor)        # at last time step, by elevation
        self.precip_byelev          = np.multiply(precip_byelev,self.conversion_factor)        # at last time step, by elevation
        self.precip_byelev_sub      = np.multiply(precip_byelev_sub,self.conversion_factor)        # at last time step, by elevation
        self.accum_byelev_sub       = np.multiply(accum_byelev_sub,self.conversion_factor)      # at last time step, by elevation
        self.snowmelt_byelev        = np.multiply(snowmelt_byelev,self.conversion_factor)       # at last time step, by elevation
        self.snowmelt_byelev_sub    = np.multiply(snowmelt_byelev_sub,self.conversion_factor)   # at last time step, by elevation
        self.state_byelev           = np.multiply(state_byelev,self.conversion_factor)          # at last time step, by elevation
        self.melt                   = np.multiply(melt,self.conversion_factor)                  # at last time step, based on cold content
        self.nonmelt                = np.multiply(nonmelt,self.conversion_factor)               # at last time step, based on cold content
        self.cold                   = np.multiply(self.cold,0.000001)                           # at last time step, spatial [MJ]
        self.state_summary          = np.multiply(state_summary,self.conversion_factor)         # daily by sub-basin
        self.accum_summary          = np.multiply(accum_summary,self.conversion_factor)         # daily by sub-basin
        self.precip_summary         = np.multiply(precip_summary,self.conversion_factor)         # daily by sub-basin
        
        # Let's print some summary info...
        mask        = self.state_summary.index.to_series().diff() > pd.Timedelta('24:10:00')
        msum        = sum(mask)
        print('%s entries in dataframe index with gaps larger than 24h '%(msum))
        
        print('done.')
       
    def accumulated(self,*arg):
        '''
        This function plots self.accum (typically SWI)
        
        '''
        
        # if this is supplied, only report accum by the subsection between psnowFile and csnowFile 
        if len(arg) != 0:
            accum           = copy.deepcopy(self.accum_sub) 
            # snowmelt_byelev = copy.deepcopy(self.snowmelt_byelev_sub) 
            accum_byelev    = copy.deepcopy(self.accum_byelev_sub) 
                                       
        else:    
            accum           = copy.deepcopy(self.accum)
            # snowmelt_byelev = copy.deepcopy(self.snowmelt_byelev)
            # rain_bg_byelev  = copy.deepcopy(self.rain_bg_byelev)
            accum_byelev    = copy.deepcopy(self.accum_byelev)
            
        # qMin,qMax       = np.percentile(accum,[0,self.acc_clmax])
        qMin,qMax       = np.percentile(accum,[0,100])
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
        mymap.set_bad('white',1.) 
        
        # Now set SWI-free
        ixf             = accum == 0
        accum[ixf]      = -1
        mymap.set_under('grey',1.) 
        
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
            
        if self.basin == 'LAKES':
            ax.set_xlim(self.imgx)
            ax.set_ylim(self.imgy)
                   
        # Do pretty stuff
        h.axes.get_xaxis().set_ticks([])
        h.axes.get_yaxis().set_ticks([])
        divider = make_axes_locatable(ax)
        
        cax     = divider.append_axes("right", size="4%", pad=0.2)
        cbar    = plt.colorbar(h, cax = cax)
        cbar.ax.tick_params() 
        cbar.set_label('[%s]'%(self.depthlbl))
 
        h.axes.set_title('Accumulated SWI \n %s to %s'%(self.dateFrom.date().strftime("%Y-%-m-%-d"),self.dateTo.date().strftime("%Y-%-m-%-d")))  
        
        # Total basin label
        sumorder        = self.plotorder[1:]
        if self.basin == 'LAKES':
            sumorder = [self.plotorder[0]]
            
        if self.dplcs == 0:
            tlbl        = '%s = %s %s'%(self.plotorder[0],str(int(accum_byelev[self.plotorder[0]].sum())),self.vollbl)
        else: 
            tlbl        = '%s = %s %s'%(self.plotorder[0],str(np.round(accum_byelev[self.plotorder[0]].sum(),self.dplcs)),self.vollbl)
              
        for iters,name in enumerate(sumorder):
            print(name)
            if self.dplcs == 0:           
                lbl = '%s = %s %s'%(name,str(int(accum_byelev[name].sum())),self.vollbl)
            else:
                lbl = '%s = %s %s'%(name,str(np.round(accum_byelev[name].sum(),self.dplcs)),self.vollbl)
            # UMMMMM, does snowmelt include re-freezing...?
            # if any(snowmelt_byelev[name] > accum_byelev[name]):
            #     ix = snowmelt_byelev[name] > accum_byelev[name]
            #     snowmelt_byelev[name][ix] = accum_byelev[name][ix]
   
            if iters == 0:
                ax1.bar(range(0,len(self.edges)),accum_byelev[name], color = self.barcolors[iters], edgecolor = 'k',label = lbl)
                # ax1.bar(range(0,len(self.edges)),snowmelt_byelev[name], color = self.barcolors[iters], edgecolor = 'k',hatch = '/////')
            elif iters == 1:   
                ax1.bar(range(0,len(self.edges)),accum_byelev[name], bottom = accum_byelev[sumorder[iters-1]], color = self.barcolors[iters], edgecolor = 'k',label = lbl)
                # ax1.bar(range(0,len(self.edges)),snowmelt_byelev[name], bottom = accum_byelev[sumorder[iters-1]], color = self.barcolors[iters], edgecolor = 'k',hatch = '/////')
            elif iters == 2:   
                ax1.bar(range(0,len(self.edges)),accum_byelev[name], bottom = (accum_byelev[sumorder[iters-1]] + accum_byelev[sumorder[iters-2]]), color = self.barcolors[iters], edgecolor = 'k',label = lbl)
                # ax1.bar(range(0,len(self.edges)),snowmelt_byelev[name], bottom = (accum_byelev[sumorder[iters-1]] + accum_byelev[sumorder[iters-2]]), color = self.barcolors[iters], edgecolor = 'k',hatch = '/////')
            elif iters == 3:   
                ax1.bar(range(0,len(self.edges)),accum_byelev[name], 
                        bottom = (accum_byelev[sumorder[iters-1]] + accum_byelev[sumorder[iters-2]] + accum_byelev[sumorder[iters-3]]), 
                        color = self.barcolors[iters], edgecolor = 'k',label = lbl)
                
               
            plt.rcParams['hatch.linewidth'] = 1
            plt.rcParams['hatch.color'] = 'k'                
            ax1.set_xlim((self.xlims[0]-0.5,self.xlims[1]))

        # Just for hatching legend entry
        # ax1.bar(range(0,1),0.1, color = self.barcolors[iters], alpha = 0, hatch = '/////', label = 'snowmelt') 
        # ax1.bar(range(0,1),0.1, color = self.barcolors[iters], alpha = 0, hatch = '/////', label = 'rain') 
        # plt.rcParams['hatch.color'] = 'k'
        
        plt.tight_layout()
        xts         = ax1.get_xticks()
        edges_lbl   = []
        for i in xts[0:len(xts)-1]:
            edges_lbl.append(str(int(self.edges[int(i)])))
        
        ax1.set_xticklabels(str(i) for i in edges_lbl)
        for tick in ax1.get_xticklabels():
            tick.set_rotation(30)        

        ax1.set_ylabel('%s - per elevation band'%(self.vollbl))
        ax1.set_xlabel('elevation [%s]'%(self.elevlbl))
       
        ax1.yaxis.set_label_position("right")
        ax1.yaxis.tick_right()
        ylims = ax1.get_ylim()
        
        # If ylims were specified
        if hasattr(self,"asylims"):
            ax1.set_ylim(self.asylims)
        else:
            ax1.set_ylim((0,ylims[1]+ylims[1]*0.3))

        plt.tight_layout()
        fig.subplots_adjust(top=0.88)
      
        if self.basin != 'LAKES':
            # more ifs for number subs...
            if len(self.plotorder) == 5:
                ax1.legend(loc= (0.01,0.68))
            elif len(self.plotorder) == 4:
                ax1.legend(loc= (0.01,0.74))
            
        if self.basin == 'BRB':
            ax1.text(0.26,0.94,tlbl,horizontalalignment='center',transform=ax1.transAxes,fontsize = 10)
        else:
            ax1.text(0.3,0.94,tlbl,horizontalalignment='center',transform=ax1.transAxes,fontsize = 10)
        
        if sum(sum(ixf)) > 1000:
            patches = [mpatches.Patch(color='grey', label='no SWI')]
            if self.basin == 'SJ':
                ax.legend(handles=patches, bbox_to_anchor=(0.3, 0.05), loc=2, borderaxespad=0. )
            else:
                ax.legend(handles=patches, bbox_to_anchor=(0.05, 0.05), loc=2, borderaxespad=0. )
             
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
        # qMin,qMax       = np.nanpercentile(state,[self.r_clmin,self.r_clmax])
        # Finding that, when forcing colors for background and snow-free pixels,
        # if qMax is set to anything less than max(state), those pixels are white. 
        # Seems like this may be a bug in matplotlib?
        qMin,qMax       = np.nanpercentile(state,[0,100])
        clims           = (qMin,qMax)
        clims2          = (-10,0) # currently this is only appropriate for cold content
        
        ixw             = state >= qMax
        state[ixw]      = qMax - 0.01
        
        # Areas outside basin
        pmask           = self.masks[self.total_lbl]['mask']
        ixo             = pmask == 0   
        
        # Prepare no-snow and outside of the basin for the colormaps   
        ixz             = state == 0
        state[ixz]      = -1
        cold[ixz]       = 1     
        
        # Colormap for self.state
        colorsbad       = plt.cm.binary(np.linspace(0., 1, 1))
        colors1         = cmocean.cm.haline_r(np.linspace(0., 1, 255))
        colors          = np.vstack((colors1,colorsbad))
        mymap           = mcolors.LinearSegmentedColormap.from_list('my_colormap', colors) 
        
        state[ixo]      = np.nan        
        mymap.set_bad('white',1.) 
        mymap.set_under('lightgrey',1)    
        
        # Colormap for cold content
        # colorsbad       = plt.cm.binary(np.linspace(0., 1, 1))
        # colors1         = plt.cm.ocean_r(np.linspace(0., 1, 128))
        # colors2         = plt.cm.YlOrRd(np.linspace(0., 1, 127))
        # colors          = np.vstack((colorsbad, colors1, colors2))
        # mymap1          = mcolors.LinearSegmentedColormap.from_list('my_colormap', colors)
        mymap1          = plt.cm.RdYlBu_r 
        cold[ixo]       = np.nan
        mymap1.set_bad('white')
        mymap1.set_over('lightgrey',1) 
                         
        sns.set_style('dark')
        sns.set_context("notebook")
        
        plt.close(3)
        fig,(ax,ax1)  = plt.subplots(num=3, figsize=self.figsize, facecolor = 'white', dpi=self.dpi, nrows = 1, ncols = 2)
        h             = ax.imshow(state, clim = clims, cmap = mymap)
        h1            = ax1.imshow(cold, clim=clims2, cmap = mymap1)
        
        if self.basin == 'LAKES':
            ax.set_xlim(self.imgx)
            ax.set_ylim(self.imgy)
            ax1.set_xlim(self.imgx)
            ax1.set_ylim(self.imgy)
        
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
        cbar    = plt.colorbar(h, cax = cax)
        
        if self.units == 'KAF':         
            cbar.set_label('[in]')
        if self.units == 'SI':
            cbar.set_label('[mm]')
         
        cbar.ax.tick_params()                
        
        # Do pretty stuff for the right plot
        h1.axes.get_xaxis().set_ticks([])
        h1.axes.get_yaxis().set_ticks([])
        h1.axes.set_title('Cold Content \n %s'%(self.dateTo.date().strftime("%Y-%-m-%-d")))
        divider = make_axes_locatable(ax1)
        cax2     = divider.append_axes("right", size="5%", pad=0.2)
        cbar1    = plt.colorbar(h1, cax = cax2)
        cbar1.set_label('Cold Content [MJ/$m^3$]')
        cbar1.ax.tick_params() 
        
        h.axes.set_title('SWE \n %s'%(self.dateTo.date().strftime("%Y-%-m-%-d")))
        fig.subplots_adjust(top=0.95,bottom=0.05,right = 0.92, left = 0.05, wspace = 0.12)
        if self.basin == 'LAKES':
            plt.tight_layout()
        
        patches = [mpatches.Patch(color='grey', label='snow free')]
        if self.basin == 'SJ':
            ax.legend(handles=patches, bbox_to_anchor=(0.3, 0.05), loc=2, borderaxespad=0. )
        else:
            ax.legend(handles=patches, bbox_to_anchor=(0.05, 0.05), loc=2, borderaxespad=0. )
        
        print('saving figure to %sresults%s.png'%(self.figs_path,self.name_append))
        plt.savefig('%sresults%s.png'%(self.figs_path,self.name_append))  
        
    def image_change(self,*args): 
        '''
        This plots self.delta_state
        
        Should add in more functionality for plotting differences between various images/runs
        
        '''  
        
        # if len(args) != 0:
        #     self.name_append = self.name_append + args[0]
              
        # Make copy so that we can add nans for the plots, but not mess up the original
        delta_state     = copy.deepcopy(self.delta_state)
        
        qMin,qMax       = np.percentile(delta_state,[self.ch_clmin,self.ch_clmax])
        clims           = (qMin,qMax) 
        
        # Override if absolute limits are provide in the config
        if hasattr(self,'ch_clminabs') and hasattr(self,'ch_clmaxabs'):
            clims       = (self.ch_clminabs,self.ch_clmaxabs)
     
        # if qMin is 0, need to change things 
        if qMax > 0:
            colorsbad   = plt.cm.Accent_r(np.linspace(0., 1, 1))           
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
        delta_state[ixf]    = -100000 # set snow-free  
        pmask               = self.masks[self.total_lbl]['mask']
        ixo                 = pmask == 0
        delta_state[ixo]    = np.nan
        cmap                = copy.copy(mymap)
        cmap.set_bad('white',1.)   
        # cmap.set_over('darkslategrey',1) 
        # cmap.set_under('pink',1)        
        
        sns.set_style('dark')
        sns.set_context("notebook")
       
        plt.close(2) 
        fig,(ax,ax1)    = plt.subplots(num=2, figsize=self.figsize, dpi=self.dpi, nrows = 1, ncols = 2)
        h               = ax.imshow(delta_state, clim = clims, interpolation='none', cmap = cmap, norm=MidpointNormalize(midpoint= 0 ,vmin = qMin,vmax = qMax))

        if self.basin == 'LAKES':    
            ax.set_xlim(self.imgx)
            ax.set_ylim(self.imgy)

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
        cbar    = plt.colorbar(h, cax = cax)
        pos     = cbar.ax.get_position()
            
        cbar.ax.tick_params() 
        
        if self.units == 'KAF': 
            cbar.set_label(r'$\Delta$ SWE [in]')           
        if self.units == 'SI':
            cbar.set_label(r'$\Delta$ SWE [mm]')   
            
        h.axes.set_title('Change in SWE \n %s to %s'%(self.dateFrom.date().strftime("%Y-%-m-%-d"),self.dateTo.date().strftime("%Y-%-m-%-d")))    
      
        # Plot the bar in order
        sumorder    = self.plotorder[1:]  
        if self.basin == 'LAKES':
            sumorder = [self.plotorder[0]]  
        if self.dplcs == 0:
            tlbl        = '%s = %s %s'%(self.plotorder[0],str(int(self.delta_state_byelev[self.plotorder[0]].sum())),self.vollbl)
        else: 
            tlbl        = '%s = %s %s'%(self.plotorder[0],str(np.round(self.delta_state_byelev[self.plotorder[0]].sum(),self.dplcs)),self.vollbl)            
        
        for iters,name in enumerate(sumorder):  
            print(name) 
                
            if self.dplcs == 0:
                lbl  = '%s = %s %s'%(name,str(int(self.delta_state_byelev[name].sum())),self.vollbl)
            else: 
                lbl  = '%s = %s %s'%(name,str(np.round(self.delta_state_byelev[name].sum(),self.dplcs)),self.vollbl)                                 
 
            if iters == 0:
                ax1.bar(range(0,len(self.edges)),self.delta_state_byelev[name], color = self.barcolors[iters], edgecolor = 'k',label = lbl)
            elif iters == 1:   
                ax1.bar(range(0,len(self.edges)),self.delta_state_byelev[name], 
                        bottom = self.delta_state_byelev[sumorder[iters-1]], 
                        color = self.barcolors[iters], edgecolor = 'k',label = lbl)
            elif iters == 2:   
                ax1.bar(range(0,len(self.edges)),self.delta_state_byelev[name], 
                        bottom = self.delta_state_byelev[sumorder[iters-1]] + self.delta_state_byelev[sumorder[iters-2]], 
                        color = self.barcolors[iters], edgecolor = 'k',label = lbl)
            elif iters == 3:   
                ax1.bar(range(0,len(self.edges)),self.delta_state_byelev[name], 
                        bottom = self.delta_state_byelev[sumorder[iters-1]] + self.delta_state_byelev[sumorder[iters-2]] + self.delta_state_byelev[sumorder[iters-3]], 
                        color = self.barcolors[iters], edgecolor = 'k',label = lbl)
        
        
        ax1.set_xlim((self.xlims[0]-0.5,self.xlims[1]))
        # ax1.set_xlim(self.xlims)
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
            if ylims[0] < 0 and ylims[1] == 0:
                ax1.set_ylim((ylims[0]+(ylims[0]*0.3),ylims[1]+ylims[1]*0.3))
            if ylims[0] < 0 and ylims[1] > 0:
                ax1.set_ylim((ylims[0]+(ylims[0]*0.3),(ylims[1]-ylims[0])*0.5))                
            if ylims[1] == 0:
                ax1.set_ylim((ylims[0]+(ylims[0]*0.3),(-ylims[0])*0.5))
            if ylims[0] == 0:
                ax1.set_ylim((ylims[0]+(ylims[0]*0.3),ylims[1]+ylims[1]*0.3))               
                         
        if self.units == 'KAF':
            ax1.set_ylabel('KAF - per elevation band')
            ax1.set_xlabel('elevation [ft]')
            ax1.axes.set_title('Change in SWE')
        if self.units == 'SI':
            ax1.set_ylabel(r'$km^3$')
            ax1.set_xlabel('elevation [m]')
            ax1.axes.set_title(r'Change in SWE')
               
        ax1.yaxis.set_label_position("right")
        ax1.tick_params(axis='x')
        ax1.tick_params(axis='y')
        ax1.yaxis.tick_right()
        
        patches = [mpatches.Patch(color='grey', label='snow free')]
        if self.basin == 'SJ':
            ax.legend(handles=patches, bbox_to_anchor=(0.3, 0.05), loc=2, borderaxespad=0. )
        else:
            ax.legend(handles=patches, bbox_to_anchor=(0.05, 0.05), loc=2, borderaxespad=0. )
            
        if self.basin != 'LAKES':
            # more ifs for number subs...
            if len(self.plotorder) == 5:
                ax1.legend(loc= (0.01,0.68))
            elif len(self.plotorder) == 4:
                ax1.legend(loc= (0.01,0.74))
            
        if self.basin == 'BRB':
            ax1.text(0.27,0.94,tlbl,horizontalalignment='center',transform=ax1.transAxes,fontsize = 10)
        else:
            ax1.text(0.3,0.94,tlbl,horizontalalignment='center',transform=ax1.transAxes,fontsize = 10)
        
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
        ylim    = np.max(lim) + np.max(lim)*0.3 
        # colors  = ['xkcd:light mustard','xkcd:dark violet']
        colors  = ['xkcd:rose red','xkcd:cool blue']
        fs      = list(self.figsize)
        fs[0]   = fs[0]*0.8
        fs      = tuple(fs)

        sns.set_style('darkgrid')
        sns.set_context("notebook")
        
        nf = len(self.masks)
        
        if nf > 1 and nf < 5:
            nr = 2
            nc = 2
        elif nf < 7:
            nr = 2
            nc = 3
            
        if nf > 1:
            plt.close(1)
            fig,ax  = plt.subplots(num=1, figsize=fs, dpi=self.dpi, nrows = nr, ncols = nc)
            axs     = ax.ravel()
            if nf == 5:
                fig.delaxes(axs[5])
                   
            for iters,name in enumerate(self.plotorder):
                # iters = 0
                # name = self.plotorder[iters]
                axs[iters].bar(range(0,len(self.edges)),self.melt[name], color = colors[0], bottom = self.nonmelt[name])
                axs[iters].bar(range(0,len(self.edges)),self.nonmelt[name], color = colors[1], label = 'unavail ')
                
                axs[iters].set_xlim(self.xlims)                          
                xts         = axs[iters].get_xticks()
                edges_lbl   = []
                for i in xts[0:len(xts)-1]:
                    edges_lbl.append(str(int(self.edges[int(i)])))
    
                axs[iters].set_xticklabels(str(i) for i in edges_lbl)
                
                if iters > nc - 1:
                    if self.units == 'KAF':
                        axs[iters].set_xlabel('elevation [ft]')
                    if self.units == 'SI':
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
                kaf    = str(np.int(sum(self.melt[name]) + sum(self.nonmelt[name]))) 
                
                if self.dplcs == 0:
                    kaf    = str(np.int(sum(self.melt[name]) + sum(self.nonmelt[name]))) 
                else: 
                    kaf    = str(np.round(sum(self.melt[name]) + sum(self.nonmelt[name]),self.dplcs))             
    
                if self.units == 'KAF':          
                    # axs[iters].axes.set_title('%s - %s KAF'%(self.masks[name]['label'],kaf))
                    axs[iters].text(0.5,0.92,'%s - %s KAF'%(self.masks[name]['label'],kaf),horizontalalignment='center',transform=axs[iters].transAxes, fontsize = 10)
                if self.units == 'SI':          
                    axs[iters].text(0.5,0.92,'%s - %s $km^3$'%(self.masks[name]['label'],kaf),horizontalalignment='center',transform=axs[iters].transAxes)                
                
                if iters == 1 and nc == 3:
                    axs[iters].set_yticklabels([])
                else:
                    axs[iters].set_ylabel(self.units)
                
                lbl = []
                for n in (0,1):
                    if n == 0:
                        
                        if self.dplcs == 0:
                            kafa = str(np.int(sum(self.melt[name]))) 
                        else: 
                            kafa = str(np.round(sum(self.melt[name]),self.dplcs)) 
                        
                        # kafa = str(np.int(sum(self.melt[name]))) 
                        if self.units == 'KAF':   
                            tmpa = ('avail = %s KAF')%(kafa)
                        if self.units == 'SI':   
                            tmpa = (r'avail = %s $km^3$')%(kafa)  
                                              
                        lbl.append(tmpa)
                    
                    # kafna = str(np.int(sum(self.nonmelt[name]))) 
                    if self.dplcs == 0:
                        kafna = str(np.int(sum(self.nonmelt[name]))) 
                    else: 
                        kafna = str(np.round(sum(self.nonmelt[name]),self.dplcs))            
    
                    if self.units == 'KAF':
                        tmpna = ('unavail = %s KAF')%(kafna)
                    if self.units == 'SI':
                        tmpna = (r'unavail = %s $km^3$')%(kafna)                    
                    lbl.append(tmpna)                     
    
                # axs[iters].legend(lbl,bbox_to_anchor=(0, 0, 0.65, 0.92),fontsize = 9)
                axs[iters].legend(lbl, loc = (0.025, 0.65),fontsize = 9)
                axs[iters].set_ylim((0,ylim)) 
                for tick in axs[iters].get_xticklabels():
                    tick.set_rotation(30) 
            
            fig.tight_layout()
            for n in range(0,len(axs)):
                axs[n].set_xticks(xts)
                axs[n].set_xlim(self.xlims)       
        
        else:
            name = self.plotorder[0]
            plt.close(1)
            fig,ax  = plt.subplots(num=1, figsize=fs, dpi=self.dpi)
                   
            ax.bar(range(0,len(self.edges)),self.melt[name], color = colors[0], bottom = self.nonmelt[name])
            ax.bar(range(0,len(self.edges)),self.nonmelt[name], color = colors[1], label = 'unavail ')
            
            ax.set_xlim(self.xlims)                          
            xts         = ax.get_xticks()
            edges_lbl   = []
            for i in xts[0:len(xts)-1]:
                edges_lbl.append(str(int(self.edges[int(i)])))

            ax.set_xticklabels(str(i) for i in edges_lbl)

            if self.units == 'KAF':
                ax.set_xlabel('elevation [ft]')
            if self.units == 'SI':
                ax.set_xlabel('elevation [m]')                    
            
            ax.tick_params(axis='x')
            ax.tick_params(axis='y')
            
            # Get basin total storage in strings for label
            kaf    = str(np.int(sum(self.melt[name]) + sum(self.nonmelt[name]))) 
            
            if self.dplcs == 0:
                kaf    = str(np.int(sum(self.melt[name]) + sum(self.nonmelt[name]))) 
            else: 
                kaf    = str(np.round(sum(self.melt[name]) + sum(self.nonmelt[name]),self.dplcs))             

            if self.units == 'KAF':          
                # axs[iters].axes.set_title('%s - %s KAF'%(self.masks[name]['label'],kaf))
                ax.text(.25,0.95,'%s - %s KAF'%(self.masks[name]['label'],kaf),horizontalalignment='center',transform=ax.transAxes)
                ax.set_ylabel('KAF')
            if self.units == 'SI':          
                ax.text(-0.88,2.1,'%s - %s $km^3$'%(self.masks[name]['label'],kaf),horizontalalignment='center',transform=ax.transAxes)
                ax.set_ylabel(r'$km^3$')                
            
            lbl = []
            for n in (0,1):
                if n == 0:
                    if self.dplcs == 0:
                        kafa = str(np.int(sum(self.melt[name]))) 
                    else: 
                        kafa = str(np.round(sum(self.melt[name]),self.dplcs)) 
                    
                    # kafa = str(np.int(sum(self.melt[name]))) 
                    if self.units == 'KAF':   
                        tmpa = ('avail = %s KAF')%(kafa)
                    if self.units == 'SI':   
                        tmpa = (r'avail = %s $km^3$')%(kafa)  
                                          
                    lbl.append(tmpa)
                
                # kafna = str(np.int(sum(self.nonmelt[name]))) 
                if self.dplcs == 0:
                    kafna = str(np.int(sum(self.nonmelt[name]))) 
                else: 
                    kafna = str(np.round(sum(self.nonmelt[name]),self.dplcs))            

                if self.units == 'KAF':
                    tmpna = ('unavail = %s KAF')%(kafna)
                if self.units == 'SI':
                    tmpna = (r'unavail = %s $km^3$')%(kafna)                    
                lbl.append(tmpna)                     

            # axs[iters].legend(lbl,loc='upper left') 
            ax.legend(lbl, loc = (0.025, 0.8),fontsize = 9)
            ax.set_ylim((0,ylim)) 
            for tick in ax.get_xticklabels():
                tick.set_rotation(30)             
        
            fig.tight_layout()
            ax.set_xticks(xts)
            ax.set_xlim(self.xlims)
  
        fig.subplots_adjust(top=0.92,wspace = 0.1)
        if self.units == 'KAF': 
            fig.suptitle('SWE, %s'%self.dateTo.date().strftime("%Y-%-m-%-d"))
        if self.units == 'SI': 
            fig.suptitle(r'SWE, %s'%self.dateTo.date().strftime("%Y-%-m-%-d"))
              
        print('saving figure to %sswe_elev%s.png'%(self.figs_path,self.name_append))   
        plt.savefig('%sswe_elev%s.png'%(self.figs_path,self.name_append))  
        
    def basin_total(self):       
             
        sns.set_style('darkgrid')
        sns.set_context("notebook")

        plt.close(4)
        fig,(ax,ax1)    = plt.subplots(num=4, figsize=self.figsize, dpi=self.dpi, nrows = 1, ncols = 2)
        # axb             = ax.twinx()
        # axb.grid()

        self.barcolors.insert(0,'black')
        
        if self.basin == 'LAKES':
            plotorder = [self.plotorder[0]]
        else:
            plotorder = self.plotorder
        
        for iters,name in enumerate(plotorder):
            self.state_summary[name].plot(ax=ax, color = self.barcolors[iters])             
            # axb.plot(self.state_summary[name] - self.state_summary[name].iloc[0], color = self.barcolors[iters], linestyle=':')
            ax1.plot(self.accum_summary[name], color = self.barcolors[iters], label='_nolegend_')
            
#             if iters == 0:
#                 lbl = 'total precipitation'
#             else:
#                 lbl = '_nolegend_'
#             ax1.plot(self.precip_summary[name], linestyle = ':', linewidth = 1.0,color = self.barcolors[iters],label = lbl)
     
        ax1.yaxis.set_label_position("right")
        ax1.set_xlim((datetime(2017, 10, 1),self.dateTo))
        ax.set_xlim((datetime(2017, 10, 1),self.dateTo))
        ax1.tick_params(axis='y')
        ax1.yaxis.tick_right()
        ax.legend(loc='upper left')
        # ax1.legend(loc='upper left')
        
        # Put on the same yaxis
        swey = ax.get_ylim()
        swiy = ax1.get_ylim()
        
        if swey < swiy:
            ax1.set_ylim(swey)
            ax.set_ylim(swey)
        else:
            ax1.set_ylim(swiy)
            ax.set_ylim(swiy)
        # axb.set_ylim(ax1.get_ylim())
        
        for tick,tick1 in zip(ax.get_xticklabels(),ax1.get_xticklabels()):
            tick.set_rotation(30) 
            tick1.set_rotation(30) 
             
        ax1.set_ylabel(r'[%s]'%(self.vollbl))
        ax.axes.set_title('Basin SWE')
        ax1.axes.set_title('Accumulated Basin SWI')
        ax.set_ylabel(r'[%s]'%(self.vollbl)) 

        plt.tight_layout()      
        
        print('saving figure to %sbasin_total%s.png'%(self.figs_path,self.name_append))   
        plt.savefig('%sbasin_total%s.png'%(self.figs_path,self.name_append))       
        
        ########################################
        #         Second figure
        ########################################
        
        if self.basin == 'BRB':
            accum_summary       = self.accum_summary
            main                = 'Boise River Basin'
            multiswe            = pd.DataFrame.from_csv('/mnt/volumes/wkspace/results/brb/brb_multiyear_summary.csv') 
            multiswi            = pd.DataFrame.from_csv('/mnt/volumes/wkspace/results/brb/brb_multiyear_swi.csv')
            multiswi            = multiswi.cumsum() 
            
            multiswe.wy17.iloc[304:] = 0
            multiswi.wy17.iloc[304:] = multiswi.wy17[303]
            
            state_summary           = self.state_summary.asfreq('D')
            
            # Put in this year
            multiswe.wy18.iloc[:len(state_summary[main])] = state_summary[main].values 
            multiswi.wy18.iloc[:len(accum_summary[main])] = accum_summary[main].values          
            
            if self.units == 'SI':
                multiswe.wy17 = np.multiply(multiswe.wy17,0.00123348)
                multiswi.wy17 = np.multiply(multiswi.wy17,0.00123348)
                multiswe.wy16 = np.multiply(multiswe.wy16,0.00123348)
                multiswi.wy16 = np.multiply(multiswi.wy16,0.00123348)
            
            plt.close(5)
            fig,(ax,ax1)    = plt.subplots(num=5, figsize=self.figsize, dpi=self.dpi, nrows = 1, ncols = 2)
    
            ax.plot(multiswe['wy16'], color = 'r',label = 'wy2016')
            ax.plot(multiswe['wy17'], color = 'k',label = 'wy2017')
            ax.plot(multiswe['wy18'], color = 'b', label = 'wy2018')
       
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
    
            # Put on the same yaxis
            # ax.set_ylim(ax1.get_ylim())
            # axb.set_ylim(ax1.get_ylim())
            
            for tick,tick1 in zip(ax.get_xticklabels(),ax1.get_xticklabels()):
                tick.set_rotation(30) 
                tick1.set_rotation(30) 
                 
            if self.units == 'KAF':
                ax.set_ylabel('[KAF]')  
                ax1.set_ylabel('[KAF]')
                ax.axes.set_title('Water Year SWE')
                ax1.axes.set_title('Accumulated Basin SWI')
            
            if self.units == 'SI':
                ax.set_ylabel(r'[$km^3$]') 
                # axb.set_ylabel(r'change during period [M $m^3$]')  
                ax1.set_ylabel(r'[$km^3$]')            
                ax.axes.set_title('Basin SWE [$km^3$]')
                ax1.axes.set_title('Accumulated Basin SWI [$km^3$]')
            
            ax.set_ylim(ax1.get_ylim())
            plt.tight_layout()      
            
            print('saving figure to %sbasin_total_multiyr%s.png'%(self.figs_path,self.name_append))   
            plt.savefig('%sbasin_total_multiyr%s.png'%(self.figs_path,self.name_append))                   
             
    def stn_validate(self,rundirs):

        if self.basin == 'BRB':
            
            stns        = ['ATAI1','BASI1','CCDI1','DHDI1','JKPI1','TRMI1']
            # stns        = ['CCDI1','COZI1','GGSI1','GLNI1','MRKI1','PRAI1']
            lbls        = ['Atlanta Summit','Banner Summit','Camas Creek','Dollarhide','Jackson Peak','Trinity Mountain']
            client      = 'BRB_2017'

        ###############
        if self.basin == 'TUOL':
            stns        = ['DDM','LELC1','LVTC1','TIOC1','TUM','VGAC1']
#             rundirs = ['/mnt/data/blizzard/tuolumne/ops/wy2018/ops/runs/run20171001_20171217/',
#                 '/mnt/data/blizzard/tuolumne/ops/wy2018/ops/runs//run20171218_20180101/',
#                 '/mnt/data/blizzard/tuolumne/ops/wy2018/runs/run20180101_20180119/',
#                 '/mnt/data/blizzard/tuolumne/ops/wy2018/runs/run20180119_20180123/',
#                 '/mnt/data/blizzard/tuolumne/ops/wy2018/runs/run20180123_20180125/',
#                 '/mnt/data/blizzard/tuolumne/ops/wy2018/ops/runs/run20180125_20180206/' ]    
            lbls = copy.deepcopy(stns)
            client = 'TUOL_2017'
        
        ##############################
        # stns        = ['CHMC1','HNTC1','KSPC1','OSTC1','POSC1','RCKC1','TMRC1','DOEC1','DPPC1','MAMC1']
        # rundirs = ['/mnt/data/blizzard/sanjoaquin/ops/wy2018/runs/run20171001_20180125/',
        #            '/mnt/data/blizzard/sanjoaquin/ops/wy2018/ops/runs/run20180125_20180207/']    
        # lbls = copy.deepcopy(stns)
        # client = 'SJ_2017'
        #########################
        
        # get metadata from the data base from snotel sites
        if self.basin == 'BRB':
            qry         = ('SELECT tbl_metadata.* FROM tbl_metadata INNER JOIN tbl_stations ON tbl_metadata.primary_id=tbl_stations.station_id'
                       ' WHERE tbl_stations.client="'"%s"'" HAVING network_name = "'"SNOTEL"'";'%client)
        else:
            qry         = ('SELECT tbl_metadata.* FROM tbl_metadata INNER JOIN tbl_stations ON tbl_metadata.primary_id=tbl_stations.station_id WHERE tbl_stations.client="'"%s"'" ;'%client)        
        cnx         = mysql.connector.connect(user='markrobertson', password='whatdystm?1',host='10.200.28.137',database='weather_db')
        meta_sno    = pd.read_sql(qry, cnx)
        # meta_sno    = meta_sno.loc[meta_sno['source'] == 'NRCS']
        # meta_sno.index = meta_sno['secondary_id']
        meta_sno.index = meta_sno['primary_id']
      
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
            
        if self.basin == 'TUOL':
            swe_meas.TIOC1 = swe_meas.TIOC1 - 300    
        
        sns.set_style('darkgrid')
        sns.set_context('notebook')
        
        plt.close(6)
        fig, axs = plt.subplots(num = 6,figsize = (10,10),nrows = 3,ncols = 2)   
        axs = axs.flatten() 
        
        # First need to combine all nc files... 
        iswe = 0    
        for rname in rundirs:
            # rname = rundirs[0]

            print(rname)
            
            ncpath  = rname
            ncf     = nc.Dataset(ncpath + 'snow.nc', 'r')    # open netcdf file
            nctvec  = ncf.variables['time'][:]
            vswe    = ncf.variables['specific_mass']            # get variable
            ncxvec  = ncf.variables['x'][:]                     # get x vec
            ncyvec  = ncf.variables['y'][:]                     # get y vec      
            
            for iters,stn in enumerate(stns):
                # iters = 0
                # stn = stns[iters]
                ll      = utm.from_latlon(meta_sno.ix[stn,'latitude'],meta_sno.ix[stn,'longitude']) # get utm coords from metadata

                xind    = np.where(abs(ncxvec-ll[0]) == min(abs(ncxvec-ll[0])))[0]  # get closest pixel index to the station
                yind    = np.where(abs(ncyvec-ll[1]) == min(abs(ncyvec-ll[1])))[0]  # get closest pixel index to the station
                
                axs[iters].plot(swe_meas[stn],'k',label='measured')
                
                px = (1,1,1,0,0,0,-1,-1,-1)
                py = (1,0,-1,1,0,-1,1,0,-1)
                
                for n,m in zip(px,py): 
                    swe     = pd.Series(vswe[:,yind+m,xind+n].flatten(),index=nctvec)  # pull out closest model pixel data
                    swe_mod.loc[iswe:(iswe + len(swe.values)),stn] = swe.values
                    axs[iters].plot(swe_mod[stn],'b',linewidth = 0.75,label='model') 
        
                
                axs[iters].set_title(lbls[iters])
                axs[iters].set_xlim((datetime(2017, 10, 1),self.dateTo))
                
                if iters == 1 or iters == 3 or iters == 5:
                    axs[iters].yaxis.tick_right()
                
                if iters == 4 or iters == 5:
                    for tick in axs[iters].get_xticklabels():
                        tick.set_rotation(30) 
                else:
                    axs[iters].set_xticklabels('') 

            ncf.close()   
            iswe = iswe + len(swe.values)
             
        # Plot
        maxm = np.nanmax(swe_meas)
        maxi = np.nanmax(swe_mod)
        
        if maxm > maxi:
            maxswe = maxm
        else:
            maxswe = maxi
            
        for iters in range(0,len(stns)):
            axs[iters].set_ylim((-0.1,maxswe + maxswe*0.05))     

        axs[0].legend(['measured','modelled'],loc='upper left')
        axs[0].set_ylabel('SWE [mm]')
        
        plt.suptitle('Validation at Measured Sites')
        plt.tight_layout()
        plt.subplots_adjust(top=0.92)
       
        print('saving figure to %svalidation%s.png'%(self.figs_path,self.name_append))   
        plt.savefig('%svalidation%s.png'%(self.figs_path,self.name_append))                  
        
    def basin_detail(self):
        import math
        import copy
        import numpy as np
        
        demcopy = copy.copy(self.dem)
        pmask           = self.masks[self.total_lbl]['mask']
        ixo             = pmask == 0
        demcopy[ixo]      = np.nan

        emin        = int(math.ceil(np.nanmin(demcopy)/10.0))*10
        emax        = int(math.ceil(np.nanmax(demcopy)/10.0))*10
        step        = 50 # 50
        edges       = np.arange(emin,emax+step,step)
        ixd         = np.digitize(demcopy,edges) 
        hypsom      = pd.DataFrame(index = edges, columns = self.masks.keys())  

        for name in self.masks:
            elevbin     = np.multiply(ixd,self.masks[name]['mask'])

            # Do it by elevation band
            for n in np.arange(0,len(edges)):
                ind        = elevbin == n    
                hypsom.loc[edges[n],name] = (np.nansum(np.nansum(ind))/np.nansum(np.nansum(self.masks[name]['mask'])))*100
            
        
        for name in self.masks:
            smin = np.nanmin(self.dem[self.dem*self.masks[name]['mask'] > 0])
            smax = np.nanmax(self.dem[self.dem*self.masks[name]['mask'] > 0])
            print(name,smin,smax)
        
        # filter out high and low hyp
        map     = plt.cm.terrain
        map.set_bad('white')
        
        sns.set_style('darkgrid')
        sns.set_context("notebook")
       
        clrs = copy.copy(self.barcolors)
        clrs.insert(0,'k')
       
        plt.close(8) 
        fig,(ax,ax1)    = plt.subplots(num=8, figsize=self.figsize, dpi=self.dpi, nrows = 1, ncols = 2)
        h               = ax.imshow(demcopy, cmap=map)
        
        # Basin boundaries
        for name in self.masks:
            ax.contour(self.masks[name]['mask'],cmap = "Greys",linewidths = 1)   

        for iters,name in enumerate(self.plotorder):
            hyp = copy.deepcopy(hypsom[name].cumsum()-hypsom[name].iloc[0])
            imin = hyp < 0.1
            hyp[imin] = np.nan
            imax = hyp > 99.9
            hyp[imax] = np.nan
            
            # ax1.plot(hypsom[name].cumsum()[::-1]-hypsom[name].iloc[0],range(0,len(edges)),color = clrs[iters],label = name) 
            ax1.plot(hyp,range(0,len(edges)),color = clrs[iters],label = name) 
         
        ax1.invert_xaxis()
        ax1.legend()
        ax1.set_xlim((102,-2))
        ax1.set_ylim((-1,len(hypsom)))
        
        ax1.set_ylabel('elevation [m]')
        ax1.yaxis.tick_right()
        ax1.yaxis.set_label_position("right")
        
        labels = ax1.get_xticklabels()
        ax1.set_xticklabels(labels[::-1])
              
        xts         = ax1.get_yticks()
        edges_lbl   = []
        for i in xts[0:len(xts)-1]:
            edges_lbl.append(str(int(edges[int(i)])))        
        
        ax1.set_yticklabels(str(i) for i in edges_lbl)
        
        ax1.set_xlabel(r'Basin area above elevation [%]')
        ax.set_title(self.plotorder[0])    

        h.axes.get_xaxis().set_ticks([])
        h.axes.get_yaxis().set_ticks([])
        divider = make_axes_locatable(ax)
        cax     = divider.append_axes("left", size="5%", pad=0.2)
        cbar    = plt.colorbar(h, cax = cax)
        cax.yaxis.set_ticks_position('left')
        cax.yaxis.set_label_position('left')
        cax.set_ylabel('elevation [m]')
        
        plt.subplots_adjust(top=0.88)
        plt.tight_layout() 
        
        if self.basin == 'SJ':
            ax.text(0.3,0.88,'Main',horizontalalignment='center',transform=ax.transAxes,fontsize = 10)  
            ax.text(0.48,0.055,'Jose Creek',horizontalalignment='center',transform=ax.transAxes,fontsize = 10)
            ax.text(0.95,0.65,'South\nFork',horizontalalignment='center',transform=ax.transAxes,fontsize = 10)
            ax.text(0.05,0.55,'Willow\nCreek',horizontalalignment='center',transform=ax.transAxes,fontsize = 10)  
        if self.basin == 'BRB':
            ax.text(0.2,0.74,'Mores Creek',horizontalalignment='center',transform=ax.transAxes,fontsize = 10)  
            ax.text(0.48,0.77,'Twin\nSprings',horizontalalignment='center',transform=ax.transAxes,fontsize = 10)
            ax.text(0.88,0.6,'Featherville',horizontalalignment='center',transform=ax.transAxes,fontsize = 10) 
        if self.basin == 'Tuol':
            ax.text(0.75,0.68,'Tuolumne',horizontalalignment='center',transform=ax.transAxes,fontsize = 10)  
            ax.text(0.2,0.84,'Cherry',horizontalalignment='center',transform=ax.transAxes,fontsize = 10)
            ax.text(0.14,0.35,'Eleanor',horizontalalignment='center',transform=ax.transAxes,fontsize = 10)                                 
        
        print('saving figure to %shypsometry%s.png'%(self.figs_path,self.name_append))   
        plt.savefig('%shypsometry%s.png'%(self.figs_path,self.name_append))   
    
    def write_summary(self,df):

        # df          = 'accum_byelev'
        
        # dataframes  = ['accum_byelev','state_byelev','delta_state_byelev','melt','nonmelt','snowmelt_byelev','state_summary','accum_summary']
        dataframes  = ['accum_summary','state_summary']
        
        if df not in dataframes:
            print('df needs to be one of: %s'%(dataframes)) 
            # break
        
        print('Writing summary file to %s%s_summary.csv'%(self.figs_path,df))   
         
        if df == 'state_summary':
            self.state_summary.to_csv('%s%s_summary.csv'%(self.figs_path,df))
        if df == 'accum_summary':
            self.accum_summary.to_csv('%s%s_summary.csv'%(self.figs_path,df))            
        