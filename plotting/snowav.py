
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
import colorcet as cc

class snowav(object):

    def __init__(self,config_file):
        '''
        Notes: 

        '''
        
        try:
            print('Reading the config file...')
            cfg = cfp.ConfigParser()
            cfg.read(config_file)
    
            ####################################################
            #             Basin section                        #
            ####################################################
            self.basin = cfg.get('Basin','basin')
            self.save_path = cfg.get('Basin','save_path')
            self.name_append = cfg.get('Basin','name_append')
            self.wy = int(cfg.get('Basin','wy'))
            self.units = cfg.get('Basin','units') 
            
            ####################################################
            #           Outputs section                        #
            ####################################################        
            # Default for snow.0000 file is 2 (SWE)
            if (cfg.has_option('Outputs','snowband')): 
                self.snowband = int(cfg.get('Outputs','snowband'))     
            else:
                self.snowband = 2
            
            # Default for em.0000 is 8 (SWI)
            if (cfg.has_option('Outputs','emband')): 
                self.snowband = int(cfg.get('Outputs','emband'))     
            else:
                self.emband = 8           

            # Default rounding decimals for volume and depth 
            if (cfg.has_option('Outputs','decimals')): 
                self.dplcs = int(cfg.get('Outputs','decimals'))   
            else:
                self.dplcs = 1   
             
            # If psnowFile and csnowFile are specified, use those, 
            # otherwise they will get defined as the first and last
            # snow.XXXX file once [Runs] has been compiled
            if (cfg.has_option('Outputs','psnowFile') and 
                cfg.has_option('Outputs','csnowFile')):   
                self.psnowFile = cfg.get('Outputs','psnowFile')
                self.csnowFile = cfg.get('Outputs','csnowFile') 
                self.cemFile = self.csnowFile.replace('snow.','em.')  
                
            # self.file_type = cfg.get('Outputs','file_type')         
            
            ####################################################
            #           Runs                                   #
            ####################################################       
            # Collect all the run directories
            runs = list(cfg.items('Runs'))             
            self.snow_files = []
            self.em_files = []
            
            for rdir in runs:

                run_files = [rdir[1] + s for s in sorted(os.listdir(rdir[1]))]
                
                self.snow_files = (self.snow_files 
                                   + [value for value in run_files 
                                   if 'snow.' in value])
                self.em_files = (self.em_files 
                                 + [value for value in run_files 
                                 if 'em.' in value])
       
            # If no psnowFile and csnowFile specified, use first and last
            if (not cfg.has_option('Outputs','csnowFile') 
                and cfg.has_option('Outputs','psnowFile')):
                self.psnowFile = self.snow_files[0] 
                self.csnowFile = self.snow_files[len(self.snow_files)-1] 
                self.cemFile = self.em_files[len(self.em_files)-1]
                
            # Check to see if they exist
            if not (os.path.isfile(self.psnowFile) 
                    and os.path.isfile(self.psnowFile)):   
                print('psnowFile and/or csnowFile do not exist!')
                return
            
            ####################################################
            #           Validate                               #
            #################################################### 
            self.val_stns = cfg.get('Validate','stations').split(',')
            self.val_lbls = cfg.get('Validate','labels').split(',')
            self.val_client = cfg.get('Validate','client')                  
            
            ####################################################
            #           Accumulated                            #
            ####################################################   
            self.acc_clmin = cfg.get('Accumulated','clmin')
            self.acc_clmax = cfg.get('Accumulated','clmax')
            
            if (cfg.has_option('Accumulated','ymin') 
                and cfg.has_option('Accumulated','ymax')):
                self.acc_ylims = (int(cfg.get('Accumulated','ymin')),
                                  int(cfg.get('Accumulated','ymax')))        
    
            ####################################################
            #           Elevation                              #
            ####################################################             
            if (cfg.has_option('Elevation','ymin') 
                and cfg.has_option('Elevation','ymax')):
                self.el_ylims = (int(cfg.get('Elevation','ymin')),
                                 int(cfg.get('Elevation','ymax')))                      
            
            ####################################################
            #           Changes                                #
            ####################################################           
            self.ch_clmin = cfg.get('Changes','clmin')
            self.ch_clmax = cfg.get('Changes','clmax')
                 
            if (cfg.has_option('Changes','clminabs') 
                and cfg.has_option('Changes','clmaxabs')):   
                self.ch_clminabs = int(cfg.get('Changes','clminabs'))
                self.ch_clmaxabs = int(cfg.get('Changes','clmaxabs'))
            if (cfg.has_option('Changes','ymin') 
                and cfg.has_option('Changes','ymax')):
                self.ch_ylims = (int(cfg.get('Changes','ymin')),
                                 int(cfg.get('Changes','ymax')))            
            
            ####################################################
            #           Results                                #
            ####################################################           
            self.r_clmin = cfg.get('Results','clmin')
            self.r_clmax = cfg.get('Results','clmax')     
            
            ####################################################
            #           DEM                                    #
            ####################################################         
            self.demPath = cfg.get('DEM','demPath')
            self.total = cfg.get('DEM','total')     
    
            ####################################################
            #          Plots                                   #
            ####################################################  
            self.figsize = (int(cfg.get('Plots','fig_length')),
                            int(cfg.get('Plots','fig_height')))
            self.dpi = int(cfg.get('Plots','dpi'))
            self.barcolors = ['xkcd:true green','palegreen',
                              'xkcd:dusty green','xkcd:vibrant green','red']     
    
            ####################################################
            #          Report                                  #
            ####################################################        
            self.env_path = cfg.get('Report','env_path')
            self.tex_file = cfg.get('Report','tex_file')
            self.rep_path = cfg.get('Report','rep_path')
            self.templ_path = cfg.get('Report','templ_path')
            
            if cfg.has_option('Report','orig_date'):
                self.orig_date = cfg.get('Report','orig_date')
            if cfg.has_option('Report','report_name_append'):
                self.rep_append = cfg.get('Report','report_name_append')                
            
            # These will later get appended with self.dateTo 
            self.report_name = cfg.get('Report','report_name')
            self.rep_title = cfg.get('Report','report_title')          
              
            # Strings for the report
            if self.units == 'KAF':
                self.reportunits = 'KAF' 
            if self.units == 'SI':
                self.reportunits = 'km^3'
                           
            ####################################################  
            # History forecast
            ####################################################  
            if cfg.has_option('Hx Forecast','adj_hours'):
                self.adj_hours = int(cfg.get('Hx Forecast','adj_hours'))
            
            ####################################################   
            
            # Get the DEM 
            # There are different formats, this will get fixed once we
            # start using netcdf
            try: 
                self.dem = np.genfromtxt(self.demPath)
            except:
                self.dem = np.genfromtxt(self.demPath,skip_header = 6)
            
            self.nrows = len(self.dem[:,0])
            self.ncols = len(self.dem[0,:]) 
            blank = np.zeros((self.nrows,self.ncols))
                           
            # Assign some basin-specific things
            if self.basin == 'BRB':
                self.pixel = 100
                sr = 0
                if self.units == 'KAF':
                    emin = 2500      
                    emax = 10500     
                    self.step = 1000
                if self.units == 'SI':
                    emin = 800       
                    emax = 3200       
                    self.step = 250                 
            if self.basin == 'TUOL':
                self.pixel = 50   
                sr = 6
                if self.units == 'KAF':
                    emin = 3000      # [ft]
                    emax = 12000    
                    self.step   = 1000
                if self.units == 'SI':
                    emin = 800       # [m]
                    emax = 3600     
                    self.step   = 500                  
            if self.basin == 'SJ':
                self.pixel = 50
                sr = 6
                if self.units == 'KAF':
                    emin = 1000      # [ft]
                    emax  = 13000   
                    self.step = 1000
                if self.units == 'SI':
                    emin = 300       # [m]
                    emax = 4000       
                    self.step = 500                 
            if self.basin == 'LAKES':
                self.pixel = 50  
                sr = 0
                self.imgx = (1200,1375)
                self.imgy = (425,225) 
                if self.units == 'KAF':
                    emin = 8000      # [ft]
                    emax = 12500      
                    self.step = 500
                if self.units == 'SI':
                    emin = 2400       # [m]
                    emax = 3800      
                    self.step = 200                    
            if self.basin == 'RCEW':
                self.pixel = 50
                sr = 0
                if self.units == 'KAF':
                    emin = 2500      # [ft]
                    emax = 7500       
                    self.step = 500
                if self.units == 'SI':
                    emin = 800       # [m]
                    emax = 2500     
                    self.step = 500  
     
            # This is for creating the elevation bins
            self.edges = np.arange(emin,emax+self.step,self.step)
            # Right now this is a placeholder, could edit by basin...
            self.xlims = (0,len(self.edges))
            
            self.subbasin1 = cfg.get('DEM','subbasin1')
            self.subbasin2 = cfg.get('DEM','subbasin2')
            self.subbasin3 = cfg.get('DEM','subbasin3')
            self.total_lbl = cfg.get('DEM','total_lbl')
            self.sub1_lbl = cfg.get('DEM','sub1_lbl')
            self.sub2_lbl = cfg.get('DEM','sub2_lbl')
            self.sub3_lbl = cfg.get('DEM','sub3_lbl') 
            
            self.plotorder = [self.total_lbl, self.sub1_lbl, 
                              self.sub2_lbl, self.sub3_lbl]  
            self.suborder = [self.sub1_lbl,self.sub2_lbl,self.sub3_lbl]  
            maskpaths = [self.total, self.subbasin1, 
                         self.subbasin2,self.subbasin3 ]

            # Add if necessary - need to generalize all this and change 
            # in the config file!
            if cfg.has_option('DEM','sub4_lbl'): 
                self.sub4_lbl = cfg.get('DEM','sub4_lbl') 
                self.subbasin4 = cfg.get('DEM','subbasin4')
                self.plotorder = self.plotorder + [self.sub4_lbl]   
                self.suborder = self.suborder + [self.sub4_lbl]                        
                maskpaths = maskpaths + [self.subbasin4]
            
            # Compile the masks 
            self.masks = dict()
            for lbl,mask in zip(self.plotorder,maskpaths):
                if lbl != self.sub4_lbl:
                    self.masks[lbl] = {'border': blank, 
                                       'mask': np.genfromtxt(mask,skip_header=sr),
                                       'label': lbl}   
                else:
                    self.masks[lbl] = {'border': blank, 
                                       'mask': np.genfromtxt(mask,skip_header=0),
                                       'label': lbl}                               
            
            # Do unit-specific things
            if self.units == 'KAF':
                self.conversion_factor = ((self.pixel**2)
                                         * 0.000000810713194*0.001) # [KAF]   
                self.depth_factor = 0.03937 # [inches]
                self.dem = self.dem * 3.28 # [ft]
                self.ixd = np.digitize(self.dem,self.edges)     
                self.depthlbl = 'in'
                self.vollbl = 'KAF' 
                self.elevlbl = 'ft'   
            
            if self.units == 'SI':
                # [km^3]
                self.conversion_factor = ((self.pixel**2)
                                          * 0.000000810713194*1233.48/1e9) 
                self.depth_factor = 1 # [m]                                                
                self.ixd = np.digitize(self.dem,self.edges)     
                self.depthlbl = 'mm'
                self.vollbl = '$km^3$' 
                self.elevlbl = 'm'     
                                 
            # Copy the config file where figs will be saved
            extf = os.path.splitext(os.path.split(config_file)[1])
            path_shr = os.path.split(self.psnowFile)
            path_ehr = os.path.split(self.csnowFile)
            ext_shr = os.path.splitext(path_shr[1])
            ext_ehr = os.path.splitext(path_ehr[1])
            self.figs_path = (self.save_path 
                             + '%s_%s/'%(ext_shr[1][1:5],ext_ehr[1][1:5]))
            
            # Check run dates, and make a new folder if necessary
            if not os.path.exists(self.figs_path):
                os.makedirs(self.figs_path)
    
            # Only need to store this name if we decide to 
            # write more to the copied config file...
            self.config_copy = (self.figs_path 
                                + extf[0]
                                + self.name_append
                                + '_%s_%s'%(ext_shr[1][1:5],ext_ehr[1][1:5]) 
                                + extf[1])
            
            # If it doesn't already exist, make it
            if not os.path.isfile(self.config_copy): 
                copyfile(config_file,self.config_copy)
                       
            print('done.')    
        
        except:
            print('error reading config file.')

    def process(self,*args):
        '''
        This function calculates everything we will need for the plots.
        
        Does not currently save to csv...

        '''
  
        print('Processing iSnobal outputs...')
        
        # If we add pre and current snow files for flights, force those
        # and add a few new things to calculate
        if len(args) != 0:
            self.psnowFile = args[0]
            self.csnowFile = args[1]
            self.cemFile = args[1].replace('snow','em')
            self.snow_files = [self.psnowFile, self.csnowFile]
            self.em_files = [self.psnowFile.replace('snow','em'), 
                             self.csnowFile.replace('snow','em')]
       
        cclimit = -5*1000*1000  # based on an average of 60 W/m2 from TL paper
        # ccM = cc./1000./1000; % cold content in MJ

        accum = np.zeros((self.nrows,self.ncols))
        accum_sub = copy.deepcopy(accum)
        snowmelt = copy.deepcopy(accum)
        rain_bg = copy.deepcopy(accum)
        precip = copy.deepcopy(accum)
        precip_sub = copy.deepcopy(accum)
        snowmelt_sub =copy.deepcopy(accum)
        state = copy.deepcopy(accum)
        depth = copy.deepcopy(accum)
        pstate = copy.deepcopy(accum)
        state_byday = np.zeros((self.nrows,self.ncols,len(self.snow_files)))
        
        accum_byelev = pd.DataFrame(index = self.edges, 
                                    columns = self.masks.keys())
        rain_bg_byelev = copy.deepcopy(accum_byelev)
        precip_byelev = copy.deepcopy(accum_byelev)
        precip_byelev_sub = copy.deepcopy(accum_byelev)
        accum_byelev_sub = copy.deepcopy(accum_byelev)  
        state_byelev = copy.deepcopy(accum_byelev)
        density_m_byelev = copy.deepcopy(accum_byelev)
        depth_byelev = copy.deepcopy(accum_byelev)
        state_mswe_byelev = copy.deepcopy(accum_byelev)
        depth_mdep_byelev = copy.deepcopy(accum_byelev)
        delta_state_byelev = copy.deepcopy(accum_byelev)
        melt = copy.deepcopy(accum_byelev)
        nonmelt = copy.deepcopy(accum_byelev)
        snowmelt_byelev = copy.deepcopy(accum_byelev)
        snowmelt_byelev_sub = copy.deepcopy(accum_byelev)
        
        state_summary = pd.DataFrame(columns = self.masks.keys())
        accum_summary = pd.DataFrame(columns = self.masks.keys())
        precip_summary = pd.DataFrame(columns = self.masks.keys())
               
        # Loop through output files
        # Currently we need to load in each em.XXXX file to 'accumulate',
        # but only the first and last snow.XXXX file for changes
        accum_sub_flag = False
        adj = 0 # this is a hack for Hx-repeats-itself forecasting
        t = 0
        for iters,(em_name,snow_name) in enumerate(zip(self.em_files,
                                                       self.snow_files)):
            # iters = 0
            # iters = iters + 1
            # em_name = self.em_files[iters]
            # snow_name = self.snow_files[iters]        

            date = wy.wyhr_to_datetime(self.wy,
                                       int(snow_name.split('.')[-1])
                                       + adj)  
            
            # It's nice to see where we are...
            print(snow_name, int(snow_name.split('.')[-1]) - t )
            print(date)
            
            # Hack for Hx-repeats-itself forecasting
            if snow_name == self.psnowFile:
                if hasattr(self,"adj_hours"):
                    print('Hacking adj_hours...')
                    adj = self.adj_hours
                   
            em_file = ipw.IPW(em_name)
            band = em_file.bands[self.emband].data
            accum = accum + band
            snowmelt = snowmelt + em_file.bands[7].data
 
            # load and calculate sub-basin total
            snow_file = ipw.IPW(snow_name)
            tmpstate = snow_file.bands[self.snowband].data  
            state_byday[:,:,iters] = tmpstate
            
            # Get rain from input data
            sf = em_name.replace('runs','data')
            sf = sf.replace('run','data')
            sf = sf.replace('output','ppt_4b')
            ppt_path = sf.split('em')[0]
            
            out_hr = snow_name.split('.')[-1]
            hrs = range(int(out_hr) - 23,int(out_hr) + 1)
            
            pFlag = False
            ppt_files = []
            for hr in hrs:
                if os.path.isfile(ppt_path +'ppt.4b_'+ str(hr)):
                    pFlag = True
                    ppt_files = ppt_files + [ppt_path +'ppt.4b_'+ str(hr)]
            
            rain_hrly = np.zeros((self.nrows,self.ncols))
            precip_hrly = np.zeros((self.nrows,self.ncols))
            
            # Load 'em in
            if pFlag:
                for pfile in ppt_files:
                    ppt = ipw.IPW(pfile)
                    pre = ppt.bands[0].data
                    percent_snow = ppt.bands[1].data
                    rain_hrly = rain_hrly + np.multiply(pre,(1-percent_snow))
                    precip_hrly = precip_hrly + pre

            rain_bg = rain_bg + rain_hrly 
            precip = precip + precip_hrly    
                
            # Currently this is grabbing the second to last
            # Being used for flight difference
            # Definitely a better way to do this...
            if iters == (len(self.snow_files) - 1):
                # Mean pixel depth [mm] on psnowFile
                self.pre_pm = (
                               np.nansum(tmpstate 
                               * self.masks[self.total_lbl]['mask'])
                               / self.masks[self.total_lbl]['mask'].sum() 
                               )  
                self.pre_swe = ( 
                                np.nansum(tmpstate 
                                * self.masks[self.total_lbl]['mask'])
                                * self.conversion_factor )      
            
            # Store daily sub-basin totals
            for name in self.masks:
                accum_summary.loc[date, name] = (np.nansum(
                                                np.multiply(accum,
                                                self.masks[name]['mask'])) )
                state_summary.loc[date, name] = (np.nansum(
                                                np.multiply(tmpstate,
                                                self.masks[name]['mask'])) )
                precip_summary.loc[date, name] = (np.nansum(
                                                np.multiply(precip,
                                                self.masks[name]['mask'])) )

            # This only add between psnowFile and csnowFile
            if accum_sub_flag:
                accum_sub = accum_sub + copy.deepcopy(band)
                snowmelt_sub = ( snowmelt_sub 
                                + copy.deepcopy(em_file.bands[7].data) ) 
                precip_sub = precip_sub + precip_hrly            
           
            # When it is the first snow file, copy
            if snow_name == self.psnowFile:
                print('psnowfile is %s'%(snow_name))
                accum_sub_flag = True
                pstate = copy.deepcopy(tmpstate)           
            
            # When it hits the current snow file, copy
            if snow_name == self.csnowFile:
                print('csnowfile is %s'%(snow_name))
                
                # Turn off, but last one will still have been added
                accum_sub_flag  = False
             
                state = copy.deepcopy(tmpstate)
                depth = snow_file.bands[0].data
                density = snow_file.bands[1].data
                self.density = copy.deepcopy(density)
                self.cold = em_file.bands[9].data   
                
                # No need to compile more files after csnowFile
                break     

            t = int(snow_name.split('.')[-1])

        self.dateFrom = wy.wyhr_to_datetime(self.wy,
                                            int(self.psnowFile.split('.')[-1]))
        self.dateTo = wy.wyhr_to_datetime(self.wy,
                                          int(self.csnowFile.split('.')[-1]) 
                                          + adj)            
 
        # Append date to report name
        parts = self.report_name.split('.')
        self.report_name = ( parts[0] 
                           + self.dateTo.date().strftime("%Y%m%d") 
                           + '.' + parts[1] )
        
        # Difference in state (SWE)
        delta_state         = state - pstate
        
        # Mask by subbasin and elevation band
        for name in self.masks:
            mask = copy.deepcopy(self.masks[name]['mask'])
            
            accum_mask = np.multiply(accum,mask)
            rain_bg_mask = np.multiply(rain_bg,mask)
            precip_mask = np.multiply(precip,mask)
            precip_sub_mask = np.multiply(precip_sub,mask)
            accum_mask_sub = np.multiply(accum_sub,mask)
            snowmelt_mask = np.multiply(snowmelt,mask)
            snowmelt_mask_sub = np.multiply(snowmelt_sub,mask)
            state_mask = np.multiply(state,mask)
            delta_state_byelev_mask = np.multiply(delta_state,mask)
            state_byelev_mask = np.multiply(state,mask)
            state_mswe_byelev_mask = np.multiply(state,mask)
            density_m_byelev_mask = np.multiply(density,mask)
            
            ix = density_m_byelev_mask == 0
            density_m_byelev_mask[ix] = np.nan
            
            depth_mdep_byelev_mask = np.multiply(depth,mask)
            elevbin = np.multiply(self.ixd,mask)
            
            # Do it by elevation band
            for n in np.arange(0,len(self.edges)):
                ind = elevbin == n
                state_bin = state_mask[ind]
                b = self.edges[n]
                
                # Cold content
                ccb = self.cold[ind]
                cind = ccb > cclimit
                
                accum_byelev.loc[b,name] = np.nansum(accum_mask[ind])
                rain_bg_byelev.loc[b,name] = np.nansum(rain_bg_mask[ind])
                precip_byelev.loc[b,name] = np.nansum(precip_mask[ind])
                precip_byelev_sub.loc[b,name] = np.nansum(precip_sub_mask[ind])
                accum_byelev_sub.loc[b,name] = np.nansum(accum_mask_sub[ind])
                snowmelt_byelev.loc[b,name] = np.nansum(snowmelt_mask[ind])
                state_byelev.loc[b,name] = np.nansum(state_byelev_mask[ind])
                depth_byelev.loc[b,name] = np.nansum(state_byelev_mask[ind])  
                melt.loc[b,name] = np.nansum(state_bin[cind])
                nonmelt.loc[b,name] = np.nansum(state_bin[~cind])   
                snowmelt_byelev_sub.loc[b,name] = np.nansum(
                                                snowmelt_mask_sub[ind])
                
                
                # Calculate mean if there is snow 
                if state_mswe_byelev_mask[ind].size:
                    state_mswe_byelev.loc[b,name] = np.nanmean(
                                                state_mswe_byelev_mask[ind])
                    depth_mdep_byelev.loc[b,name] = np.nanmean(
                                                depth_mdep_byelev_mask[ind])
                    delta_state_byelev.loc[b,name] = np.nansum(
                                                delta_state_byelev_mask[ind])
                else:
                    state_mswe_byelev.loc[b,name] = np.nan
                    depth_mdep_byelev.loc[b,name] = np.nan
                    delta_state_byelev.loc[b,name] = np.nan
                
            self.masks[name]['SWE'] = ( (melt[name].sum() + nonmelt[name].sum())
                                       * self.conversion_factor )      
        
        # Convert to desired units
        
        # Sum over all time steps, spatial
        self.precip = np.multiply(precip,self.depth_factor)  
        self.accum = np.multiply(accum,self.depth_factor)     
        self.accum_sub = np.multiply(accum_sub,self.depth_factor)   
        
        # At last time step
        self.state_byelev = np.multiply(state_byelev,self.conversion_factor)
        self.depth_byelev = np.multiply(depth_byelev,self.conversion_factor)         
        self.state = np.multiply(state,self.depth_factor) 
        self.depth = np.multiply(np.multiply(depth,1000),self.depth_factor) 
        self.delta_state = np.multiply(delta_state,self.depth_factor)  
        
        # Sum, by elevation
        self.accum_byelev = np.multiply(accum_byelev,self.conversion_factor) 
        self.precip_byelev = np.multiply(precip_byelev,self.conversion_factor) 
        self.rain_bg_byelev = np.multiply(rain_bg_byelev,
                                          self.conversion_factor)
        self.snowmelt_byelev = np.multiply(snowmelt_byelev,
                                           self.conversion_factor)  
        
        # Sub-set of time defined by psnowfile and csnowfile, sum, by elevation
        self.snowmelt_byelev_sub = np.multiply(snowmelt_byelev_sub,
                                               self.conversion_factor)  
        self.accum_byelev_sub = np.multiply(accum_byelev_sub,
                                            self.conversion_factor)      
        self.precip_byelev_sub = np.multiply(precip_byelev_sub,
                                             self.conversion_factor)
        
        # Change, by elevation
        self.delta_state_byelev = np.multiply(delta_state_byelev,
                                              self.conversion_factor) 
        
        # At first time step, spatial
        self.pstate = np.multiply(pstate,self.conversion_factor)
        
        # Daily 
        self.state_summary = np.multiply(state_summary,self.conversion_factor)       
        self.accum_summary = np.multiply(accum_summary,self.conversion_factor)         
        self.precip_summary = np.multiply(precip_summary,self.conversion_factor)                                
        self.state_byday = np.multiply(state_byday,self.depth_factor)
        self.state_summary = np.multiply(state_summary,self.conversion_factor)        
        self.accum_summary = np.multiply(accum_summary,self.conversion_factor)     
        self.precip_summary = np.multiply(precip_summary,self.conversion_factor)   
    
        # Mean      
        self.state_mswe_byelev = np.multiply(state_mswe_byelev,
                                             self.depth_factor)
        self.depth_mdep_byelev = np.multiply(np.multiply(depth_mdep_byelev,1000)
                                             ,self.depth_factor)         
        self.density_m_byelev = density_m_byelev
        
        self.melt = np.multiply(melt,self.conversion_factor)
        self.nonmelt = np.multiply(nonmelt,self.conversion_factor)               
        self.cold = np.multiply(self.cold,0.000001)                                    
        
        # Let's print some summary info...
        mask        = self.state_summary.index.to_series().diff() > pd.Timedelta('24:10:00')
        msum        = sum(mask)
        print('%s entries in dataframe index with gaps larger than 24h '%(msum))
        print('done.')
       
    def accumulated(self,*arg):
        '''
        This function plots self.accum (typically SWI)
        
        '''
        
        # This is outdated and probably needs to be re-done...
        # If this is supplied, only report accum by the subsection between 
        # psnowFile and csnowFile 
        if arg == 'sub':
            accum = copy.deepcopy(self.accum_sub) 
            accum_byelev = copy.deepcopy(self.accum_byelev_sub) 
                                       
        else:    
            accum = copy.deepcopy(self.accum)
            accum_byelev = copy.deepcopy(self.accum_byelev)
            
        qMin,qMax = np.percentile(accum,[0,99.8])
        clims = (0,qMax)
        colors1 = cmocean.cm.dense(np.linspace(0., 1, 255))
        colors2 = plt.cm.binary(np.linspace(0, 1, 1))
        colors = np.vstack((colors2, colors1))
        mymap = mcolors.LinearSegmentedColormap.from_list('my_colormap', colors)

        sns.set_style('darkgrid')
        sns.set_context("notebook")
        
        # White background
        pmask = self.masks[self.total_lbl]['mask']
        ixo = pmask == 0
        accum[ixo] = np.nan
        mymap.set_bad('white',1.) 
        
        # Now set SWI-free to some color
        ixf = accum == 0
        accum[ixf] = -1
        mymap.set_under('grey',1.) 
        
        plt.close(0)
        fig,(ax,ax1) = plt.subplots(num=0, 
                                    figsize = self.figsize, 
                                    dpi=self.dpi, 
                                    nrows = 1, ncols = 2)      
        h = ax.imshow(accum, clim = clims, cmap = mymap)
        
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
        
        cax = divider.append_axes("right", size="4%", pad=0.2)
        cbar = plt.colorbar(h, cax = cax)
        # cbar.ax.tick_params() 
        cbar.set_label('[%s]'%(self.depthlbl))
 
        h.axes.set_title('Accumulated SWI \n %s to %s'
                         %(self.dateFrom.date().strftime("%Y-%-m-%-d"),
                           self.dateTo.date().strftime("%Y-%-m-%-d")))  
        
        # Total basin label
        sumorder = self.plotorder[1:]
        if self.basin == 'LAKES' or self.basin == 'RCEW':
            sumorder = [self.plotorder[0]]
            
        if self.dplcs == 0:
            tlbl = '%s = %s %s'%(self.plotorder[0],
                                 str(int(accum_byelev[self.plotorder[0]].sum())),
                                 self.vollbl)
        else: 
            tlbl = '%s = %s %s'%(self.plotorder[0],
                                 str(np.round(accum_byelev[self.plotorder[0]].sum(),
                                self.dplcs)),self.vollbl)
              
        # Plot the bars
        for iters,name in enumerate(sumorder):
            if self.dplcs == 0:           
                lbl = '%s = %s %s'%(name,str(int(accum_byelev[name].sum())),
                                    self.vollbl)
            else:
                lbl = '%s = %s %s'%(name,str(np.round(accum_byelev[name].sum(),
                                    self.dplcs)),self.vollbl)

            if iters == 0:
                ax1.bar(range(0,len(self.edges)),accum_byelev[name], 
                        color = self.barcolors[iters], 
                        edgecolor = 'k',label = lbl)
            elif iters == 1:   
                ax1.bar(range(0,len(self.edges)),accum_byelev[name], 
                        bottom = accum_byelev[sumorder[iters-1]], 
                        color = self.barcolors[iters], edgecolor = 'k',label = lbl)
              
            elif iters == 2:   
                ax1.bar(range(0,len(self.edges)),accum_byelev[name], 
                        bottom = (accum_byelev[sumorder[iters-1]] + accum_byelev[sumorder[iters-2]]), 
                        color = self.barcolors[iters], edgecolor = 'k',label = lbl)

            elif iters == 3:   
                ax1.bar(range(0,len(self.edges)),accum_byelev[name], 
                        bottom = (accum_byelev[sumorder[iters-1]] + accum_byelev[sumorder[iters-2]] + accum_byelev[sumorder[iters-3]]), 
                        color = self.barcolors[iters], edgecolor = 'k',label = lbl)
                
               
            plt.rcParams['hatch.linewidth'] = 1
            plt.rcParams['hatch.color'] = 'k'                
            ax1.set_xlim((self.xlims[0]-0.5,self.xlims[1]))

        
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
        
        ax1.set_ylim((0,ylims[1] + ylims[1]*0.2))

        plt.tight_layout()
        fig.subplots_adjust(top=0.88)
      
        if self.basin != 'LAKES' and self.basin != 'RCEW':
            # more ifs for number subs...
            if len(self.plotorder) == 5:
                ax1.legend(loc= (0.01,0.68))
            elif len(self.plotorder) == 4:
                ax1.legend(loc= (0.01,0.74))
            
        if self.basin == 'BRB':
            ax1.text(0.26,0.94,tlbl,horizontalalignment='center',
                     transform=ax1.transAxes,fontsize = 10)
        else:
            ax1.text(0.3,0.94,tlbl,horizontalalignment='center',
                     transform=ax1.transAxes,fontsize = 10)
        
        # Make SWI-free legend if we need one
        if sum(sum(ixf)) > 1000:
            patches = [mpatches.Patch(color='grey', label='no SWI')]
            if self.basin == 'SJ':
                ax.legend(handles=patches, bbox_to_anchor=(0.3, 0.05), 
                          loc=2, borderaxespad=0. )
            else:
                ax.legend(handles=patches, bbox_to_anchor=(0.05, 0.05), 
                          loc=2, borderaxespad=0. )
             
        print('saving figure to %sswi%s.png'%(self.figs_path,self.name_append))
        plt.savefig('%sswi%s.png'%(self.figs_path,self.name_append))   
        
    def current_image(self):
        '''
        This plots self.state (typically SWE) and self.cold
        
        edits: make flexible input arguments for depth, density, etc
        
        '''
        # Make a copy so we can edit for plots

        state = copy.deepcopy(self.state)
        cold = copy.deepcopy(self.cold)

        qMin,qMax = np.nanpercentile(state,[0,99.9])    
        clims = (qMin,qMax)
        clims2 = (-5,0) 
        
        # Areas outside basin
        pmask = self.masks[self.total_lbl]['mask']
        ixo = pmask == 0   
        
        # Prepare no-snow and outside of the basin for the colormaps   
        ixz = state == 0
        state[ixz] = -1
        cold[ixz]  = 1     
        
        # Colormap for self.state
        colorsbad = plt.cm.Set2_r(np.linspace(0., 1, 1))
        colors1 = cmocean.cm.haline_r(np.linspace(0., 1, 254))
        colors = np.vstack((colorsbad,colors1))
        mymap = mcolors.LinearSegmentedColormap.from_list('my_colormap', colors) 
        
        state[ixo] = np.nan        
        mymap.set_bad('white',1.) 
        # mymap.set_under('lightgrey',1)    
        
        # Colormap for cold content
        # colorsbad       = plt.cm.binary(np.linspace(0., 1, 1))
        # colors1         = plt.cm.ocean_r(np.linspace(0., 1, 128))
        # colors2         = plt.cm.YlOrRd(np.linspace(0., 1, 127))
        # colors          = np.vstack((colorsbad, colors1, colors2))
        # mymap1          = mcolors.LinearSegmentedColormap.from_list('my_colormap', colors)
        # mymap1          = plt.cm.RdYlBu_r 
        # mymap1          = cc.m_diverging_linear_bjy_30_90_c45
        mymap1 = plt.cm.Spectral_r
        cold[ixo] = np.nan
        mymap1.set_bad('white')
        mymap1.set_over('lightgrey',1) 
                         
        sns.set_style('dark')
        sns.set_context("notebook")
        
        plt.close(1)
        fig,(ax,ax1) = plt.subplots(num=1, figsize=self.figsize, 
                                    facecolor = 'white', dpi=self.dpi, 
                                    nrows = 1, ncols = 2)
        h = ax.imshow(state, cmap = mymap, clim=clims)
        h1 = ax1.imshow(cold, clim=clims2, cmap = mymap1)
        
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
        cax = divider.append_axes("right", size="5%", pad=0.2)
        cbar = plt.colorbar(h, cax = cax)
        
        if self.units == 'KAF':         
            cbar.set_label('[in]')
        if self.units == 'SI':
            cbar.set_label('[mm]')
         
        # cbar.ax.tick_params()                
        
        # Do pretty stuff for the right plot
        h1.axes.get_xaxis().set_ticks([])
        h1.axes.get_yaxis().set_ticks([])
        h1.axes.set_title('Cold Content \n %s'%(self.dateTo.date().strftime("%Y-%-m-%-d")))
        divider = make_axes_locatable(ax1)
        cax2 = divider.append_axes("right", size="5%", pad=0.2)
        cbar1 = plt.colorbar(h1, cax = cax2)
        cbar1.set_label('Cold Content [MJ/$m^3$]')
        cbar1.ax.tick_params() 
        
        h.axes.set_title('SWE \n %s'%(self.dateTo.date().strftime("%Y-%-m-%-d")))
        fig.subplots_adjust(top=0.95,bottom=0.05,
                            right = 0.92, left = 0.05, wspace = 0.12)
        if self.basin == 'LAKES':
            plt.tight_layout()
        
        patches = [mpatches.Patch(color='grey', label='snow free')]
        if self.basin == 'SJ':
            ax.legend(handles=patches, bbox_to_anchor=(0.3, 0.05), loc=2, borderaxespad=0. )
        if self.basin == 'RCEW':
            ax.legend(handles=patches, bbox_to_anchor=(-0.2, 0.05), loc=2, borderaxespad=0. )            
        else:
            ax.legend(handles=patches, bbox_to_anchor=(0.05, 0.05), loc=2, borderaxespad=0. )
        
        print('saving figure to %sresults%s.png'%(self.figs_path,self.name_append))
        plt.savefig('%sresults%s.png'%(self.figs_path,self.name_append))  
        
    def pixel_swe(self):
        
        
        plt.close(2) 
        fig,(ax,ax1) = plt.subplots(num=2, figsize=self.figsize,
                                    dpi=self.dpi, nrows = 1, ncols = 2)
   
        sumorder = self.plotorder[1:]
        if self.basin == 'LAKES' or self.basin == 'RCEW':
            sumorder = [self.plotorder[0]]  
            swid = 0.45
        else:
            sumorder = self.plotorder[1::]
            swid = 0.25
        
        wid = np.linspace(-0.3,0.3,len(sumorder))
        
        for iters,name in enumerate(sumorder): 
            # iters = 0
            # name = sumorder[iters]                             
            ax.bar(range(0,len(self.edges))-wid[iters],
                    self.depth_mdep_byelev[name], 
                    color = self.barcolors[iters], width = swid, edgecolor = 'k',label = name) 
                        
            ax1.bar(range(0,len(self.edges))-wid[iters],
                    self.state_mswe_byelev[name], 
                    color = self.barcolors[iters], width = swid, edgecolor = 'k',label = name)  
            
        ax.set_xlim((0,len(self.edges)))    
        ax1.set_xlim((0,len(self.edges)))        
        
        plt.tight_layout()                            
        xts         = ax1.get_xticks()
        edges_lbl   = []
        for i in xts[0:len(xts)-1]:
            edges_lbl.append(str(int(self.edges[int(i)])))
  
        ax.set_xticklabels(str(i) for i in edges_lbl)
        ax1.set_xticklabels(str(i) for i in edges_lbl)
        for tick,tick1 in zip(ax.get_xticklabels(),ax1.get_xticklabels()):
            tick.set_rotation(30)
            tick1.set_rotation(30)  
             
        ylims = ax1.get_ylim()
        ax1.set_ylim((ylims[0],ylims[1]+ylims[1]*0.2))  
        
        ylims = ax.get_ylim()
        ax.set_ylim((ylims[0],ylims[1]+ylims[1]*0.2))               

        ax.set_ylabel('mean depth [%s]'%(self.depthlbl))
        ax.set_xlabel('elevation [ft]')
        ax1.set_ylabel('mean SWE [%s]'%(self.depthlbl))
        ax1.set_xlabel('elevation [ft]')
        ax.legend(loc='upper left')
        ax.grid('on')
        ax1.grid('on')
        
        ax1.yaxis.set_label_position("right")
        ax1.yaxis.tick_right()       
        
        ax.set_title('Mean Depth, %s'%(self.dateTo.date().strftime("%Y-%-m-%-d")))
        ax1.set_title('Mean SWE, %s'%(self.dateTo.date().strftime("%Y-%-m-%-d")))
        
        plt.tight_layout()
        print('saving figure to %smean_swe_depth%s.png'%(self.figs_path,self.name_append))
        plt.savefig('%smean_swe_depth%s.png'%(self.figs_path,self.name_append))   
        
    def density(self): 
        '''
        None of these are really finished...
        '''
        
        value = copy.deepcopy(self.density_m_byelev)
        lim = np.max(value[self.total_lbl]) 
        ylim = (0,600) 
        color = 'xkcd:windows blue'

        sns.set_style('darkgrid')
        sns.set_context("notebook")        

        nf = len(self.masks)
        
        if nf > 1 and nf < 5:
            nr = 2
            nc = 2
        elif nf < 7:
            nr = 2
            nc = 3
            
        plt.close(3)
        fig,ax  = plt.subplots(num=3, figsize=self.figsize, dpi=self.dpi, nrows = nr, ncols = nc)
        axs     = ax.ravel()
        if nf == 5:
            fig.delaxes(axs[5])
               
        for iters,name in enumerate(self.plotorder):
            # iters = 0
            # name = self.plotorder[iters]
            axs[iters].bar(range(0,len(self.edges)),value[name], color = color)
            
            axs[iters].set_xlim(self.xlims)                          
            xts         = axs[iters].get_xticks()
            if len(xts) < 6:
                dxt = xts[1] - xts[0]
                xts = np.arange(xts[0],xts[-1] + 1 ,dxt/2)
                axs[iters].set_xticks(xts)                 
    
            edges_lbl   = []
            for i in xts[0:len(xts)-1]:
                edges_lbl.append(str(int(self.edges[int(i)])))

            axs[iters].set_xticklabels(str(i) for i in edges_lbl)
            
            if iters == 0:
                axs[iters].set_ylabel(r'$\rho$ [kg/$m^3$]') 
            
            if iters > nc - 1:
                axs[iters].set_xlabel('elevation [%s]'%(self.elevlbl))
   
            # Put yaxis on right 
            if iters == nc - 1 or iters == nf -1 :
                axs[iters].yaxis.set_label_position("right")
                axs[iters].yaxis.tick_right()
                axs[iters].set_ylabel(r'$\rho$ [kg/$m^3$]') 
            
            if iters == 1 and nc == 3:
                axs[iters].set_yticklabels([])
                
            if iters <= nc - 1:
                axs[iters].set_xticklabels([])
            
            axs[iters].set_ylim((ylim)) 
            for tick in axs[iters].get_xticklabels():
                tick.set_rotation(30) 
                
            axs[iters].text(0.5,0.92,name,horizontalalignment='center',transform=axs[iters].transAxes,fontsize = 10)
        
        for n in range(0,len(axs)):
            axs[n].set_xticks(xts)
            axs[n].set_xlim(self.xlims)       

        fig.tight_layout()
        fig.subplots_adjust(top=0.92,wspace = 0.1)
        fig.suptitle(r'Density, %s'%(self.dateTo.date().strftime("%Y-%-m-%-d")) )       
        
        print('saving figure to %sdensity_subs%s.png'%(self.figs_path,self.name_append))
        plt.savefig('%sdensity_subs%s.png'%(self.figs_path,self.name_append))     
        
        ###############################
        # 2nd density fig
        ###############################

        cvalue = copy.deepcopy(self.density)
        
        colors1 = cmocean.cm.speed(np.linspace(0., 1, 255))
        colors2 = plt.cm.binary(np.linspace(0, 1, 1))
        colors = np.vstack((colors2, colors1))
        mymap = mcolors.LinearSegmentedColormap.from_list('my_colormap', colors)

        sns.set_style('darkgrid')
        sns.set_context("notebook")
        
        # This is to get the background white
        pmask = self.masks[self.total_lbl]['mask']
        ixo = pmask == 0
        cvalue[ixo] = np.nan
        mymap.set_bad('white',1.) 
        
        ixf = cvalue == 0
        cvalue[ixf] = -1
        mymap.set_under('lightgrey',1.) 
        
        plt.close(4)
        fig,(ax,ax1) = plt.subplots(num=4, figsize = self.figsize, 
                                    dpi=self.dpi, nrows = 1, ncols = 2)      
        h = ax.imshow(cvalue, clim = (50,550), interpolation='none', cmap = mymap)
        
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
        
        cax = divider.append_axes("right", size="4%", pad=0.2)
        cbar = plt.colorbar(h, cax = cax)
        cbar.ax.tick_params() 
        cbar.set_label('[kg/$m^3$]')
 
        h.axes.set_title('Density\n%s'%(self.dateTo.date().strftime("%Y-%-m-%-d")))  

        ax1.bar(range(0,len(self.edges)),value[self.total_lbl], 
                color = 'g', edgecolor = 'k')
      
        plt.rcParams['hatch.linewidth'] = 1
        plt.rcParams['hatch.color'] = 'k'                
        ax1.set_xlim((self.xlims[0]-0.5,self.xlims[1]))
        
        plt.tight_layout()
        xts = ax1.get_xticks()
        edges_lbl = []
        for i in xts[0:len(xts)-1]:
            edges_lbl.append(str(int(self.edges[int(i)])))
        
        ax1.set_xticklabels(str(i) for i in edges_lbl)
        for tick in ax1.get_xticklabels():
            tick.set_rotation(30)        

        ax1.set_ylabel('density - per elevation band')
        ax1.set_xlabel('elevation [%s]'%(self.elevlbl))
       
        ax1.yaxis.set_label_position("right")
        ax1.yaxis.tick_right()
        ax1.set_ylim((0,600))

        plt.tight_layout()
        fig.subplots_adjust(top=0.88)
    
        if sum(sum(ixf)) > 1000:
            patches = [mpatches.Patch(color='grey', label='snow free')]
            if self.basin == 'SJ':
                ax.legend(handles=patches, bbox_to_anchor=(0.3, 0.05), loc=2, borderaxespad=0. )
            else:
                ax.legend(handles=patches, bbox_to_anchor=(0.05, 0.05), loc=2, borderaxespad=0. )
        
        print('saving figure to %sdensity%s.png'%(self.figs_path,self.name_append))
        plt.savefig('%sdensity%s.png'%(self.figs_path,self.name_append))         
        
        
        ###############################
        # 3rd density fig
        ###############################
        
        nsub = 100
        depth  = copy.deepcopy(self.depth)
        density = cvalue
        swe = copy.deepcopy(self.state)
        
        depths = np.reshape(depth,(1,len(depth[:,0])*len(depth[0,:])))
        densitys = np.reshape(density,(1,len(depth[:,0])*len(density[0,:])))
        swesats = np.reshape(swe,(1,len(swe[:,0])*len(swe[0,:])))
        densitys = densitys[0,:]
        depths = depths[0,:]
        swesats = swesats[0,:]
        depthsub = depths[::nsub]
        densitysub = densitys[::nsub]
        swesub = swesats[::nsub]
        
        z = self.dem
        zs = np.reshape(z,(1,len(z[:,0])*len(z[0,:])))
        zs = zs[0,:]
        zsub = zs[::100]
        
        cm = cc.m_bgy
        
        plt.close(20)
        fig = plt.figure(num=20, figsize = (6,4), dpi=self.dpi)      
        ax = plt.gca()       
        h = ax.scatter(swesub,densitysub,c=zsub,vmin=5000,vmax=12500 ,cmap=cm,s=5)
        
        ax.set_ylabel(r'$\rho$ [kg/$m^3$]')
        ax.set_xlabel('SWE [in]')
        
        divider = make_axes_locatable(ax)
        cax = divider.append_axes("right", size="4%", pad=0.2)
        cbar = plt.colorbar(h, cax = cax)
        cbar.ax.tick_params() 
        cbar.set_label('elevation [ft]')
        plt.tight_layout()
              
        print('saving figure to %sdensity_swe%s.png'%(self.figs_path,self.name_append))
        plt.savefig('%sdensity_swe%s.png'%(self.figs_path,self.name_append))          
    
    
    def image_change(self,*args): 
        '''
        This plots self.delta_state
        
        Should add in more functionality for plotting differences between various images/runs
        
        '''  
        
        # if len(args) != 0:
        #     self.name_append = self.name_append + args[0]
              
        # Make copy so that we can add nans for the plots, but not mess up the original
        delta_state = copy.deepcopy(self.delta_state)
        qMin,qMax = np.percentile(delta_state,[1,99.5])
        
        ix = np.logical_and(delta_state < qMin, delta_state >= np.nanmin(np.nanmin(delta_state)))
        delta_state[ix] = qMin + qMin*0.2
        vMin,vMax = np.percentile(delta_state,[1,99])
        # clims = (qMin,qMax )
           
        # Override if absolute limits are provide in the config
        if hasattr(self,'ch_clminabs') and hasattr(self,'ch_clmaxabs'):
            clims       = (self.ch_clminabs,self.ch_clmaxabs)
     
        # if qMin is 0, need to change things 
        if qMax > 0:
            colorsbad = plt.cm.Accent_r(np.linspace(0., 1, 1))           
            colors1 = cmocean.cm.matter_r(np.linspace(0., 1, 127))
            colors2 = plt.cm.Blues(np.linspace(0, 1, 128))
            colors = np.vstack((colorsbad,colors1, colors2))
            mymap = mcolors.LinearSegmentedColormap.from_list('my_colormap', colors)

        else:
            colors1 = cmocean.cm.matter_r(np.linspace(0., 1, 255))
            colors2 = plt.cm.Set2_r(np.linspace(0, 1, 1))
            colors = np.vstack((colors1, colors2))
            mymap = mcolors.LinearSegmentedColormap.from_list('my_colormap', colors)
    
        ixf = delta_state == 0
        delta_state[ixf] = -100000 # set snow-free  
        pmask = self.masks[self.total_lbl]['mask']
        ixo = pmask == 0
        delta_state[ixo] = np.nan
        cmap = copy.copy(mymap)
        cmap.set_bad('white',1.)   
             
        sns.set_style('darkgrid')
        sns.set_context("notebook")
       
        plt.close(6) 
        fig,(ax,ax1) = plt.subplots(num=6, figsize=self.figsize, 
                                    dpi=self.dpi, nrows = 1, ncols = 2)
        h = ax.imshow(delta_state, interpolation='none', 
            cmap = cmap, norm=MidpointNormalize(midpoint=0,
                                                vmin = vMin-0.01,vmax=vMax+0.01))

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
        cax = divider.append_axes("right", size="5%", pad=0.2)
        cbar = plt.colorbar(h, cax = cax)     
        # cbar.ax.tick_params() 
        
        if self.units == 'KAF': 
            cbar.set_label(r'$\Delta$ SWE [in]')           
        if self.units == 'SI':
            cbar.set_label(r'$\Delta$ SWE [mm]')   
            
        h.axes.set_title('Change in SWE \n %s to %s'
                         %(self.dateFrom.date().strftime("%Y-%-m-%-d"),
                           self.dateTo.date().strftime("%Y-%-m-%-d")))    
      
        # Plot the bar in order
        sumorder  = self.plotorder[1:]  
        if self.basin == 'LAKES' or self.basin == 'RCEW':
            sumorder = [self.plotorder[0]]  
        if self.dplcs == 0:
            tlbl = '%s = %s %s'%(self.plotorder[0],
                                 str(int(self.delta_state_byelev[self.plotorder[0]].sum())),
                                 self.vollbl)
        else: 
            tlbl = '%s = %s %s'%(self.plotorder[0],
                                 str(np.round(self.delta_state_byelev[self.plotorder[0]].sum(),
                                              self.dplcs)),self.vollbl)            
        
        for iters,name in enumerate(sumorder):  
                
            if self.dplcs == 0:
                lbl = '%s = %s %s'%(name,
                                    str(int(self.delta_state_byelev[name].sum())),
                                    self.vollbl)
            else: 
                lbl = '%s = %s %s'%(name,
                                    str(np.round(self.delta_state_byelev[name].sum(),
                                    self.dplcs)),self.vollbl)                                 
 
            if iters == 0:
                ax1.bar(range(0,len(self.edges)),self.delta_state_byelev[name],
                        color = self.barcolors[iters],
                        edgecolor = 'k',label = lbl)
            elif iters == 1:   
                ax1.bar(range(0,len(self.edges)),self.delta_state_byelev[name], 
                        bottom = self.delta_state_byelev[sumorder[iters-1]], 
                        color = self.barcolors[iters], edgecolor = 'k',label = lbl)
            elif iters == 2:   
                ax1.bar(range(0,len(self.edges)),self.delta_state_byelev[name], 
                        bottom = self.delta_state_byelev[sumorder[iters-1]]
                        + self.delta_state_byelev[sumorder[iters-2]], 
                        color = self.barcolors[iters], edgecolor = 'k',label = lbl)
            elif iters == 3:   
                ax1.bar(range(0,len(self.edges)),self.delta_state_byelev[name], 
                        bottom = self.delta_state_byelev[sumorder[iters-1]]
                        + self.delta_state_byelev[sumorder[iters-2]]
                        + self.delta_state_byelev[sumorder[iters-3]], 
                        color = self.barcolors[iters], edgecolor = 'k',label = lbl)
               
        ax1.set_xlim((self.xlims[0]-0.5,self.xlims[1]))
        plt.tight_layout()                            
        xts = ax1.get_xticks()
        edges_lbl = []
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
                ax1.set_ylim((ylims[0]+(ylims[0]*0.3),(ylims[1] + ylims[1]*0.9)))  
                if (ylims[1] + ylims[1]*0.9) < abs(ylims[0]):
                    ax1.set_ylim((ylims[0]+(ylims[0]*0.3),(-(ylims[0]*0.6))))
                        
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
            ax.legend(handles=patches, bbox_to_anchor=(0.3, 0.05),
                      loc=2, borderaxespad=0. )
        elif self.basin == 'RCEW':
            ax.legend(handles=patches, bbox_to_anchor=(-0.1, 0.05),
                      loc=2, borderaxespad=0. )            
        else:
            ax.legend(handles=patches, bbox_to_anchor=(0.05, 0.05),
                      loc=2, borderaxespad=0. )
            
        if self.basin != 'LAKES' and self.basin != 'RCEW':
            # more ifs for number subs...
            if len(self.plotorder) == 5:
                ax1.legend(loc= (0.01,0.68))
            elif len(self.plotorder) == 4:
                ax1.legend(loc= (0.01,0.76))
            
        if self.basin == 'BRB':
            ax1.text(0.26,0.96,tlbl,horizontalalignment='center',
                     transform=ax1.transAxes,fontsize = 10)
        
        if self.basin == 'TUOL' or self.basin == 'SJ':
            ax1.text(0.3,0.94,tlbl,horizontalalignment='center',
                     transform=ax1.transAxes,fontsize = 10)
      
        plt.tight_layout() 
        fig.subplots_adjust(top=0.88)
        
        print('saving figure to %sswe_change%s.png'%(self.figs_path,self.name_append))
        plt.savefig('%sswe_change%s.png'%(self.figs_path,self.name_append))    
        
    def state_by_elev(self): 
        '''
        Plots SWE by elevation, delineated by melt/nonmelt
        
        '''
            
        lim = np.max(self.melt[self.total_lbl]) + np.max(self.nonmelt[self.total_lbl])
        ylim = np.max(lim) + np.max(lim)*0.3 
        colors = ['xkcd:rose red','xkcd:cool blue']
        fs = list(self.figsize)
        fs[0] = fs[0]*0.9
        fs = tuple(fs)

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
            plt.close(7)
            fig,ax  = plt.subplots(num=7, figsize=fs, dpi=self.dpi,
                                   nrows = nr, ncols = nc)
            axs = ax.ravel()
            if nf == 5:
                fig.delaxes(axs[5])
                   
            for iters,name in enumerate(self.plotorder):
                # iters = 0
                # name = self.plotorder[iters]
                axs[iters].bar(range(0,len(self.edges)),self.melt[name], 
                               color = colors[0], bottom = self.nonmelt[name])
                axs[iters].bar(range(0,len(self.edges)),self.nonmelt[name], 
                               color = colors[1], label = 'unavail ')
                
                axs[iters].set_xlim(self.xlims)                          
                xts = axs[iters].get_xticks()
                if len(xts) < 6:
                    dxt = xts[1] - xts[0]
                    xts = np.arange(xts[0],xts[-1] + 1 ,dxt/2)
                    axs[iters].set_xticks(xts)                 
        
                edges_lbl = []
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
                kaf = str(np.int(sum(self.melt[name])
                                 + sum(self.nonmelt[name]))) 
                
                if self.dplcs == 0:
                    kaf = str(np.int(sum(self.melt[name])
                                     + sum(self.nonmelt[name]))) 
                else: 
                    kaf = str(np.round(sum(self.melt[name])
                                       + sum(self.nonmelt[name]),self.dplcs))             
    
                if self.units == 'KAF':          
                    axs[iters].text(0.5,0.92,'%s - %s KAF'
                                    %(self.masks[name]['label'],kaf),
                                    horizontalalignment='center',
                                    transform=axs[iters].transAxes, fontsize = 10)
                if self.units == 'SI':          
                    axs[iters].text(0.5,0.92,'%s - %s $km^3$'
                                    %(self.masks[name]['label'],kaf),
                                    horizontalalignment='center',
                                    transform=axs[iters].transAxes)                
                
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

                        tmpa = (r'avail = %s')%(kafa)                         
                        lbl.append(tmpa)

                    if self.dplcs == 0:
                        kafna = str(np.int(sum(self.nonmelt[name]))) 
                    else: 
                        kafna = str(np.round(sum(self.nonmelt[name]),self.dplcs))            

                    tmpna = ('unavail = %s')%(kafna)
                    lbl.append(tmpna)                     

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
                   
            ax.bar(range(0,len(self.edges)),self.melt[name],
                   color = colors[0], bottom = self.nonmelt[name])
            ax.bar(range(0,len(self.edges)),self.nonmelt[name],
                   color = colors[1], label = 'unavail ')
            
            ax.set_xlim(self.xlims)                          
            xts = ax.get_xticks()
            edges_lbl = []
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
            kaf = str(np.int(sum(self.melt[name]) + sum(self.nonmelt[name]))) 
            
            if self.dplcs == 0:
                kaf = str(np.int(sum(self.melt[name]) 
                                 + sum(self.nonmelt[name]))) 
            else: 
                kaf = str(np.round(sum(self.melt[name]) 
                                   + sum(self.nonmelt[name]),self.dplcs))             

            if self.units == 'KAF':          
                ax.text(.25,0.95,'%s - %s KAF'
                        %(self.masks[name]['label'],kaf),
                        horizontalalignment='center',transform=ax.transAxes)
                ax.set_ylabel('KAF')
            if self.units == 'SI':          
                ax.text(-0.88,2.1,'%s - %s $km^3$'
                        %(self.masks[name]['label'],kaf),
                        horizontalalignment='center',transform=ax.transAxes)
                ax.set_ylabel(r'$km^3$')                
            
            lbl = []
            for n in (0,1):
                if n == 0:
                    if self.dplcs == 0:
                        kafa = str(np.int(sum(self.melt[name]))) 
                    else: 
                        kafa = str(np.round(sum(self.melt[name]),self.dplcs)) 
                    
                    tmpa = ('avail = %s')%(kafa)
                                          
                    lbl.append(tmpa)
                
                # kafna = str(np.int(sum(self.nonmelt[name]))) 
                if self.dplcs == 0:
                    kafna = str(np.int(sum(self.nonmelt[name]))) 
                else: 
                    kafna = str(np.round(sum(self.nonmelt[name]),self.dplcs))            

                tmpna = ('unavail = %s')%(kafna)                     
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
        fig.suptitle('SWE, %s'%self.dateTo.date().strftime("%Y-%-m-%-d"))
              
        print('saving figure to %sswe_elev%s.png'%(self.figs_path,self.name_append))   
        plt.savefig('%sswe_elev%s.png'%(self.figs_path,self.name_append))  
        
    def basin_total(self):       
             
        sns.set_style('darkgrid')
        sns.set_context("notebook")
        
        plt.close(8)
        fig,(ax,ax1) = plt.subplots(num=8, figsize=self.figsize, 
                                    dpi=self.dpi, nrows = 1, ncols = 2)
      
        self.barcolors.insert(0,'black')
        
        if self.basin == 'LAKES' or self.basin == 'RCEW':
            plotorder = [self.plotorder[0]]
        else:
            plotorder = self.plotorder
        
        for iters,name in enumerate(plotorder):
            # name = self.plotorder[0]
            # iters = 0
            self.state_summary[name].plot(ax=ax, color = self.barcolors[iters])             
            ax1.plot(self.accum_summary[name], 
                     color = self.barcolors[iters], label='_nolegend_')
        
        ax1.yaxis.set_label_position("right")
        ax1.set_xlim((datetime(self.wy -1 , 10, 1),self.dateTo))
        ax.set_xlim((datetime(self.wy - 1, 10, 1),self.dateTo))
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
             
        ax1.set_ylabel(r'[%s]'%(self.vollbl))
        ax.axes.set_title('Basin SWE')
        ax1.axes.set_title('Accumulated Basin SWI')
        ax.set_ylabel(r'[%s]'%(self.vollbl)) 
        
        # plt.tight_layout()   
        del self.barcolors[0] 
        
        print('saving figure to %sbasin_total%s.png'%(self.figs_path,
                                                      self.name_append))   
        plt.savefig('%sbasin_total%s.png'%(self.figs_path,self.name_append))       

        ########################################
        #         Second figure
        ########################################
        
        # This needs to be improved...
        
        if self.basin == 'BRB':
            accum_summary = self.accum_summary
            main = 'Boise River Basin'
            multiswe = pd.DataFrame.from_csv(
                '/mnt/volumes/wkspace/results/brb/brb_multiyear_summary.csv') 
            multiswi = pd.DataFrame.from_csv(
                '/mnt/volumes/wkspace/results/brb/brb_multiyear_swi.csv')
            
            multiswe.wy17.iloc[304:] = 0
            multiswi.wy17 = np.cumsum(multiswi.wy17)
            multiswi.wy17.iloc[304:] = multiswi.wy17[303]
            
            state_summary = self.state_summary.asfreq('D')
            
            # Put in this year
            multiswe.wy18.iloc[:len(state_summary[main])] = (
                                                state_summary[main].values) 
            multiswi.wy18.iloc[:len(accum_summary[main])] = (
                                                accum_summary[main].values)          
            
            if self.units == 'SI':
                multiswe.wy17 = np.multiply(multiswe.wy17,0.00123348)
                multiswi.wy17 = np.multiply(multiswi.wy17,0.00123348)
                multiswe.wy16 = np.multiply(multiswe.wy16,0.00123348)
                multiswi.wy16 = np.multiply(multiswi.wy16,0.00123348)
                multiswe.wy15 = np.multiply(multiswe.wy15,0.00123348)
                multiswi.wy15 = np.multiply(multiswi.wy15,0.00123348)                
            
            plt.close(8)
            fig,(ax,ax1)    = plt.subplots(num=8, figsize=self.figsize,
                                           dpi=self.dpi, nrows = 1, ncols = 2)
    
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
                 
            if self.units == 'KAF':
                ax.set_ylabel('[KAF]')  
                ax1.set_ylabel('[KAF]')
                ax.axes.set_title('Water Year SWE')
                ax1.axes.set_title('Accumulated Basin SWI')
            
            if self.units == 'SI':
                ax.set_ylabel(r'[$km^3$]') 
                ax1.set_ylabel(r'[$km^3$]')            
                ax.axes.set_title('Basin SWE [$km^3$]')
                ax1.axes.set_title('Accumulated Basin SWI [$km^3$]')
            
            ax.set_ylim((0,ax1.get_ylim()[1]))
            plt.tight_layout()      
            
            print('saving figure to %sbasin_total_multiyr%s.png'%(self.figs_path,self.name_append))   
            plt.savefig('%sbasin_total_multiyr%s.png'%(self.figs_path,self.name_append))                   
             
    def stn_validate(self,initswe):
        
        rundirs = self.nc_dir
        stns = self.val_stns
        lbls = self.val_lbls
        client = self.val_client
           
        # get metadata from the data base from snotel sites
        if self.basin == 'BRB':
            qry = ('SELECT tbl_metadata.* FROM tbl_metadata '
                   + 'INNER JOIN tbl_stations ON tbl_metadata.primary_id = '
                   + 'tbl_stations.station_id WHERE tbl_stations.client = '
                   + ' "'"%s"'" HAVING network_name = "'"SNOTEL"'";'%client)
        else:
            qry = ('SELECT tbl_metadata.* FROM tbl_metadata ' 
                   + 'INNER JOIN tbl_stations ON tbl_metadata.primary_id ='
                   + 'tbl_stations.station_id WHERE tbl_stations.client = '
                   + '"'"%s"'" ;'%client)        
        cnx = mysql.connector.connect(user='markrobertson', 
                                      password='whatdystm?1',
                                      host='10.200.28.137',
                                      database='weather_db')
        
        meta_sno = pd.read_sql(qry, cnx)
        meta_sno.index = meta_sno['primary_id']
        swe_meas    = pd.DataFrame(index = pd.date_range(datetime(2017,10,1), 
                                                         self.dateTo, 
                                                         freq='D'),columns = stns)  
        swe_mod     = pd.DataFrame(index = pd.date_range(datetime(2017,10,1),
                                                         self.dateTo, 
                                                         freq='D'),columns = stns)   
        tbl         = 'tbl_level1'
        var         = 'snow_water_equiv'
        st_time     = '2017-10-1 00:00:00'
        end_time    = self.dateTo.date().strftime("%Y-%-m-%-d")
        
        # Get Snotel station results
        for iters,stn in enumerate(stns): 
            cnx     = mysql.connector.connect(user='markrobertson',
                                              password='whatdystm?1',
                                              host='10.200.28.137',
                                              port='32768',database='weather_db')
            var_qry = ('SELECT weather_db.%s.date_time, weather_db.%s.%s ' % (tbl,tbl,var) +
                        'FROM weather_db.%s ' % tbl +
                        "WHERE weather_db.%s.date_time between '" % tbl + st_time+ "' and '"+end_time+"'"
                        "AND weather_db.%s.station_id IN ('" % tbl + stn + "');")
            
            data = pd.read_sql(var_qry, cnx, index_col='date_time')
            dind = pd.date_range(st_time,end_time,freq='D')
            swe_meas[stn] = data.reindex(dind)
            
        if self.basin == 'TUOL':
            swe_meas.TIOC1 = swe_meas.TIOC1 - 300  
            
        if 'AGP' in swe_meas:
            swe_meas.AGP = swe_meas.AGP - 40   
        if 'VLC' in swe_meas:
            swe_meas.VLC = swe_meas.VLC + 250  
        if 'UBC' in swe_meas:
            swe_meas.UBC = swe_meas.UBC - 50      
        
        sns.set_style('darkgrid')
        sns.set_context('notebook')
        
        plt.close(9)
        fig, axs = plt.subplots(num = 9,figsize = (10,10),nrows = 3,ncols = 2)   
        axs = axs.flatten() 
        
        ### sdwitching the loop!###
        
        # First need to combine all nc files... 
        px = (1,1,1,0,0,0,-1,-1,-1)
        py = (1,0,-1,1,0,-1,1,0,-1)       
        for iters,stn in enumerate(stns):
            # iters = 0
            # stn = stns[iters]
            for n,m in zip(px,py): 
                # n = 0
                # m = 0
                # iswe = 0  
                # iswe = 31+19 # brb
                iswe = initswe
        
                for rname in rundirs:
                    # rname = rundirs[1]
                    
                    try:
                        ncpath  = rname.split('output')[0]
                        ncf     = nc.Dataset(ncpath + 'snow.nc', 'r')    # open netcdf file
                    except:
                        ncpath  = rname
                        ncf     = nc.Dataset(ncpath + 'snow.nc', 'r')    # open netcdf file
                        
                    nctvec  = ncf.variables['time'][:]
                    vswe    = ncf.variables['specific_mass']            # get variable
                    ncxvec  = ncf.variables['x'][:]                     # get x vec
                    ncyvec  = ncf.variables['y'][:]                     # get y vec                      
                    ll      = utm.from_latlon(meta_sno.ix[stn,'latitude'],meta_sno.ix[stn,'longitude']) # get utm coords from metadata
                    # ll      = utm.from_latlon(37.641922,-119.055443)
                    # ll      = utm.from_latlon(37.655201,-119.060783)
        
                    xind    = np.where(abs(ncxvec-ll[0]) == min(abs(ncxvec-ll[0])))[0]  # get closest pixel index to the station
                    yind    = np.where(abs(ncyvec-ll[1]) == min(abs(ncyvec-ll[1])))[0]  # get closest pixel index to the station
                    # print(xind,yind)
                    swe     = pd.Series(vswe[:,yind+m,xind+n].flatten(),index=nctvec)  # pull out closest model pixel data
             
                    try:
                        swe_mod.loc[iswe:(iswe + len(swe.values)),stn] = swe.values  
                    except:
                        sv = swe_mod[stn].values
                        lx = len(sv[iswe::])
                        swe_mod.loc[iswe:(iswe + lx),stn] = swe.values[0:lx]
                        
                    ncf.close()   
                    iswe = iswe + len(swe.values)
                    
                z = self.dem[yind,xind]    
               
                axs[iters].plot(swe_meas[stn],'k',label='measured')
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
             
        # Plot
        maxm = np.nanmax(swe_meas)
        maxi = np.nanmax(swe_mod.max().values)
        
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
        
        self.swe_meas = swe_meas
        self.swe_mod = swe_mod
        
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
       
        plt.close(10) 
        fig,(ax,ax1)    = plt.subplots(num=10, figsize=self.figsize, dpi=self.dpi, nrows = 1, ncols = 2)
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

    
    def distribution_detail(self,*args): 
        '''
        Maybe get rid of this one, not using at this point...
        ''' 
                
        colors  = ['xkcd:rose red','xkcd:cool blue']
        
        if len(args) != 0:
            if args[0] == 'depth':
                name = self.plotorder[args[1]]
                print(name)
                
                value = self.state_mswe_byelev
                lim     = np.max(value[self.total_lbl]) 
                ylim    = np.max(lim) + np.max(lim)*0.15 

            elif args[0] == 'volume':
                value   =  self.melt
                value2  =  self.nonmelt 
                lim     = np.max(value[self.total_lbl]) + np.max(value[self.total_lbl])
                ylim    = np.max(lim) + np.max(lim)*0.1                        
        
        sns.set_style('darkgrid')
        
        plt.close(5)
        fig  = plt.figure(num=5, figsize=(self.figsize[0]*0.9,self.figsize[1]), dpi=self.dpi)
        ax = plt.gca()
     
        if 'value2' in locals():
            ax.bar(range(0,len(self.edges)),value[name], color = colors[0], bottom = self.nonmelt[name])
            ax.bar(range(0,len(self.edges)),value2[name], color = colors[1], label = 'unavail ') 
            
            # Get basin total storage in strings for label
            kaf    = str(np.int(sum(self.melt[name]) + sum(self.nonmelt[name]))) 
                
            if self.dplcs == 0:
                kaf    = str(np.int(sum(self.melt[name]) + sum(self.nonmelt[name]))) 
            else: 
                kaf    = str(np.round(sum(self.melt[name]) + sum(self.nonmelt[name]),self.dplcs))             
    
            ax.text(.25,0.95,'%s - %s %s'%(self.masks[name]['label'],kaf,self.vollbl),horizontalalignment='center',transform=ax.transAxes)
            ax.set_ylabel('%s'%(self.vollbl))
                  
                
            lbl = []
            for n in (0,1):
                if n == 0:
                    if self.dplcs == 0:
                        kafa = str(np.int(sum(self.melt[name]))) 
                    else: 
                        kafa = str(np.round(sum(self.melt[name]),self.dplcs)) 
                    
                    tmpa = ('avail = %s')%(kafa)
                                          
                    lbl.append(tmpa)
                    
                # kafna = str(np.int(sum(self.nonmelt[name]))) 
                if self.dplcs == 0:
                    kafna = str(np.int(sum(self.nonmelt[name]))) 
                else: 
                    kafna = str(np.round(sum(self.nonmelt[name]),self.dplcs))            
    
                tmpna = ('unavail = %s')%(kafna)                     
                lbl.append(tmpna) 
                
            ax.legend(lbl, loc = (0.025, 0.8),fontsize = 9)    
            fig.suptitle('SWE, %s'%self.dateTo.date().strftime("%Y-%-m-%-d"))                                                   
        
        else:
            ax.bar(range(0,len(self.edges)),value[name], color = 'xkcd:peacock blue')  
               
            ax.set_ylabel('Mean SWE [%s]'%(self.depthlbl)) 
            fig.suptitle('%s, Mean SWE, %s'%(name,self.dateTo.date().strftime("%Y-%-m-%-d")) ) 
      
        
        ax.set_xlim(self.xlims)                          
        xts         = ax.get_xticks()
        edges_lbl   = []
        for i in xts[0:len(xts)-1]:
            edges_lbl.append(str(int(self.edges[int(i)])))

        ax.set_xticklabels(str(i) for i in edges_lbl)
        ax.set_xlabel('elevation [%s]'%(self.elevlbl))
        
        ax.set_ylim((0,ylim)) 
        for tick in ax.get_xticklabels():
            tick.set_rotation(30)             
        
        fig.tight_layout()
        ax.set_xticks(xts)
        ax.set_xlim(self.xlims)
  
        fig.subplots_adjust(top=0.92,wspace = 0.1)  

        print('saving figure to %smean_detail%s.png'%(self.figs_path,self.name_append))
        plt.savefig('%smean_detail%s.png'%(self.figs_path,self.name_append))  
    
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
        