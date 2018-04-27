
import numpy as np
from matplotlib import pyplot as plt
from mpl_toolkits.axes_grid1 import make_axes_locatable
import seaborn as sns
from smrf import ipw
from shutil import copyfile
import os
import copy
import pandas as pd
import ConfigParser as cfp
import snowav.methods.wyhr_to_datetime as wy


class SNOWAV(object):

    def __init__(self,config_file):
        '''
        Notes: 

        '''
        
        try:
            if not os.path.isfile(config_file):
                print('Config file does not exist!')
                self.error = True
                return
                
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
            
            ####################################################
            #           Runs                                   #
            ####################################################       
            # Collect all the run directories
            self.run_dirs = list(cfg.items('Runs'))             
            self.snow_files = []
            self.em_files = []
            
            for rdir in self.run_dirs:

                run_files = [rdir[1] + s for s in sorted(os.listdir(rdir[1]))]
                
                self.snow_files = (self.snow_files 
                                   + [value for value in run_files 
                                   if ( ('snow.' in value) and not ('.nc' in value))])
                self.em_files = (self.em_files 
                                 + [value for value in run_files 
                                 if ( ('em.' in value) and not ('.nc' in value))])
            
            while '*snow.nc' in self.snow_files:
                self.snow_files.remove('*snow.nc')
            
            while '*em.nc' in self.em_files:
                self.em_files.remove('*em.nc')                
       
            # If no psnowFile and csnowFile specified, use first and last
            if not cfg.has_option('Outputs','csnowFile'):
                print('Using first and last outputs as start and end times.')
                self.psnowFile = self.snow_files[0] 
                self.csnowFile = self.snow_files[len(self.snow_files)-1] 
                self.cemFile = self.em_files[len(self.em_files)-1]
                
            # Check to see if they exist
            if not (os.path.isfile(self.psnowFile)):   
                print('psnowFile %s does not exist!'%(self.psnowFile))
                self.error = True
                return            
            
            if not (os.path.isfile(self.csnowFile)):   
                print('csnowFile does not exist!')
                self.error = True
                return
            
            ####################################################
            #           Validate                               #
            #################################################### 
            if (cfg.has_option('Validate','stations')
                and (cfg.has_option('Validate','labels'))
                and (cfg.has_option('Validate','client'))
                ):
                self.val_stns = cfg.get('Validate','stations').split(',')
                self.val_lbls = cfg.get('Validate','labels').split(',')
                self.val_client = cfg.get('Validate','client')   
                self.valid_flag = True
            else:
                self.valid_flag = False
                print('No validation stations listed, will not generate figure')
            
            # This is being used to combine 2017 HRRR data
            if cfg.has_option('Validate','offset'):
                self.offset = int(cfg.get('Validate','offset'))   
            else:
                self.offset = 0                    
            ####################################################
            #           Accumulated                            #
            ####################################################   
            # Right now we aren't using these because matplotlib and multiple
            # colormaps don't get along well with clims...
            # self.acc_clmin = cfg.get('Accumulated','clmin')
            # self.acc_clmax = cfg.get('Accumulated','clmax') 
            if (cfg.has_option('Accumulated','ymin') 
                and cfg.has_option('Accumulated','ymax')):
    
                self.acc_ylims = (int(cfg.get('Accumulated','ymin')),
                                  int(cfg.get('Accumulated','ymax')))  
            
            if cfg.has_option('Accumulated','save_fig'):
                self.acc_flag = cfg.get('Accumulated','save_fig')
            else:
                self.acc_flag = True

            if cfg.has_option('Accumulated','min_swi'):
                self.min_swi = cfg.get('Accumulated','min_swi')   
                       
    
            ####################################################
            #           Elevation                              #
            ####################################################             
            if (cfg.has_option('Elevation','ymin') 
                and cfg.has_option('Elevation','ymax')):
                self.elv_ylims = (int(cfg.get('Elevation','ymin')),
                                 int(cfg.get('Elevation','ymax')))  
            
            if cfg.has_option('Elevation','save_fig'):
                self.elv_flag = cfg.get('Elevation','save_fig')
            else:
                self.elv_flag = True                                   
            
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
            #           Basin Total                            #
            ####################################################          
            if (cfg.has_option('Basin Total','summary_swe') and 
                cfg.has_option('Basin Total','summary_swi')):
                self.summary_swe = cfg.get('Basin Total','summary_swe')
                self.summary_swi = cfg.get('Basin Total','summary_swi')  
            if not (os.path.isfile(self.summary_swe) and 
                    os.path.isfile(self.summary_swi) ):   
                print('Failed reading in Basin Total section!')
                self.error = True
                return 
            
            #
            if cfg.has_option('Basin Total','netcdf'):
                self.ncvars = cfg.get('Basin Total','netcdf').split(',')     
                self.nc_flag = True
            else:
                self.nc_flag = False    
            
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
            self.report_flag = cfg.get('Report','report')  
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

            if not (os.path.isfile(self.templ_path + self.tex_file)):
                print('Error reading in Reports section!')
                self.error = True
                return                  
              
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
            
            try:
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
                # HACK FOR SJ SUB4!
                self.masks = dict()
                for lbl,mask in zip(self.plotorder,maskpaths):
                    if (self.basin == 'SJ' and lbl == self.sub4_lbl):
                        self.masks[lbl] = {'border': blank, 
                                           'mask': np.genfromtxt(mask,skip_header=0),
                                           'label': lbl}   
                    # hack
                    else:
                        self.masks[lbl] = {'border': blank, 
                                           'mask': np.genfromtxt(mask,skip_header=sr),
                                           'label': lbl}                              
            except:
                print('Error creating mask dicts!')
                self.error = True
                return
                  
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
        
        except:
            print('Error reading config file.')

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
        evap = copy.deepcopy(accum)
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
        evap_byelev = copy.deepcopy(accum_byelev)
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
        evap_summary = pd.DataFrame(columns = self.masks.keys())
               
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
            evap = evap + em_file.bands[6].data
 
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
                if hr < 100:
                    shr = '00' + str(hr)
                    
                if hr >=100 and hr < 1000:
                    shr = '0' + str(hr)  
                    
                if hr >= 1000:
                    shr = str(hr)     
                
                if os.path.isfile(ppt_path +'ppt.4b_'+ shr):
                    pFlag = True
                    ppt_files = ppt_files + [ppt_path +'ppt.4b_'+ shr]
            
            rain_hrly = np.zeros((self.nrows,self.ncols))
            precip_hrly = np.zeros((self.nrows,self.ncols))
            
            # Load 'em in
            if pFlag and accum_sub_flag == True:
                print('adding precip for %s'%(snow_name))
                for pfile in ppt_files:
                    ppt = ipw.IPW(pfile)
                    pre = ppt.bands[0].data
                    percent_snow = ppt.bands[1].data
                    rain_hrly = rain_hrly + np.multiply(pre,(1-percent_snow))
                    precip_hrly = precip_hrly + pre
            else:
                print(pFlag,accum_sub_flag)
            
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
                evap_summary.loc[date, name] = (np.nansum(
                                                np.multiply(evap,
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
            evap_mask = np.multiply(evap,mask)
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
                evap_byelev.loc[b,name] = np.nansum(evap_mask[ind])
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
        self.evap_byelev = np.multiply(evap_byelev,self.conversion_factor) 
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
        self.evap_summary = np.multiply(evap_summary,self.conversion_factor)                                
        self.state_byday = np.multiply(state_byday,self.depth_factor)  
    
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
        