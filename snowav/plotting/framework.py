
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
import datetime


class SNOWAV(object):

    def __init__(self,config_file = None):
        '''
        Initializing reporting object, and read in configuration file.
        '''

        try:
            if not os.path.isfile(config_file):
                print('SNOWAV config file does not exist!')
                self.error = True
                return

            print('Reading SNOWAV config file...')
            cfg = cfp.ConfigParser()
            cfg.read(config_file)

            ####################################################
            #             Basin section                        #
            ####################################################
            self.basin = cfg.get('Basin','basin')
            if (cfg.has_option('Basin','save_path') and
                os.path.exists(cfg.get('Basin','save_path'))):
                self.save_path = cfg.get('Basin','save_path')
            else:
                print('Figures save_path either not specified in config file'
                      ' or does not exist, using ./snowav/data/'
                      + '\nFigures will be saved but report will fail')
                self.save_path = './snowav/data/'
                
            if cfg.has_option('Basin','name_append'):
                self.name_append = cfg.get('Basin','name_append')
            else:
                self.name_append = '_gen_' + datetime.datetime.now().strftime("%Y-%-m-%-d")
            
            self.wy = int(cfg.get('Basin','wy'))
            
            if cfg.has_option('Basin','units'):
                self.units = cfg.get('Basin','units')
            else:
                self.units = 'KAF'
                
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
                
                # Check to see if they exist
                if not (os.path.isfile(cfg.get('Outputs','psnowFile'))):
                    print('psnowFile does not exist!')
                    self.error = True
                    return

                if not (os.path.isfile(cfg.get('Outputs','csnowFile'))):
                    print('csnowFile does not exist!')
                    self.error = True
                    return                  
                
                self.psnowFile = cfg.get('Outputs','psnowFile')
                self.csnowFile = cfg.get('Outputs','csnowFile')
                self.cemFile = self.csnowFile.replace('snow.','em.')    
                
            # Flights
            if (cfg.has_option('Outputs','fltpsnowFile') and
                cfg.has_option('Outputs','fltcsnowFile')): 
                                                        # Check to see if they exist
                if not (os.path.isfile(cfg.get('Outputs','fltpsnowFile'))):
                    print('fltpsnowFile does not exist!')
                    self.error = True
                    return

                if not (os.path.isfile(cfg.get('Outputs','fltcsnowFile'))):
                    print('fltcsnowFile does not exist!')
                    self.error = True
                    return 
                
                self.fltpsnowFile = cfg.get('Outputs','fltpsnowFile')
                self.fltcsnowFile = cfg.get('Outputs','fltcsnowFile')
                self.flt_flag = True                   

            ####################################################
            #           Runs                                   #
            ####################################################
            # Collect all the run directories
            self.run_dirs = list(cfg.items('Runs'))

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
            if ((cfg.has_option('Changes','clmin'))
                and (cfg.has_option('Changes','clmax'))):
                self.ch_clmin = cfg.get('Changes','clmin')
                self.ch_clmax = cfg.get('Changes','clmax')
            else:
                self.ch_clmin = 0.01
                self.ch_clmax = 99.9

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
            if ((cfg.has_option('Results','clmin'))
                and (cfg.has_option('Results','clmax'))):
                self.ch_clmin = cfg.get('Results','clmin')
                self.ch_clmax = cfg.get('Results','clmax')
            else:
                self.ch_clmin = 0.01
                self.ch_clmax = 99.9

            ####################################################
            #           Basin Total                            #
            ####################################################
            if cfg.has_option('Basin Total','summary_swe'):
                self.summary_swe = cfg.get('Basin Total','summary_swe')
                self.summary_swi = cfg.get('Basin Total','summary_swi')
                if not (os.path.isfile(self.summary_swe) and
                        os.path.isfile(self.summary_swi) ):
                    print('Failed reading in Basin Total section!')
                    self.error = True
                    return

            if cfg.has_option('Basin Total','netcdf'):
                self.ncvars = cfg.get('Basin Total','netcdf').split(',')
                self.nc_flag = True
            else:
                self.nc_flag = False

            ####################################################
            #           DEM                                    #
            ####################################################
            self.dempath = cfg.get('Masks','dempath')
            self.total = cfg.get('Masks','total')

            ####################################################
            #          Plots                                   #
            ####################################################
            if cfg.has_option('Plots','fig_length'):
                self.figsize = (int(cfg.get('Plots','fig_length')),
                                int(cfg.get('Plots','fig_height')))
                self.dpi = int(cfg.get('Plots','dpi'))
                self.barcolors = ['xkcd:true green','palegreen',
                              'xkcd:dusty green','xkcd:vibrant green','red']
            else:
                self.figsize = (10,5)
                self.dpi = 200
                self.barcolors = ['xkcd:true green','palegreen',
                              'xkcd:dusty green','xkcd:vibrant green','red']                

            ####################################################
            #          Report                                  #
            ####################################################
            try:
                # Report defaults to True
                if cfg.has_option('Report','report'):
                    self.report_flag = cfg.getboolean('Report','report')
                                
                if self.save_path == './snowav/data/': 
                    self.report_flag = False
                    print('Report flag being set to false becaues of relative'
                          + ' figure path ./snowav/data/')
                else:
                    self.report_flag = True
                    
                # Add date if necessary
                if cfg.has_option('Report','orig_date'):
                    self.orig_date = cfg.get('Report','orig_date')
                    
                if cfg.has_option('Report','exclude_figs'):
                    self.exclude_figs = cfg.get('Report','exclude_figs').split(',') 
                    options = ['CHANGES','SWI','RESULTS','ELEV','TOTALS',
                               'MEAN','VALID']
                    for name in self.exclude_figs:
                        if name not in options:
                            print('[Report] exclude_fig options are: %s'%(options))
                            self.error = True
                            return
                    
                # These will later get appended with self.dateTo
                if cfg.has_option('Report','report_name'):
                    self.report_name = cfg.get('Report','report_name')
                else:
                    self.report_name = 'SnowpackSummary.pdf'
                
                if cfg.has_option('Report','report_title'):
                    self.rep_title = cfg.get('Report','report_title')     
                else: 
                    self.rep_title = 'Snowpack Summary'                                                                    
                                      
                if cfg.has_option('Report','rep_path'):    
                    self.rep_path = cfg.get('Report','rep_path')
                else:
                    if os.path.exists('./snowav/data/'):
                        print('No save path for report given in config file,'
                              + ' using ./snowav/data/')
                        self.rep_path = './snowav/data/'
                    else:
                        print('No save path for report given in config file,'
                              + ' either list in [Report] rep_path or run'
                              + ' in /SNOWAV')
                        self.error = True
                        return
                
                if cfg.has_option('Report','env_path'):    
                    self.env_path = cfg.get('Report','env_path')
                else:
                    if os.path.exists('./snowav/report/template/section_text/'):
                        print('No environment path for report given in config '
                              + 'file, using ./snowav/report/template/section_text/')
                        self.env_path = './snowav/report/template/section_text/'
                    else:
                        print('No environment path for report given in config file,'
                              + ' either list in [Report] env_path or run'
                              + ' in /SNOWAV')
                        self.error = True
                        return  
                                  
                if cfg.has_option('Report','templ_path'):    
                    self.templ_path = cfg.get('Report','templ_path')
                else:
                    if os.path.exists('./snowav/report/template/'):
                        print('No environment path for report given in config '
                              + 'file, using ./snowav/report/template/')
                        self.templ_path = './snowav/report/template/'
                    else:
                        print('No environment path for report given in config file,'
                              + ' either list in [Report] templ_path or run'
                              + ' in /SNOWAV')
                        self.error = True
                        return  

                if cfg.has_option('Report','tex_file'):    
                    self.tex_file = cfg.get('Report','tex_file')
                else:
                    if os.path.isfile('./snowav/report/template/snowav_report.tex'):
                        print('No LaTeX file for report given in config '
                              + 'file, using ./snowav/report/template/snowav_report.tex')
                        self.tex_file = 'snowav_report.tex'
                    else:
                        print('No LaTeX file for report given in config file,'
                              + ' either list in [Report] tex_file or run'
                              + ' in /SNOWAV')
                        self.error = True
                        return  

                if cfg.has_option('Report','summary_file'):    
                    self.summary_file = cfg.get('Report','summary_file')
                else:
                    if os.path.isfile('./snowav/report/template/section_text/report_summary.txt'):
                        print('No summary file for report given in config '
                              + 'file, using ./snowav/report/template/section_text/report_summary.txt')
                        self.summary_file = './snowav/report/template/section_text/report_summary.txt'
                    else:
                        print('No LaTeX file for report given in config file,'
                              + ' either list in [Report] summary_file or run'
                              + ' in /SNOWAV')
                        self.error = True
                        return  
                
                if cfg.has_option('Report','figs_tpl_path'):    
                    self.figs_tpl_path = cfg.get('Report','figs_tpl_path')
                else:
                    if os.path.exists('./snowav/report/figs/'):
                        print('No figs template path for report given in config '
                              + 'file, using ./snowav/report/figs/')
                        self.figs_tpl_path = './snowav/report/figs/'
                    else:
                        print('No figs template path for report given in config file,'
                              + ' either list in [Report] figs_tpl_path or run'
                              + ' in /SNOWAV')
                        self.error = True
                        return                      
           
            except:
                print('Error reading in Reports section!')
                self.error = True
                return
            
            # Strings for the report
            if self.units == 'KAF':
                self.reportunits = 'KAF'
            if self.units == 'SI':
                self.reportunits = 'km^3'

            ####################################################
            #           History forecast
            ####################################################
            if cfg.has_option('Hx Forecast','adj_hours'):
                self.adj_hours = int(cfg.get('Hx Forecast','adj_hours'))

            ####################################################
            #           Masks
            ####################################################
            self.subbasin1 = cfg.get('Masks','subbasin1')
            self.subbasin2 = cfg.get('Masks','subbasin2')
            self.subbasin3 = cfg.get('Masks','subbasin3')
            self.total_lbl = cfg.get('Masks','total_lbl')
            self.sub1_lbl = cfg.get('Masks','sub1_lbl')
            self.sub2_lbl = cfg.get('Masks','sub2_lbl')
            self.sub3_lbl = cfg.get('Masks','sub3_lbl')

            self.plotorder = [self.total_lbl, self.sub1_lbl,
                                 self.sub2_lbl, self.sub3_lbl]

            self.suborder = [self.sub1_lbl,self.sub2_lbl,self.sub3_lbl]
            maskpaths = [self.total, self.subbasin1,
                            self.subbasin2,self.subbasin3 ]

            # Add if necessary - need to generalize all this and change
            # in the config file!
            if cfg.has_option('Masks','sub4_lbl'):
                self.sub4_lbl = cfg.get('Masks','sub4_lbl')
                self.subbasin4 = cfg.get('Masks','subbasin4')
                self.plotorder = self.plotorder + [self.sub4_lbl]
                self.suborder = self.suborder + [self.sub4_lbl]
                maskpaths = maskpaths + [self.subbasin4]

            # Collect the run directories
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
            if not hasattr(self,'csnowFile'):
                self.psnowFile = self.snow_files[0]
                self.csnowFile = self.snow_files[-1]
                self.cemFile = self.em_files[-1] 
                print('psnowFile and/or csnowFile not specified, using:' 
                      + ' \n%s and \n%s'%(self.psnowFile,self.csnowFile))             

            # Get the DEM
            # There are different formats, this will get fixed once we
            # start using netcdf
            try:
                self.dem = np.genfromtxt(self.dempath)
            except:
                self.dem = np.genfromtxt(self.dempath,skip_header = 6)

            self.nrows = len(self.dem[:,0])
            self.ncols = len(self.dem[0,:])
            blank = np.zeros((self.nrows,self.ncols))

            # Assign some basin-specific things, also needs to be generalized
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
            print('Error reading SNOWAV config file.')

    def process(self):
        '''
        This function calculates everything we will need for the plots.

        Does not currently save to csv...

        '''

        print('SNOWAV processing iSnobal outputs...')

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
        flt_delta_state_byelev = copy.deepcopy(accum_byelev)
        delta_swe_byelev = copy.deepcopy(accum_byelev)
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

            # starting empty string for debug statement
            pf = ''

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
            if pFlag:
                # print('adding precip for %s'%(snow_name))
                for pfile in ppt_files:
                    ppt = ipw.IPW(pfile)
                    pre = ppt.bands[0].data
                    percent_snow = ppt.bands[1].data
                    rain_hrly = rain_hrly + np.multiply(pre,(1-percent_snow))
                    precip_hrly = precip_hrly + pre

            rain_bg = rain_bg + rain_hrly
            precip = precip + precip_hrly
            if pFlag:
                pfs = ', precip added'
            else:
                pfs = ''

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
                pf = ', psnowFile'
                accum_sub_flag = True
                pstate = copy.deepcopy(tmpstate)
                
            # When it is the first flt snow file, copy
            if ((hasattr(self,'flt_flag')) and (snow_name == self.fltpsnowFile)):
                fltpstate = copy.deepcopy(tmpstate)   
            
            # When it is the second flt snow file, copy
            if ((hasattr(self,'flt_flag')) and (snow_name == self.fltcsnowFile)):
                fltcstate = copy.deepcopy(tmpstate) 
                flt_delta_state = fltpstate - fltcstate                               

            # When it hits the current snow file, copy
            if snow_name == self.csnowFile:
                # print('csnowfile is %s'%(snow_name))

                # Run debug statement before ending the process
                pf = ', csnowFile'
                debug = 'snow file: %s, hours: %s, date: %s%s%s'%(
                                    snow_name.split('runs')[1],
                                    str(int(snow_name.split('.')[-1]) - t),
                                    date.date().strftime("%Y-%-m-%-d"),pfs,pf)
                print(debug)

                # Turn off, but last one will still have been added
                accum_sub_flag  = False

                state = copy.deepcopy(tmpstate)
                depth = snow_file.bands[0].data
                density = snow_file.bands[1].data
                self.density = copy.deepcopy(density)
                self.cold = em_file.bands[9].data

                # No need to compile more files after csnowFile
                break

            # It's nice to see where we are...
            debug = 'snow file: %s, hours: %s, date: %s%s%s'%(
                                snow_name.split('runs')[1],
                                str(int(snow_name.split('.')[-1]) - t),
                                date.date().strftime("%Y-%-m-%-d"),pfs,pf)
            print(debug)

            # Step this along so that we can see how many hours between outputs
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
        delta_state = state - pstate

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
            if hasattr(self,'flt_flag'):
                flt_delta_state_byelev_mask = np.multiply(flt_delta_state,mask)
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
                    if hasattr(self,'flt_flag'):
                        flt_delta_state_byelev.loc[b,name] = np.nansum(
                                                flt_delta_state_byelev_mask[ind])                    
                    delta_swe_byelev.loc[b,name] = np.nanmean(
                                                delta_state_byelev_mask[ind])
                else:
                    state_mswe_byelev.loc[b,name] = np.nan
                    depth_mdep_byelev.loc[b,name] = np.nan
                    delta_state_byelev.loc[b,name] = np.nan
                    if hasattr(self,'flt_flag'):
                        flt_delta_state_byelev.loc[b,name] = np.nan
                    delta_swe_byelev.loc[b,name] = np.nan

            self.masks[name]['SWE'] = ( (melt[name].sum() + nonmelt[name].sum())
                                       * self.conversion_factor )

        # Convert to desired units

        # Sum over all time steps, spatial
        self.precip = np.multiply(precip,self.depth_factor)
        self.accum = np.multiply(accum,self.depth_factor)
        self.evap = np.multiply(evap,self.depth_factor)
        self.accum_sub = np.multiply(accum_sub,self.depth_factor)

        # At last time step
        self.state_byelev = np.multiply(state_byelev,self.conversion_factor)
        self.depth_byelev = np.multiply(depth_byelev,self.conversion_factor)
        self.state = np.multiply(state,self.depth_factor)
        self.depth = np.multiply(np.multiply(depth,1000),self.depth_factor)
        self.delta_state = np.multiply(delta_state,self.depth_factor)
        if hasattr(self,'flt_flag'):
            self.flt_delta_state = np.multiply(flt_delta_state,self.depth_factor)

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
        self.flt_delta_state_byelev = np.multiply(flt_delta_state_byelev,
                                              self.conversion_factor)        
        self.delta_swe_byelev = np.multiply(delta_swe_byelev,
                                              self.depth_factor)

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

        # Print some summary info...
        mask = self.state_summary.index.to_series().diff() > pd.Timedelta('24:10:00')
        msum = sum(mask)
        if int(msum) >= 1:
            print('%s entries in dataframe index with gaps larger than 24h '%(msum))
