
import numpy as np
from matplotlib import pyplot as plt
from mpl_toolkits.axes_grid1 import make_axes_locatable
import seaborn as sns
from smrf import ipw
from shutil import copyfile
import os
import copy
import pandas as pd
import sys
import datetime
from snowav.utils.OutputReader import iSnobalReader
from inicheck.tools import get_user_config
from inicheck.output import generate_config


class SNOWAV(object):

    def __init__(self,config_file = None):
        '''
        Initialize snowav object and read in configuration file.
        '''

        if not os.path.isfile(config_file):
            print('SNOWAV config file does not exist!')
            self.error = True
            return

        print('Reading SNOWAV config file...')
        ucfg = get_user_config(config_file, modules = 'snowav')
    
        ####################################################
        #             basin                                #
        ####################################################
        self.basin = ucfg.cfg['basin']['basin']
        self.save_path = ucfg.cfg['basin']['save_path']

        if ucfg.cfg['basin']['name_append'] != None:
            self.name_append = ucfg.cfg['basin']['name_append']
        else:
            self.name_append = '_gen_' + \
                               datetime.datetime.now().strftime("%Y-%-m-%-d")

        self.wy = ucfg.cfg['basin']['wy']

        self.units = ucfg.cfg['basin']['units']
        self.filetype = ucfg.cfg['basin']['filetype']

        ####################################################
        #           outputs                                #
        ####################################################
        self.snowband = ucfg.cfg['outputs']['snowband']
        self.emband = ucfg.cfg['outputs']['emband']

        # Rounding decimals for volume and depth
        self.dplcs = ucfg.cfg['outputs']['decimals']

        # If previous snowFile and current snowFile are specified, use those,
        # otherwise they will get defined as the first and last
        # snow.XXXX file once [Runs] has been compiled
        if (ucfg.cfg['outputs']['psnowfile'] != None and
            ucfg.cfg['outputs']['csnowfile'] != None):

            self.psnowFile = ucfg.cfg['outputs']['psnowfile']
            self.csnowFile = ucfg.cfg['outputs']['csnowfile']
            self.cemFile = self.csnowFile.replace('snow.','em.')

        # Check for forced flight comparison images
        if (ucfg.cfg['outputs']['fltpsnowfile'] != None and
            ucfg.cfg['outputs']['fltcsnowfile'] != None):

            self.fltpsnowFile = ucfg.cfg['outputs']['fltpsnowFile']
            self.fltcsnowFile = ucfg.cfg['outputs']['fltcsnowFile']
            self.flt_flag = True

        self.summary = ucfg.cfg['outputs']['summary']
        if type(self.summary) != list:
            self.summary = [self.summary]

        ####################################################
        #           runs                                   #
        ####################################################
        self.run_dirs = ucfg.cfg['runs']['run_dirs']
        if type(self.run_dirs) != list:
            self.run_dirs = [self.run_dirs]

        ####################################################
        #           validate                               #
        ####################################################
        if (ucfg.cfg['validate']['stations'] != None and
            ucfg.cfg['validate']['labels'] != None and
            ucfg.cfg['validate']['client'] != None):

            self.val_stns = ucfg.cfg['validate']['stations']
            self.val_lbls = ucfg.cfg['validate']['labels']
            self.val_client = ucfg.cfg['validate']['client']
            self.valid_flag = True

        else:
            self.valid_flag = False
            print('No validation stations listed, will not generate figure')

        # This is being used to combine 2017 HRRR data
        self.offset = int(ucfg.cfg['validate']['offset'])

        ####################################################
        #           basin total                            #
        ####################################################
        if ucfg.cfg['basin total']['summary_swe'] != None:
            self.summary_swe = ucfg.cfg['basin total']['summary_swe']
            self.summary_swi = ucfg.cfg['basin total']['summary_swi']
        else:
            self.basin_total_flag = False

        if ucfg.cfg['basin total']['netcdf']:
            self.ncvars = ucfg.cfg['Basin Total']['netcdf'].split(',')
            self.nc_flag = True
        else:
            self.nc_flag = False

        ####################################################
        #           masks                                  #
        ####################################################
        self.dempath = ucfg.cfg['masks']['dempath']
        self.total = ucfg.cfg['masks']['basin_masks'][0]

        ####################################################
        #          plots                                   #
        ####################################################
        self.figsize = (ucfg.cfg['plots']['fig_length'],
                        ucfg.cfg['plots']['fig_height'])
        self.dpi = ucfg.cfg['plots']['dpi']
        self.barcolors = ['xkcd:true green','palegreen', 'xkcd:dusty green',
                          'xkcd:vibrant green','red']

        ####################################################
        #          report                                  #
        ####################################################
        self.report_flag = ucfg.cfg['report']['report']

        self.exclude_figs = ucfg.cfg['report']['exclude_figs']
        if type(self.exclude_figs) != list and self.exclude_figs != None:
            self.exclude_figs = [self.exclude_figs]

        self.report_name = ucfg.cfg['report']['report_name']
        self.rep_title = ucfg.cfg['report']['report_title']
        self.rep_path = ucfg.cfg['report']['rep_path']
        self.env_path = ucfg.cfg['report']['env_path']
        self.templ_path = ucfg.cfg['report']['templ_path']
        self.tex_file = ucfg.cfg['report']['tex_file']
        self.summary_file = ucfg.cfg['report']['summary_file']
        self.figs_tpl_path = ucfg.cfg['report']['figs_tpl_path']

        ####################################################
        #           hx forecast
        ####################################################
        self.adj_hours = ucfg.cfg['hx forecast']['adj_hours']

        ####################################################
        #           masks
        ####################################################
        for item in ['basin_masks', 'mask_labels']:
            if type(ucfg.cfg['masks'][item]) != list:
                ucfg.cfg['masks'][item] = [ucfg.cfg['Masks'][item]]

        self.plotorder = []
        maskpaths = []

        # List of paths to masks. The first should always be the total basin
        masks = ucfg.cfg['masks']['basin_masks']
        for idx, m in enumerate(masks):
            maskpaths.append(m)
            self.plotorder.append(ucfg.cfg['masks']['mask_labels'][idx])

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

        if self.filetype == 'ipw':
            self.snow_files = []
            self.em_files = []
            for rdir in self.run_dirs:
                run_files = [rdir + s for s in sorted(os.listdir(rdir))]

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
                print('psnowfile and/or csnowfile not specified, using:'
                      + ' \n%s and \n%s'%(self.psnowFile,self.csnowFile))

        # Currently this is a dead end, and process() will fail
        if self.filetype == 'netcdf':
            self.swi = []
            self.evap = []
            self.swe = []
            self.depth = []
            self.dates = []
            self.time = []
            self.rho = []

            for rd in self.run_dirs:
                output = iSnobalReader(rd.split('output')[0],
                                       'netcdf',
                                       snowbands = [0,1,2],
                                       embands = [7,8])
                self.dates = np.append(self.dates,output.dates)
                self.time = np.append(self.time,output.time)

                for n in range(0,len(output.em_data[8])):
                    self.swi.append(output.em_data[8][n,:,:])
                    self.evap.append(output.em_data[7][n,:,:])
                    self.swe.append(output.snow_data[2][n,:,:])
                    self.depth.append(output.snow_data[0][n,:,:])
                    self.rho.append(output.snow_data[1][n,:,:])

        # Assign some basin-specific things, also needs to be generalized
        if self.basin == 'BRB':
            self.pixel = 100
            sr = 0
            if self.units == 'KAF':
                emin = 2000
                emax = 11000
                self.step = 1000
            if self.units == 'SI':
                emin = 800
                emax = 3200
                self.step = 250
        if self.basin == 'TUOL':
            self.pixel = 50
            sr = 6
            if self.units == 'KAF':
                emin = 3000      
                emax = 12000
                self.step   = 1000
            if self.units == 'SI':
                emin = 800       
                emax = 3600
                self.step   = 500
        if self.basin == 'SJ':
            self.pixel = 50
            sr = 6
            if self.units == 'KAF':
                emin = 1000      
                emax  = 13000
                self.step = 1000
            if self.units == 'SI':
                emin = 300       
                emax = 4000
                self.step = 500
        if self.basin == 'LAKES':
            self.pixel = 50
            sr = 0
            self.imgx = (1200,1375)
            self.imgy = (425,225)
            if self.units == 'KAF':
                emin = 8000     
                emax = 12500
                self.step = 500
            if self.units == 'SI':
                emin = 2400       
                emax = 3800
                self.step = 200
        if self.basin == 'RCEW':
            self.pixel = 50
            sr = 0
            if self.units == 'KAF':
                emin = 2500     
                emax = 7500
                self.step = 500
            if self.units == 'SI':
                emin = 800      
                emax = 2500
                self.step = 500

        # Create the elevation bins
        self.edges = np.arange(emin,emax+self.step,self.step)
        # Right now this is a placeholder, could edit by basin...
        self.xlims = (0,len(self.edges))

        # Compile the masks, need to streamline handling for Willow Creek
        try:
            self.masks = dict()
            for lbl,mask in zip(self.plotorder,maskpaths):
                if (self.basin == 'SJ' and lbl == 'Willow Creek'):
                    self.masks[lbl] = {'border': blank,
                                       'mask': np.genfromtxt(mask,skip_header=0),
                                       'label': lbl}
                else:
                    self.masks[lbl] = {'border': blank,
                                       'mask': np.genfromtxt(mask,skip_header=sr),
                                       'label': lbl}
        except:
            print('Failed creating mask dicts..')
            self.error = True
            return

        # Assign unit-specific things
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
        self.figs_path = os.path.join(self.save_path,
                                      '%s_%s/'%(ext_shr[1][1:5],ext_ehr[1][1:5]))

        if not os.path.exists(self.figs_path):
            os.makedirs(self.figs_path)

        # Only need to store this name if we decide to
        # write more to the copied config file...
        self.config_copy = (self.figs_path
                            + extf[0]
                            + self.name_append
                            + '_%s_%s'%(ext_shr[1][1:5],ext_ehr[1][1:5])
                            + extf[1])

        if not os.path.isfile(self.config_copy):
            generate_config(ucfg,self.config_copy)

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
