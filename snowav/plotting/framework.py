
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
import snowav.methods.wyhr_to_datetime as wy
import snowav.utils.get_topo_stats as ts
from snowav.utils.utilities import get_snowav_path
from snowav.utils.OutputReader import iSnobalReader
from inicheck.tools import get_user_config, check_config
from inicheck.output import generate_config, print_config_report
from inicheck.config import MasterConfig
import logging
import coloredlogs
import math


class SNOWAV(object):

    def __init__(self,config_file = None, external_logger=None):
        '''
        Initialize snowav object and read in configuration file.
        '''

        if not os.path.isfile(config_file):
            print('SNOWAV config file does not exist!')
            self.error = True
            return

        print('Reading SNOWAV config file...')
        ucfg = get_user_config(config_file, modules = 'snowav')

        # find path to snowav code directory
        self.snowav_path = get_snowav_path()

        # create blank log and error log because logger is not initialized yet
        self.tmp_log = []
        self.tmp_err = []
        self.tmp_warn = []

        # Check the user config file for errors and report issues if any
        self.tmp_log.append("Checking config file for issues...")
        warnings, errors = check_config(ucfg)
        # print_config_report(warnings, errors)

        ####################################################
        #             snowav system                        #
        ####################################################
        self.loglevel = ucfg.cfg['snowav system']['log_level'].upper()
        self.log_to_file = ucfg.cfg['snowav system']['log_to_file']
        self.basin = ucfg.cfg['snowav system']['basin']
        self.save_path = ucfg.cfg['snowav system']['save_path']
        if self.save_path is None:
            self.save_path = os.path.join(self.snowav_path, 'snowav/data/')
        self.wy = ucfg.cfg['snowav system']['wy']
        self.units = ucfg.cfg['snowav system']['units']
        self.filetype = ucfg.cfg['snowav system']['filetype']
        self.elev_bins = ucfg.cfg['snowav system']['elev_bins']
        if ucfg.cfg['snowav system']['name_append'] != None:
            self.name_append = ucfg.cfg['snowav system']['name_append']
        else:
            self.name_append = '_gen_' + \
                               datetime.datetime.now().strftime("%Y-%-m-%-d")

        ####################################################
        #           outputs                                #
        ####################################################
        self.snowband = ucfg.cfg['outputs']['snowband']
        self.emband = ucfg.cfg['outputs']['emband']
        self.dplcs = ucfg.cfg['outputs']['decimals']
        self.phour = ucfg.cfg['outputs']['phour']
        self.chour = ucfg.cfg['outputs']['chour']

        if (self.phour is not None and self.chour is not None):
            if self.phour >= self.chour:
                self._logger.info('phour > chour, needs to be fixed in config file,'
                                  + ' exiting...')
                print('phour > chour, needs to be fixed in config file, exiting...')
                return

        # Check for forced flight comparison images
        self.fltphour = ucfg.cfg['outputs']['fltphour']
        self.fltchour = ucfg.cfg['outputs']['fltchour']
        if self.fltchour is not None:
            self.flt_flag = True
            if self.fltphour >= self.fltchour:
                self._logger.info('fltphour > fltchour, fix in config file,'
                              + ' exiting...')
                print('fltphour > fltchour, fix in config file, exiting...')
                return
        else:
            self.flt_flag = False

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
            self.basin_total_flag = True
        else:
            self.basin_total_flag = False

        if ucfg.cfg['basin total']['netcdf']:
            self.ncvars = ucfg.cfg['Basin Total']['netcdf'].split(',')
            self.nc_flag = True
        else:
            self.nc_flag = False

        self.flight_dates = ucfg.cfg['basin total']['flights']
        if (self.flight_dates is not None) and (type(self.flight_dates) != list):
            self.flight_dates = [self.flight_dates]

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

        # check paths to see if they need default snowav path
        if self.rep_path is None:
            self.rep_path = os.path.join(self.snowav_path,'snowav/data/')
        if self.env_path is None:
            self.env_path = os.path.join(self.snowav_path,
                                         'snowav/report/template/section_text/')
        if self.templ_path is None:
            self.templ_path = os.path.join(self.snowav_path,'snowav/report/template/')
        if self.summary_file is None:
            self.summary_file = os.path.join(self.snowav_path,
                                             'snowav/report/template/section_text/report_summary.txt')
        if self.tex_file is None:
            self.tex_file = os.path.join(self.snowav_path,
                                         'snowav/report/template/snowav_report.tex')
        if self.figs_tpl_path is None:
            self.figs_tpl_path = os.path.join(self.snowav_path,'snowav/report/figs/')

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

        masks = ucfg.cfg['masks']['basin_masks']
        for idx, m in enumerate(masks):
            maskpaths.append(m)
            self.plotorder.append(ucfg.cfg['masks']['mask_labels'][idx])

        try:
            self.dem = np.genfromtxt(self.dempath)
        except:
            self.dem = np.genfromtxt(self.dempath,skip_header = 6)

        self.nrows = len(self.dem[:,0])
        self.ncols = len(self.dem[0,:])
        blank = np.zeros((self.nrows,self.ncols))

        # Set up processed information
        self.swi = []
        self.evap = []
        self.snowmelt = []
        self.swe = []
        self.depth = []
        self.dates = []
        self.time = []
        self.rho = []
        self.cold = []
        self.rundirs_dict = {}

        for rd in self.run_dirs:
            output = iSnobalReader(rd.split('output')[0],
                                   'netcdf',
                                   snowbands = [0,1,2],
                                   embands = [6,7,8,9],
                                   wy = self.wy)
            self.dates = np.append(self.dates,output.dates)
            self.time = np.append(self.time,output.time)

            # Make a dict for wyhr-rundir lookup
            for t in output.time:
                self.rundirs_dict[t] = rd

            for n in range(0,len(output.em_data[8])):
                self.swi.append(output.em_data[8][n,:,:])
                self.snowmelt.append(output.em_data[7][n,:,:])
                self.evap.append(output.em_data[6][n,:,:])
                self.cold.append(output.em_data[9][n,:,:])
                self.swe.append(output.snow_data[2][n,:,:])
                self.depth.append(output.snow_data[0][n,:,:])
                self.rho.append(output.snow_data[1][n,:,:])

        self.em_files = np.ndarray.tolist(self.time.astype('int'))
        self.snow_files = np.ndarray.tolist(self.time.astype('int'))

        # If no psnowFile and csnowFile specified, use first and last
        if self.chour is not None:
            self.psnowFile = self.phour
            self.csnowFile = self.chour

        else:
            self.phour = self.snow_files[0]
            self.chour = self.snow_files[-1]
            self.psnowFile = self.phour
            self.csnowFile = self.chour
            print('phour and/or chour not specified, using:'
                  + ' %s and %s'%(self.psnowFile,self.csnowFile))

        # Pixel size and elevation bins
        fp = os.path.abspath(self.run_dirs[0].split('output/')[0] + 'snow.nc')
        topo = ts.get_topo_stats(fp,filetype = 'netcdf')
        self.pixel = int(topo['dv'])
        self.edges = np.arange(self.elev_bins[0],
                               self.elev_bins[1]+self.elev_bins[2],
                               self.elev_bins[2])

        # A few remaining basin-specific things
        if self.basin == 'TUOL' or self.basin == 'SJ':
            sr = 6
        else:
            sr = 0

        if self.basin == 'LAKES':
            self.imgx = (1200,1375)
            self.imgy = (425,225)

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
            self.vollbl = self.units
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
        ext_shr = str(self.phour)
        ext_ehr = str(self.chour)
        self.figs_path = os.path.join(self.save_path,
                                    '%s_%s/'%(str(self.phour),str(self.chour)))

        if not os.path.exists(self.figs_path):
            os.makedirs(self.figs_path)

        ####################################################
        #             log file                             #
        ####################################################
        if external_logger == None:
            self.createLog()
        else:
            self._logger = external_logger

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

        self._logger.info('SNOWAV processing iSnobal outputs...')

        cclimit = -5*1000*1000  # based on an average of 60 W/m2 from TL paper

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
        accum_sub_flag = False
        adj = 0 # this is a hack for Hx-repeats-itself forecasting
        t = 0
        for iters,(em_name,snow_name) in enumerate(zip(self.em_files,
                                                       self.snow_files)):
            # iters = 0
            # em_name = self.em_files[iters]
            # snow_name = self.snow_files[iters]
            date = self.dates[iters]

            # starting empty string for debug statement
            pf = ''

            # Hack for Hx-repeats-itself forecasting
            if snow_name == self.psnowFile:
                self.dateFrom = self.dates[iters]
                if self.adj_hours != None:
                    self._logger.debug('Hacking adj_hours...')
                    adj = self.adj_hours

            band = self.swi[iters]
            accum = accum + self.swi[iters]
            daily_snowmelt = self.snowmelt[iters]
            snowmelt = snowmelt + daily_snowmelt
            evap = evap + self.evap[iters]
            tmpstate = self.swe[iters]
            state_byday[:,:,iters] = tmpstate

            # Get rain from input data
            sf = self.rundirs_dict[snow_name].replace('runs','data')
            sf = sf.replace('run','data')
            sf = sf.replace('output','ppt_4b')
            ppt_path = sf.split('em')[0]
            out_hr = snow_name
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

            if pFlag:
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
                               * self.masks[self.plotorder[0]]['mask'])
                               / self.masks[self.plotorder[0]]['mask'].sum()
                               )
                self.pre_swe = (
                                np.nansum(tmpstate
                                * self.masks[self.plotorder[0]]['mask'])
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
                                + copy.deepcopy(daily_snowmelt) )
                precip_sub = precip_sub + precip_hrly

            # When it is the first snow file, copy
            if snow_name == self.psnowFile:
                pf = ', psnowFile'
                accum_sub_flag = True
                pstate = copy.deepcopy(tmpstate)

            # When it is the first flt snow file, copy
            if (self.flt_flag is True) and (snow_name == self.fltphour):
                fltpstate = copy.deepcopy(tmpstate)

            # When it is the second flt snow file, copy
            if (self.flt_flag is True) and (snow_name == self.fltchour):
                fltcstate = copy.deepcopy(tmpstate)
                flt_delta_state = fltpstate - fltcstate

            # When it hits the current snow file, copy
            if snow_name == self.csnowFile:
                self.dateTo = self.dates[iters]

                # Run debug statement before ending the process
                pf = ', csnowFile'

                debug = 'snow file: %s, hours: %s, date: %s%s%s'%(
                                    self.rundirs_dict[int(snow_name)],
                                    str(int(snow_name) - t),
                                    date.date().strftime("%Y-%-m-%-d"),pfs,pf)
                self._logger.debug(debug)

                # Turn off, but last one will still have been added
                accum_sub_flag  = False

                state = copy.deepcopy(tmpstate)
                depth = self.depth[iters]
                self.density = copy.deepcopy(self.rho[iters])
                self.cold = self.cold[iters]

                # No need to compile more files after csnowFile
                break

            # It's nice to see where we are...
            debug = 'snow file: %s, hours: %s, date: %s%s%s'%(
                                    self.rundirs_dict[int(snow_name)],
                                    str(int(snow_name) - t),
                                    date.date().strftime("%Y-%-m-%-d"),pfs,pf)

            self._logger.debug(debug)

            # Step this along so that we can see how many hours between outputs
            t = int(snow_name)

        # Append date to report name
        parts = self.report_name.split('.')
        self.report_name = ( parts[0]
                           + self.dateTo.date().strftime("%Y%m%d")
                           + '.' + parts[1] )

        # Difference in state (SWE)
        delta_state = state - pstate

        # Mask by subbasin and elevation band
        for name in self.masks:
            # name = self.plotorder[0]
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
            if self.flt_flag is True:
                flt_delta_state_byelev_mask = np.multiply(flt_delta_state,mask)
            state_byelev_mask = np.multiply(state,mask)
            state_mswe_byelev_mask = np.multiply(state,mask)
            density_m_byelev_mask = np.multiply(self.density,mask)

            ix = density_m_byelev_mask == 0
            density_m_byelev_mask[ix] = np.nan

            depth_mdep_byelev_mask = np.multiply(depth,mask)
            elevbin = np.multiply(self.ixd,mask)

            # Do it by elevation band
            for n in np.arange(0,len(self.edges)):
                # n = 0
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
                    if self.flt_flag is True:
                        flt_delta_state_byelev.loc[b,name] = np.nansum(
                                                flt_delta_state_byelev_mask[ind])
                    delta_swe_byelev.loc[b,name] = np.nanmean(
                                                delta_state_byelev_mask[ind])
                else:
                    state_mswe_byelev.loc[b,name] = np.nan
                    depth_mdep_byelev.loc[b,name] = np.nan
                    delta_state_byelev.loc[b,name] = np.nan
                    if self.flt_flag is True:
                        flt_delta_state_byelev.loc[b,name] = np.nan
                    delta_swe_byelev.loc[b,name] = np.nan

                # There are apparently side cases with density and not swe
                if np.sum(np.sum(density_m_byelev_mask[ind])) > 0:
                    density_m_byelev.loc[b,name] = np.nanmean(
                                                    copy.deepcopy(density_m_byelev_mask[ind]))
                else:
                    density_m_byelev.loc[b,name] = np.nan

            self.masks[name]['SWE'] = ( (melt[name].sum() + nonmelt[name].sum())
                                       * self.conversion_factor )


        # First time step, volume
        self.pstate = np.multiply(pstate,self.conversion_factor)

        # At last time step, volume
        self.state_byelev = np.multiply(state_byelev,self.conversion_factor)
        # self.depth_byelev = np.multiply(depth_byelev,self.conversion_factor)

        # Sum over all time steps, depth
        self.precip = np.multiply(precip,self.depth_factor)
        self.accum = np.multiply(accum,self.depth_factor)
        self.evap = np.multiply(evap,self.depth_factor)
        self.accum_sub = np.multiply(accum_sub,self.depth_factor)

        # Sum over all time steps, volume
        self.accum_byelev = np.multiply(accum_byelev,self.conversion_factor)
        self.evap_byelev = np.multiply(evap_byelev,self.conversion_factor)
        self.precip_byelev = np.multiply(precip_byelev,self.conversion_factor)
        self.rain_bg_byelev = np.multiply(rain_bg_byelev,
                                          self.conversion_factor)
        self.snowmelt_byelev = np.multiply(snowmelt_byelev,
                                           self.conversion_factor)

        # At last time step, depth
        self.depth = np.multiply(np.multiply(depth,1000),self.depth_factor)
        self.state = np.multiply(state,self.depth_factor)
        self.delta_state = np.multiply(delta_state,self.depth_factor)
        if self.flt_flag is True:
            self.flt_delta_state = np.multiply(flt_delta_state,self.depth_factor)

        # Sub-set of time defined by psnowfile and csnowfile
        self.snowmelt_byelev_sub = np.multiply(snowmelt_byelev_sub,
                                               self.conversion_factor)
        self.accum_byelev_sub = np.multiply(accum_byelev_sub,
                                            self.conversion_factor)
        self.precip_byelev_sub = np.multiply(precip_byelev_sub,
                                             self.conversion_factor)

        # Change
        self.delta_state_byelev = np.multiply(delta_state_byelev,
                                              self.conversion_factor)
        self.flt_delta_state_byelev = np.multiply(flt_delta_state_byelev,
                                              self.conversion_factor)
        self.delta_swe_byelev = np.multiply(delta_swe_byelev,
                                              self.depth_factor)

        # Daily
        self.state_summary = np.multiply(state_summary,self.conversion_factor)
        self.accum_summary = np.multiply(accum_summary,self.conversion_factor)
        self.precip_summary = np.multiply(precip_summary,self.conversion_factor)
        self.evap_summary = np.multiply(evap_summary,self.conversion_factor)
        # self.state_byday = np.multiply(state_byday,self.depth_factor)

        # Mean
        self.state_mswe_byelev = np.multiply(state_mswe_byelev,
                                             self.depth_factor)
        self.depth_mdep_byelev = np.multiply(np.multiply(depth_mdep_byelev,1000),
                                             self.depth_factor)
        self.density_m_byelev = density_m_byelev

        # odds and ends
        self.melt = np.multiply(melt,self.conversion_factor)
        self.nonmelt = np.multiply(nonmelt,self.conversion_factor)
        self.cold = np.multiply(self.cold,0.000001)

        # Print some summary info...
        mask = self.state_summary.index.to_series().diff() > pd.Timedelta('24:10:00')
        msum = sum(mask)
        if int(msum) >= 1:
            self._logger.debug('%s entries in dataframe index with gaps larger than 24h '%(msum))

    def createLog(self):
        '''
        Now that the directory structure is done, create log file and print out
        saved logging statements.
        '''

        level_styles = {'info': {'color': 'white'},
                        'notice': {'color': 'magenta'},
                        'verbose': {'color': 'blue'},
                        'success': {'color': 'green', 'bold': True},
                        'spam': {'color': 'green', 'faint': True},
                        'critical': {'color': 'red', 'bold': True},
                        'error': {'color': 'red'},
                        'debug': {'color': 'green'},
                        'warning': {'color': 'yellow'}}

        field_styles =  {'hostname': {'color': 'magenta'},
                         'programname': {'color': 'cyan'},
                         'name': {'color': 'white'},
                         'levelname': {'color': 'white', 'bold': True},
                         'asctime': {'color': 'green'}}

        # start logging
        loglevel = self.loglevel

        numeric_level = getattr(logging, loglevel, None)
        if not isinstance(numeric_level, int):
            raise ValueError('Invalid log level: %s' % loglevel)

        # setup the logging
        logfile = None
        if self.log_to_file:
            logfile = os.path.join(self.figs_path, 'log_snowav.out')
            # let user know
            print('Logging to file: {}'.format(logfile))

        fmt = '%(levelname)s:%(name)s:%(message)s'
        if logfile is not None:
            logging.basicConfig(filename=logfile,
                                filemode='w',
                                level=numeric_level,
                                format=fmt)
        else:
            logging.basicConfig(level=numeric_level)
            coloredlogs.install(level=numeric_level,
                                fmt=fmt,
                                level_styles=level_styles,
                                field_styles=field_styles)

        self._loglevel = numeric_level

        self._logger = logging.getLogger(__name__)

        # print title and mountains
        # title, mountain = self.title()
        # for line in mountain:
        #     self._logger.info(line)
        # for line in title:
        #     self._logger.info(line)
        # dump saved logs
        if len(self.tmp_log) > 0:
            for l in self.tmp_log:
                self._logger.info(l)
        if len(self.tmp_warn) > 0:
            for l in self.tmp_warn:
                self._logger.warning(l)
        if len(self.tmp_err) > 0:
            for l in self.tmp_err:
                self._logger.error(l)
