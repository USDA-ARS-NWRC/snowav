
import numpy as np
from smrf import ipw
from shutil import copyfile
import os
import copy
import pandas as pd
import sys
import datetime
import snowav.utils.wyhr_to_datetime as wy
import snowav.utils.get_topo_stats as ts
from snowav.utils.utilities import get_snowav_path
from snowav.utils.OutputReader import iSnobalReader
import logging
import coloredlogs
import math
import netCDF4 as nc
import mysql.connector

def process(self):
    '''
    This function processes and summarizes iSnobal outputs.

    '''

    self._logger.info('SNOWAV processing iSnobal outputs...')

    # Is this a good idea?
    pd.options.mode.chained_assignment = None

    # based on an average of 60 W/m2 from TL paper
    cclimit = -5*1000*1000

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
    state_byday = np.zeros((self.nrows,self.ncols,len(self.outputs['dates'])))

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
    for iters,snow_name in enumerate(self.outputs['dates']):

        date = self.outputs['dates'][iters]

        # starting empty string for debug statement
        pf = ''

        # Hack for Hx-repeats-itself forecasting
        if snow_name == self.psnowFile:
            self.dateFrom = self.outputs['dates'][iters]
            if self.adj_hours != None:
                self._logger.debug('Hacking adj_hours...')
                adj = self.adj_hours

        band = self.outputs['swi'][iters]
        ixn = np.isnan(band)
        band[ixn] = 0
        accum = accum + band

        daily_snowmelt = self.outputs['snowmelt'][iters]
        ixn = np.isnan(daily_snowmelt)
        daily_snowmelt[ixn] = 0
        snowmelt = snowmelt + daily_snowmelt

        daily_evap = self.outputs['evap'][iters]
        ixn = np.isnan(daily_evap)
        daily_evap[ixn] = 0
        evap = evap + daily_evap

        tmpstate = self.outputs['swe'][iters]
        state_byday[:,:,iters] = tmpstate

        # Get rain from input data
        hr = int(self.outputs['time'][iters])

        sf = self.rundirs_dict[hr].replace('runs','data')
        sf = sf.replace('run','data')
        sf = sf.replace('output','ppt_4b')
        ppt_path = sf.split('em')[0]
        out_hr = hr
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
            self.fltdateFrom = self.outputs['dates'][iters]

        # When it is the second flt snow file, copy
        if (self.flt_flag is True) and (snow_name == self.fltchour):
            fltcstate = copy.deepcopy(tmpstate)
            flt_delta_state = fltpstate - fltcstate
            self.fltdateTo = self.outputs['dates'][iters]

        # When it hits the current snow file, copy
        if snow_name == self.csnowFile:
            self.dateTo = self.outputs['dates'][iters]

            # Run debug statement before ending the process
            pf = ', csnowFile'

            debug = 'snow file: %s, hours: %s, date: %s%s%s'%(
                                self.rundirs_dict[hr],
                                str(hr - t),
                                date.date().strftime("%Y-%-m-%-d"),pfs,pf)
            self._logger.debug(debug)

            # Turn off, but last one will still have been added
            accum_sub_flag  = False

            state = copy.deepcopy(tmpstate)
            depth = self.outputs['depth'][iters]
            self.density = copy.deepcopy(self.outputs['rho'][iters])
            self.cold = self.outputs['coldcont'][iters]

            # No need to compile more files after csnowFile
            break

        # It's nice to see where we are...
        debug = 'snow file: %s, hours: %s, date: %s%s%s'%(
                                self.rundirs_dict[hr],
                                str(hr - t),
                                date.date().strftime("%Y-%-m-%-d"),pfs,pf)

        self._logger.debug(debug)

        # Step this along so that we can see how many hours between outputs
        t = int(hr)

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
