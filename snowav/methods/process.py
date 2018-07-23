
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
from snowav import database

def process(self):
    '''
    This function processes and summarizes iSnobal outputs.

    to-dos:
    -round
    - basin total on database for depths isn't correct

    images:
    - swi, period summary
    - cold content, end of period
    - swe, end of period
    - swe, change during period
    - swe, flt change
    - density, end of period

    '''

    self._logger.info('SNOWAV processing iSnobal outputs...')

    # Is this a good idea?
    pd.options.mode.chained_assignment = None

    # based on an average of 60 W/m2 from TL paper
    cclimit = -5*1000*1000
    # self.cold = np.multiply(self.cold,0.000001)

    # Daily, by elevation
    dz = pd.DataFrame(np.nan, index = self.edges, columns = self.masks.keys())

    # These get treated different for basin totals than volumes do
    mean_fields = ['swe_z','swi_z','precip_z','depth','density','rain_z',
                   'evap_z','coldcont']

    adj = 0
    t = 0
    for iters, out_date in enumerate(self.outputs['dates']):
        daily_outputs = {'swe_vol':dz.copy(),
                         'swe_avail':dz.copy(),
                         'swe_unavail':dz.copy(),
                         'swe_z':dz.copy(),
                         'swi_vol':dz.copy(),
                         'swi_z':dz.copy(),
                         'precip_vol':dz.copy(),
                         'precip_z':dz.copy(),
                         'depth':dz.copy(),
                         'density':dz.copy(),
                         'rain_z':dz.copy(),
                         'evap_z':dz.copy(),
                         'coldcont':dz.copy()}

        # starting empty string for debug statement
        pf = ''
        pfs = ''

        # Hack for Hx-repeats-itself forecasting
        if out_date == self.psnowFile:
            self.dateFrom = self.outputs['dates'][iters]
            if self.adj_hours != None:
                self._logger.debug('Hacking adj_hours...')
                adj = self.adj_hours

        self.cold = np.multiply(self.outputs['coldcont'][iters],0.000001)

        # Get daily rain from hourly input data
        hr = int(self.outputs['time'][iters])
        sf = self.rundirs_dict[hr].replace('runs','data')
        sf = sf.replace('run','data')
        sf = sf.replace('output','ppt_4b')
        ppt_path = sf.split('em')[0]
        hrs = range(hr - 23, hr + 1)

        pFlag = False
        ppt_files = []
        for hr in hrs:
            # Make 0-padded strings
            if hr < 100:
                shr = '00' + str(hr)
            if hr >= 100 and hr < 1000:
                shr = '0' + str(hr)
            if hr >= 1000:
                shr = str(hr)

            if os.path.isfile(ppt_path + 'ppt.4b_' + shr):
                pFlag = True
                ppt_files = ppt_files + [ppt_path + 'ppt.4b_' + shr]

        precip = np.zeros((self.nrows,self.ncols))
        rain = np.zeros((self.nrows,self.ncols))

        if pFlag:
            for pfile in ppt_files:
                ppt = ipw.IPW(pfile)
                pre = ppt.bands[0].data
                percent_snow = ppt.bands[1].data
                rain = rain + np.multiply(pre,(1-percent_snow))
                precip = precip + pre

        if pFlag:
            pfs = ', precip added'

        # When it is the first snow file, copy
        if out_date == self.psnowFile:
            pf = ', psnowFile'
            pstate = copy.deepcopy(self.outputs['swe_z'][iters])

        # When it is the first flt snow file, copy
        if (self.flt_flag is True) and (out_date == self.fltphour):
            fltpstate = copy.deepcopy(self.outputs['swe_z'][iters])
            self.fltdateFrom = self.outputs['dates'][iters]

        # When it is the second flt snow file, copy
        if (self.flt_flag is True) and (out_date == self.fltchour):
            fltcstate = copy.deepcopy(self.outputs['swe_z'][iters])
            flt_delta_state = fltpstate - fltcstate
            self.fltdateTo = self.outputs['dates'][iters]

        # When it hits the current snow file, copy
        if out_date == self.csnowFile:
            self.dateTo = self.outputs['dates'][iters]

            depth = self.outputs['depth'][iters]
            self.density = copy.deepcopy(self.outputs['density'][iters])
            self.cold = self.outputs['coldcont'][iters]

            # Run debug statement before ending the process
            pf = ', csnowFile'
            debug = 'snow file: %s, hours: %s, date: %s%s%s'%(
                                self.rundirs_dict[hr],
                                str(hr - t),
                                out_date.date().strftime("%Y-%-m-%-d"),pfs,pf)
            self._logger.debug(debug)

        # Get a snow-free mask ready
        swe = copy.deepcopy(self.outputs['swe_z'][iters])

        # Go through each outputs; some are available from self.outputs,
        # some are calculated
        for k in list(daily_outputs.keys()):
            # Mask by subbasin

            for name in self.masks:

                mask = copy.deepcopy(self.masks[name]['mask'])
                elevbin = np.multiply(self.ixd,mask)

                # If the key is something we've already calculated (i.e., not
                # volume), get the masked subbasin output, otherwise it will be
                # calculated below
                if k in self.outputs.keys():
                    mask_out = np.multiply(self.outputs[k][iters],mask)

                # make a subbasin, by elevation, snow-free mask
                swe_mask_sub = np.multiply(copy.deepcopy(swe),mask)

                # only need this for swe and cold content
                if k in ['swe_z','coldcont','density']:
                    cold = np.multiply(self.outputs[k][iters],mask)

                # Now mask the subbasins by elevation band
                for n in np.arange(0,len(self.edges)):
                    ind = elevbin == n
                    b = self.edges[n]

                    # index for pixels with SWE, in subbasin, in elevation band
                    ix = swe_mask_sub[ind] > 0

                    # Assign values
                    # swe_z and derivatives
                    if k == 'swe_z':
                        ccb = cold[ind]
                        cind = ccb > cclimit

                        # masked out by pixels with snow -> [ix]
                        swe_bin = mask_out[ind]
                        if swe_bin.size:
                            r = np.nansum(swe_bin[cind])
                            ru = np.nansum(swe_bin[~cind])

                            daily_outputs['swe_unavail'].loc[b,name] = (
                                        np.multiply(ru,self.conversion_factor) )
                            daily_outputs['swe_avail'].loc[b,name] = (
                                        np.multiply(r,self.conversion_factor) )
                            daily_outputs['swe_z'].loc[b,name] = (
                                        np.multiply(np.nanmean(swe_bin),
                                                    self.depth_factor) )
                            daily_outputs['swe_vol'].loc[b,name] = (
                                        np.multiply(np.nansum(swe_bin),
                                                    self.conversion_factor) )

                    elif k == 'swi_z':
                        # Not masked out by pixels with snow
                        r = np.nansum(mask_out[ind])

                        daily_outputs['swi_z'].loc[b,name] = (
                                        np.multiply(np.nanmean(mask_out[ind]),
                                                    self.depth_factor) )
                        daily_outputs['swi_vol'].loc[b,name] = (
                                        np.multiply(r,self.conversion_factor) )

                    elif k in ['evap_z','depth']:
                        # masked out by pixels with snow -> [ix]
                        r = np.nanmean(mask_out[ind][ix])

                        # Depth is in m (swi/precip/etc are mm)
                        if k == 'depth':
                            daily_outputs[k].loc[b,name] = (
                                        np.multiply(r,self.depth_factor) )

                        # evap_z
                        daily_outputs[k].loc[b,name] = (
                                        np.multiply(r,self.depth_factor) )

                    elif k in ['coldcont','density']:
                        # masked out by pixels with snow -> [ix]
                        daily_outputs[k].loc[b,name] = np.nanmean(mask_out[ind][ix])

                    elif k == 'precip_z':
                        # not masked out by pixels with snow
                        mask_out = np.multiply(precip,mask)
                        r = np.nanmean(mask_out[ind])
                        rv = np.nansum(mask_out[ind])

                        daily_outputs[k].loc[b,name] = (
                                        np.multiply(r,self.depth_factor) )
                        daily_outputs['precip_vol'].loc[b,name] = (
                                        np.multiply(rv,self.conversion_factor) )

                    elif k == 'rain_z':
                        # not masked out by pixels with snow
                        mask_out = np.multiply(rain,mask)
                        r = np.nanmean(mask_out[ind])
                        daily_outputs[k].loc[b,name] = (
                                        np.multiply(r,self.depth_factor) )

                # calculate basin total mean/volume
                # previously was looping over mean_fields...
                # if k in mean_fields:
                ixs = swe_mask_sub > 0

                if k in ['swe_z','evap_z','coldcont']:
                    # Mask by snow-free
                    out = np.multiply(self.outputs[k][iters],mask)
                    total = np.multiply(np.nanmean(out[ixs]), self.depth_factor)
                    daily_outputs[k].loc['total',name] = copy.deepcopy(total)

                    if k == 'swe_z':
                        # calc swe_z derived values
                        s = np.nansum(daily_outputs['swe_vol'][name].values)
                        s1 = np.nansum(daily_outputs['swe_avail'][name].values)
                        s2 = np.nansum(daily_outputs['swe_unavail'][name].values)
                        daily_outputs['swe_vol'].loc['total',name] = copy.deepcopy(s)
                        daily_outputs['swe_avail'].loc['total',name] = copy.deepcopy(s1)
                        daily_outputs['swe_unavail'].loc['total',name] = copy.deepcopy(s2)

                if k == 'depth':
                    # Mask by snow-free
                    out = np.multiply(self.outputs[k][iters],mask)
                    total = np.multiply(np.nanmean(out[ixs]), self.depth_factor)
                    daily_outputs[k].loc['total',name] = copy.deepcopy(total)

                if k == 'density':
                    # Mask by snow-free
                    out = np.multiply(self.outputs[k][iters],mask)
                    total = np.nanmean(out[ixs])
                    daily_outputs[k].loc['total',name] = copy.deepcopy(total)

                if k == 'swi_z':
                    # Do not mask by snow free
                    out = np.multiply(self.outputs[k][iters],mask)
                    total = np.multiply(np.nanmean(out), self.depth_factor)
                    daily_outputs[k].loc['total',name] = copy.deepcopy(total)

                    s = np.nansum(daily_outputs['swi_vol'][name].values)
                    daily_outputs['swi_vol'].loc['total',name] = copy.deepcopy(s)

                if k == 'precip_z':
                    out = np.multiply(precip,mask)
                    total = np.multiply(np.nanmean(out), self.depth_factor)
                    daily_outputs[k].loc['total',name] = copy.deepcopy(total)

                    s = np.nansum(daily_outputs['precip_vol'][name].values)
                    daily_outputs['precip_vol'].loc['total',name] = copy.deepcopy(s)

                if k == 'rain_z':
                    out = np.multiply(rain,mask)
                    total = np.multiply(np.nanmean(out), self.depth_factor)
                    daily_outputs[k].loc['total',name] = copy.deepcopy(total)

                # k is never called for these!
                # if k in ['swe_vol','swe_avail','swe_unavail','swi_vol','precip_vol']:
                #     s = np.nansum(daily_outputs[k][name].values)
                #     daily_outputs[k].loc['total',name] = copy.deepcopy(s)

            # Send daily results to database
            if self.location == 'database':
                df = daily_outputs[k].copy()
                database.package_results.package_results(self, df, k,
                                                         self.outputs['dates'][iters])

        # It's nice to see where we are...
        debug = 'snow file: %s, hours: %s, date: %s%s%s'%(
                                self.rundirs_dict[hr],
                                str(hr - t),
                                out_date.date().strftime("%Y-%-m-%-d"),pfs,pf)

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

    # Print some summary info...
    mask = self.state_summary.index.to_series().diff() > pd.Timedelta('24:10:00')
    msum = sum(mask)
    if int(msum) >= 1:
        self._logger.debug('%s entries in dataframe index with gaps larger than 24h '%(msum))
