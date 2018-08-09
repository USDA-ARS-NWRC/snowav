
import numpy as np
from smrf import ipw
import os
import copy
import pandas as pd
import datetime
import logging
import coloredlogs
import netCDF4 as nc
from snowav import database
import warnings

def process(self):
    '''
    This function processes and summarizes iSnobal outputs.

    It loops over each day in the specified directory and:
        - sums total precip and rain from the ipw ppt ppt files
        - loops over each daily_outputs ('swe_z', 'swe_vol', etc)
        - loops over sub-basins
        - loops over elevation band, and calculates volumes and mean depths
            (masked where there is SWE greater than 0) for each daily_outputs
            value in that elevation band
        - once all elevation bands for each subbasin have been calculated,
            'total' fields (either mean depth for the subbasin or total volume
            for the subbasin) are calculated
        - finally, that is sent to the database

    to-dos:
    - Hx-Repeats-Itself forecasting probably no longer works in this format

    '''

    self._logger.info('SNOWAV processing iSnobal outputs...')

    # Is this a good idea?
    pd.options.mode.chained_assignment = None

    # This suppresses the "mean of empty slice" warning
    warnings.filterwarnings('ignore')

    # based on an average of 60 W/m2 from TL paper
    cclimit = -5*1000*1000

    # Daily, by elevation
    dz = pd.DataFrame(0, index = self.edges, columns = self.masks.keys())

    # Order matters!
    vars = ['coldcont','evap_z','rain_z','density','depth','precip_z',
            'precip_vol','swi_z','swi_vol','swe_z','swe_vol','swe_avail',
            'swe_unavail']

    adj = 0
    t = 0
    for iters, out_date in enumerate(self.outputs['dates']):

        # Initialize output dataframes for every day
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

        # Hack for Hx-repeats-itself forecasting
        if out_date == self.start_date:
            self.dateFrom = self.outputs['dates'][iters]
            if self.adj_hours != None:
                self._logger.debug('Hacking adj_hours...')
                adj = self.adj_hours

        # Get daily rain from hourly input data
        hr = int(self.outputs['time'][iters])
        sf = self.rundirs_dict[hr].replace('runs','data')
        sf = sf.replace('run','data')
        ppt_path = sf.replace('output','ppt_4b')
        hrs = range(hr - 23, hr + 1)

        pFlag = False
        ppt_files = []
        # Get all the ppt_files for that day and make 0-padded strings
        for hr in hrs:
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

        # Get a snow-free mask ready
        swe = copy.deepcopy(self.outputs['swe_z'][iters])
        cold = copy.deepcopy(self.outputs['coldcont'][iters])

        # Loop over outputs (depths are copied, volumes are calculated)
        for k in vars:

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
                cold_sub = np.multiply(cold,mask)

                # Now mask the subbasins by elevation band
                for n in np.arange(0,len(self.edges)):
                    ind = elevbin == n
                    b = self.edges[n]

                    # index for pixels with SWE, in subbasin, in elevation band,
                    # needed for mean values
                    ix = swe_mask_sub[ind] > 0

                    # Assign values
                    if k == 'swe_z':
                        mask_out = np.multiply(self.outputs['swe_z'][iters],mask)
                        ccb = cold_sub[ind][ix]
                        cind = ccb > cclimit

                        # masked out by pixels with snow -> [ix]
                        swe_bin = mask_out[ind][ix]
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

                        else:
                            daily_outputs['swe_unavail'].loc[b,name] = 0
                            daily_outputs['swe_avail'].loc[b,name] = 0
                            daily_outputs['swe_z'].loc[b,name] = 0
                            daily_outputs['swe_vol'].loc[b,name] = 0

                    elif k == 'swi_z':
                        # Not masked out by pixels with snow
                        mask_out = np.multiply(self.outputs['swi_z'][iters],mask)
                        r = np.nansum(mask_out[ind])
                        daily_outputs['swi_z'].loc[b,name] = (
                                        np.multiply(np.nanmean(mask_out[ind]),
                                                    self.depth_factor) )
                        daily_outputs['swi_vol'].loc[b,name] = (
                                        np.multiply(r,self.conversion_factor) )

                    elif k in ['evap_z','depth']:
                        # masked out by pixels with snow -> [ix]
                        r = np.nanmean(mask_out[ind][ix])

                        if k == 'depth':
                            daily_outputs[k].loc[b,name] = (
                                        np.multiply(r,self.depth_factor*1000) )

                        else:
                            daily_outputs[k].loc[b,name] = (
                                            np.multiply(r,self.depth_factor) )

                    elif k in ['coldcont','density']:
                        # masked out by pixels with snow -> [ix]
                        daily_outputs[k].loc[b,name] = np.nanmean(mask_out[ind][ix])

                    elif k == 'precip_z':
                        # not masked out by pixels with snow
                        nanmask = mask == 0
                        mask_out = np.multiply(precip,mask)
                        mask_out[nanmask] = np.nan
                        r = np.nanmean(mask_out[ind])
                        rv = np.nansum(mask_out[ind])

                        daily_outputs[k].loc[b,name] = (
                                        np.multiply(r,self.depth_factor) )
                        daily_outputs['precip_vol'].loc[b,name] = (
                                        np.multiply(rv,self.conversion_factor) )

                    elif k == 'rain_z':
                        # not masked out by pixels with snow
                        nanmask = mask == 0
                        mask_out = np.multiply(rain,mask)
                        mask_out[nanmask] = np.nan

                        r = np.nanmean(mask_out[ind])
                        daily_outputs[k].loc[b,name] = (
                                        np.multiply(r,self.depth_factor) )

                # Add basin total mean/volume field once all elevations are done
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
                    total = np.multiply(np.nanmean(out[ixs]), self.depth_factor*1000)
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
                    nanmask = mask == 0
                    out = np.multiply(precip,mask)
                    out[nanmask] = np.nan
                    total = np.multiply(np.nanmean(out), self.depth_factor)
                    daily_outputs[k].loc['total',name] = copy.deepcopy(total)

                    s = np.nansum(daily_outputs['precip_vol'][name].values)
                    daily_outputs['precip_vol'].loc['total',name] = copy.deepcopy(s)

                if k == 'rain_z':
                    nanmask = mask == 0
                    out = np.multiply(rain,mask)
                    out[nanmask] = np.nan

                    total = np.multiply(np.nanmean(out), self.depth_factor)
                    daily_outputs[k].loc['total',name] = copy.deepcopy(total)

            # Send daily results to database
            # if self.location == 'database':
            df = daily_outputs[k].copy()
            df = df.round(decimals = 3)
            database.package_results.package(self, df, k,
                                             self.outputs['dates'][iters])

        # It's nice to see where we are...
        debug = 'snow file: %s, hours: %s, date: %s'%(
                                self.rundirs_dict[hr],
                                str(hr - t),
                                out_date.date().strftime("%Y-%-m-%-d"))

        self._logger.debug(debug)

        # Step this along so that we can see how many hours between outputs
        t = int(hr)

    # Append date to report name
    parts = self.report_name.split('.')
    self.report_name = ( parts[0]
                       + self.outputs['dates'][-1].date().strftime("%Y%m%d")
                       + '.' + parts[1] )
