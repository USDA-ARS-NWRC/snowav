
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
from datetime import timedelta

def process(self):
    '''
    This function processes and summarizes iSnobal outputs.

    It loops over each day in the specified directory and:
        - sums total precip and rain from the ipw ppt ppt files
        - loops over each daily_outputs ('swe_z', 'swe_vol', etc)
        - loops over sub-basins
        - loops over elevation band, and calculates volumes and mean depths
            for each daily_outputs value in that elevation band
        - once all elevation bands for each subbasin have been calculated,
            'total' fields (either mean depth for the subbasin or total volume
            for the subbasin) are calculated
        - finally, that is sent to the database

    '''

    # Suppress warnings - empty slices
    pd.options.mode.chained_assignment = None
    warnings.filterwarnings('ignore')

    # based on an average of 60 W/m2 from TL paper
    cclimit = -5*1000*1000

    # Daily, by elevation
    dz = pd.DataFrame(0, index = self.edges, columns = self.masks.keys())

    adj = 0
    t = 0
    self.precip_total = np.zeros((self.nrows,self.ncols))
    self.rain_total = np.zeros((self.nrows,self.ncols))

    for iters, out_date in enumerate(self.outputs['dates']):

        # end processing at self.ixe
        if iters < self.ixs:
            continue

        if iters > self.ixe:
            break

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
        ppt_path = sf.split('output')[0] + '/smrfOutputs/precip.nc'
        percent_snow_path = ppt_path.replace('precip','percent_snow')

        precip = np.zeros((self.nrows,self.ncols))
        rain = np.zeros((self.nrows,self.ncols))

        if os.path.isfile(ppt_path):

            ppt = nc.Dataset(ppt_path, 'r')
            percent_snow = nc.Dataset(percent_snow_path, 'r')

            # For the wy2019 daily runs, precip.nc always has an extra hour...
            if len(ppt.variables['time'][:]) <= 24:
                nb = 24
            else:
                nb = len(ppt.variables['time'][:])


            for nb in range(0,nb):
                pre = ppt['precip'][nb]
                ps = percent_snow['percent_snow'][nb]

                rain = rain + np.multiply(pre,(1-ps))
                precip = precip + copy.deepcopy(pre)

            ppt.close()
            percent_snow.close()

        # self.precip_total = np.nansum(np.dstack((self.precip_total,copy.deepcopy(precip))),2)
        # self.rain_total = np.nansum(np.dstack((self.rain_total,copy.deepcopy(rain))),2)
        self.precip_total = self.precip_total + copy.deepcopy(precip)
        self.rain_total = self.rain_total + copy.deepcopy(rain)

        # Get a snow-free mask ready
        swe = copy.deepcopy(self.outputs['swe_z'][iters])
        cold = copy.deepcopy(self.outputs['coldcont'][iters])

        # Loop over outputs (depths are copied, volumes are calculated)
        for k in self.vars.keys():

            # Mask by subbasin
            for name in self.masks:
                mask = copy.deepcopy(self.masks[name]['mask'])
                mask = mask.astype(float)

                mask[mask < 1] = np.nan
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
                        ccb = cold_sub[ind]
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

                        else:
                            daily_outputs['swe_unavail'].loc[b,name] = np.nan
                            daily_outputs['swe_avail'].loc[b,name] = np.nan
                            daily_outputs['swe_z'].loc[b,name] = np.nan
                            daily_outputs['swe_vol'].loc[b,name] = np.nan

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
                        r = np.nanmean(mask_out[ind])

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

                if k in ['evap_z','coldcont']:
                    # Mask by snow-free
                    out = np.multiply(self.outputs[k][iters],mask)
                    total = np.multiply(np.nanmean(out), self.depth_factor)
                    daily_outputs[k].loc['total',name] = copy.deepcopy(total)

                if k == 'swe_z':
                    # calc swe_z derived values
                    out = np.multiply(self.outputs[k][iters],mask)
                    total = np.multiply(np.nanmean(out), self.depth_factor)

                    s = np.nansum(daily_outputs['swe_vol'][name].values)
                    s1 = np.nansum(daily_outputs['swe_avail'][name].values)
                    s2 = np.nansum(daily_outputs['swe_unavail'][name].values)
                    daily_outputs['swe_vol'].loc['total',name] = copy.deepcopy(s)
                    daily_outputs['swe_avail'].loc['total',name] = copy.deepcopy(s1)
                    daily_outputs['swe_unavail'].loc['total',name] = copy.deepcopy(s2)

                    daily_outputs[k].loc['total',name] = copy.deepcopy(total)

                if k == 'depth':
                    # Mask by snow-free
                    out = np.multiply(self.outputs[k][iters],mask)
                    total = np.multiply(np.nanmean(out), self.depth_factor*1000)
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
            df = daily_outputs[k].copy()
            df = df.round(decimals = 3)

            # If there are already results on database, only insert if overwriting
            if self.pflag is False:
                database.package_results.package(self, df, k,
                                                 self.outputs['dates'][iters])

            if (self.pflag is True) and (self.write_db is True):
                database.package_results.package(self, df, k,
                                                 self.outputs['dates'][iters])

        # Add water year totals for swi and evap
        if (out_date.date() != datetime.datetime(self.wy-1,10,1).date() and
            (self.pflag is False) ):
            database.package_results.post_process(self, out_date)

        if (out_date.date() != datetime.datetime(self.wy-1,10,1).date() and
            (self.pflag is True) and (self.write_db is True) ):
            database.package_results.post_process(self, out_date)

        # It's nice to see where we are...
        debug = 'snow file: %s, hours: %s, date: %s'%(
                                self.rundirs_dict[hr],
                                str(hr - t),
                                out_date.date().strftime("%Y-%-m-%-d"))

        self._logger.debug(debug)

        # Step this along so that we can see how many hours between outputs
        t = int(hr)

    # Save images
    swiz = np.zeros((self.nrows,self.ncols))
    for n in range(0,len(self.outputs['swi_z'])):
        swiz = swiz + self.outputs['swi_z'][n]

    # self.end_date = self.end_date + timedelta(hours=1)

    np.save('{}{}'.format(self.figs_path,'swi_z.npy'),swiz)
    np.save('{}{}'.format(self.figs_path,'precip.npy'),self.precip_total)
    np.save('{}{}'.format(self.figs_path,'rain.npy'),self.rain_total)
    np.save('{}{}'.format(self.figs_path,'swe_z.npy'),self.outputs['swe_z'][self.ixe])
