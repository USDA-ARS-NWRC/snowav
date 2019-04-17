
from snowav.utils.OutputReader import iSnobalReader
import numpy as np
from spatialnc import ipw
import os
import copy
import pandas as pd
from datetime import datetime, timedelta
import logging
import coloredlogs
import netCDF4 as nc
from snowav import database
import warnings
from collections import OrderedDict
from snowav.database.tables import Basins

class day(object):

    depth_factor = 0.03937
    step = 1000
    depthlbl = 'in'
    vollbl = 'TAF'
    elevlbl = 'ft'
    filetype = 'netcdf'
    snowbands = [0,1,2]
    embands = [6,7,8,9]

    def process_day(day, nc_path, value, basin=None, dem=None, wy=None):
        '''

        '''

        # Suppress warnings - empty slices
        pd.options.mode.chained_assignment = None
        warnings.filterwarnings('ignore')

        vars = OrderedDict([('coldcont','cold content'),
                     ('evap_z','evaporation depth'),
                     ('rain_z','rain depth'),
                     ('density','density'),
                     ('depth','depth'),
                     ('precip_z','precipitation depth'),
                     ('precip_vol','precipitation volume'),
                     ('swi_z','surface water input depth'),
                     ('swi_vol','surface water input volume'),
                     ('swe_z','snow water equivalent depth'),
                     ('swe_vol','snow water equivalent volume'),
                     ('swe_avail','snow water equivalent available for melt'),
                     ('swe_unavail','snow water equivalent unavailable for melt')])

        barcolors = ['xkcd:cobalt',
                     'xkcd:mustard green',
                     'xkcd:lichen',
                     'xkcd:pale green',
                     'xkcd:blue green',
                     'xkcd:bluish purple',
                     'xkcd:lightish purple',
                     'xkcd:deep magenta']

        plotorder = Basins.basins[basin]['defaults']['plotorder']
        dem_path = Basins.basins[basin]['defaults']['dem_path']
        edges = Basins.basins[basin]['defaults']['edges']
        results = {}

        # if wy is None:
        date_time = datetime.now()

        if date_time.month in (10,11,12):
            wy = date_time.year-1

        else:
            wy = date_time.year

        ncf = nc.Dataset(dem_path, 'r')
        dem = ncf.variables['dem'][:] * 3.28
        nrows = len(dem[:,0])
        ncols = len(dem[0,:])
        mask = ncf.variables['mask'][:]
        x = ncf.variables['x'][:]
        pixel = x[1] - x[0]
        conversion_factor = ((pixel**2) * 0.000000810713194*0.001)

        precip_total = np.zeros((nrows,ncols))
        rain_total = np.zeros((nrows,ncols))

        masks = dict()

        for lbl in plotorder:

            # Exceptions
            if lbl == 'Cherry Creek':
                nclbl = 'Cherry'
            else:
                nclbl = lbl

            if lbl != plotorder[0]:
                masks[lbl] = {'mask': ncf[nclbl + ' mask'][:],
                              'label': lbl}

            else:
                masks[lbl] = {'mask': ncf['mask'][:],
                              'label': nclbl}

        ncf.close()

        # min_dem = np.min(np.min(dem))
        # max_dem = np.max(np.max(dem))
        # min_dem -= min_dem % + 1000
        # max_dem -= max_dem % - 1000
        ixd = np.digitize(dem,edges)
        xlims = (0,len(edges))

        # pre-processing
        outputs = {'swi_z':[], 'evap_z':[], 'snowmelt':[], 'swe_z':[],
                        'depth':[], 'dates':[], 'time':[], 'density':[],
                        'coldcont':[] }

        # will need to look over this
        if type(nc_path) != list:
            nc_path = [nc_path]

        rundirs_dict = {}
        for ncp in nc_path:
            path = ncp.split('snow.nc')[0]
            output = iSnobalReader(path,
                                   day.filetype,
                                   snowbands = day.snowbands,
                                   embands = day.embands,
                                   wy = wy)

            for n in range(0,len(output.em_data[8])):
                outputs['swi_z'].append(output.em_data[8][n,:,:])
                outputs['snowmelt'].append(output.em_data[7][n,:,:])
                outputs['evap_z'].append(output.em_data[6][n,:,:])
                outputs['coldcont'].append(output.em_data[9][n,:,:])
                outputs['swe_z'].append(output.snow_data[2][n,:,:])
                outputs['depth'].append(output.snow_data[0][n,:,:])
                outputs['density'].append(output.snow_data[1][n,:,:])

            outputs['dates'] = np.append(outputs['dates'],output.dates)
            outputs['time'] = np.append(outputs['time'],output.time)

            for t in output.time:
                rundirs_dict[int(t)] = path

        # based on an average of 60 W/m2 from TL paper
        cclimit = -5*1000*1000

        # Daily, by elevation
        dz = pd.DataFrame(0, index = edges, columns = masks.keys())

        for iters, out_date in enumerate(outputs['dates']):
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

            # Get daily rain from hourly input data
            hr = int(outputs['time'][iters])
            sf = rundirs_dict[hr].replace('runs','data')
            sf = sf.replace('run','data')
            # ppt_path = sf.split('output')[0] + '/smrfOutputs/precip.nc'
            ppt_path = sf + '/smrfOutputs/precip.nc'
            percent_snow_path = ppt_path.replace('precip','percent_snow')

            precip = np.zeros((nrows,ncols))
            rain = np.zeros((nrows,ncols))

            if os.path.isfile(ppt_path):

                ppt = nc.Dataset(ppt_path, 'r')
                percent_snow = nc.Dataset(percent_snow_path, 'r')

                # For the wy2019 daily runs, precip.nc always has an extra hour...
                if len(ppt.variables['precip'][:]) > 24:
                    nb = 24
                else:
                    nb = len(ppt.variables['precip'][:])

                for nb in range(0,nb):
                    pre = ppt['precip'][nb]
                    ps = percent_snow['percent_snow'][nb]

                    rain = rain + np.multiply(pre,(1-ps))
                    precip = precip + copy.deepcopy(pre)

                ppt.close()
                percent_snow.close()

            # self.precip_total = np.nansum(np.dstack((self.precip_total,copy.deepcopy(precip))),2)
            # self.rain_total = np.nansum(np.dstack((self.rain_total,copy.deepcopy(rain))),2)
            precip_total = precip_total + copy.deepcopy(precip)
            rain_total = rain_total + copy.deepcopy(rain)

            # Get a snow-free mask ready
            swe = copy.deepcopy(outputs['swe_z'][iters])
            cold = copy.deepcopy(outputs['coldcont'][iters])

            # Loop over outputs (depths are copied, volumes are calculated)
            for k in vars.keys():

                # Mask by subbasin
                for name in masks:
                    mask = copy.deepcopy(masks[name]['mask'])
                    mask = mask.astype(float)

                    mask[mask < 1] = np.nan
                    elevbin = np.multiply(ixd,mask)

                    # If the key is something we've already calculated (i.e., not
                    # volume), get the masked subbasin output, otherwise it will be
                    # calculated below
                    if k in outputs.keys():
                        mask_out = np.multiply(outputs[k][iters],mask)

                    # make a subbasin, by elevation, snow-free mask
                    swe_mask_sub = np.multiply(copy.deepcopy(swe),mask)
                    cold_sub = np.multiply(cold,mask)

                    # Now mask the subbasins by elevation band
                    for n in np.arange(0,len(edges)):
                        ind = elevbin == n
                        b = edges[n]

                        # index for pixels with SWE, in subbasin, in elevation band,
                        # needed for mean values
                        ix = swe_mask_sub[ind] > 0

                        # Assign values
                        if k == 'swe_z':
                            mask_out = np.multiply(outputs['swe_z'][iters],mask)
                            ccb = cold_sub[ind]
                            cind = ccb > cclimit

                            # masked out by pixels with snow -> [ix]
                            swe_bin = mask_out[ind]
                            if swe_bin.size:
                                r = np.nansum(swe_bin[cind])
                                ru = np.nansum(swe_bin[~cind])

                                daily_outputs['swe_unavail'].loc[b,name] = (
                                            np.multiply(ru,conversion_factor) )
                                daily_outputs['swe_avail'].loc[b,name] = (
                                            np.multiply(r,conversion_factor) )
                                daily_outputs['swe_z'].loc[b,name] = (
                                            np.multiply(np.nanmean(swe_bin),
                                                        day.depth_factor) )
                                daily_outputs['swe_vol'].loc[b,name] = (
                                            np.multiply(np.nansum(swe_bin),
                                                        conversion_factor) )

                            else:
                                daily_outputs['swe_unavail'].loc[b,name] = np.nan
                                daily_outputs['swe_avail'].loc[b,name] = np.nan
                                daily_outputs['swe_z'].loc[b,name] = np.nan
                                daily_outputs['swe_vol'].loc[b,name] = np.nan

                        elif k == 'swi_z':
                            # Not masked out by pixels with snow
                            mask_out = np.multiply(outputs['swi_z'][iters],mask)
                            r = np.nansum(mask_out[ind])
                            daily_outputs['swi_z'].loc[b,name] = (
                                            np.multiply(np.nanmean(mask_out[ind]),
                                                        day.depth_factor) )
                            daily_outputs['swi_vol'].loc[b,name] = (
                                            np.multiply(r,conversion_factor) )

                        elif k in ['evap_z','depth']:
                            # masked out by pixels with snow -> [ix]
                            r = np.nanmean(mask_out[ind])

                            if k == 'depth':
                                daily_outputs[k].loc[b,name] = (
                                            np.multiply(r,day.depth_factor*1000) )

                            else:
                                daily_outputs[k].loc[b,name] = (
                                                np.multiply(r,day.depth_factor) )

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
                                            np.multiply(r,day.depth_factor) )
                            daily_outputs['precip_vol'].loc[b,name] = (
                                            np.multiply(rv,conversion_factor) )

                        elif k == 'rain_z':
                            # not masked out by pixels with snow
                            nanmask = mask == 0
                            mask_out = np.multiply(rain,mask)
                            mask_out[nanmask] = np.nan

                            r = np.nanmean(mask_out[ind])
                            daily_outputs[k].loc[b,name] = (
                                            np.multiply(r,day.depth_factor) )

                    # Add basin total mean/volume field once all elevations are done
                    isub = swe_mask_sub > 0

                    if k in ['evap_z','coldcont']:
                        # Mask by snow-free
                        out = np.multiply(outputs[k][iters],mask)
                        total = np.multiply(np.nanmean(out), day.depth_factor)
                        daily_outputs[k].loc['total',name] = copy.deepcopy(total)

                    if k == 'swe_z':
                        # calc swe_z derived values
                        out = np.multiply(outputs[k][iters],mask)
                        total = np.multiply(np.nanmean(out), day.depth_factor)

                        s = np.nansum(daily_outputs['swe_vol'][name].values)
                        s1 = np.nansum(daily_outputs['swe_avail'][name].values)
                        s2 = np.nansum(daily_outputs['swe_unavail'][name].values)
                        daily_outputs['swe_vol'].loc['total',name] = copy.deepcopy(s)
                        daily_outputs['swe_avail'].loc['total',name] = copy.deepcopy(s1)
                        daily_outputs['swe_unavail'].loc['total',name] = copy.deepcopy(s2)

                        daily_outputs[k].loc['total',name] = copy.deepcopy(total)

                    if k == 'depth':
                        # Mask by snow-free
                        out = np.multiply(outputs[k][iters],mask)
                        total = np.multiply(np.nanmean(out), day.depth_factor*1000)
                        daily_outputs[k].loc['total',name] = copy.deepcopy(total)

                    if k == 'density':
                        # Mask by snow-free
                        out = np.multiply(outputs[k][iters],mask)
                        total = np.nanmean(out[isub])
                        daily_outputs[k].loc['total',name] = copy.deepcopy(total)

                    if k == 'swi_z':
                        # Do not mask by snow free
                        out = np.multiply(outputs[k][iters],mask)
                        total = np.multiply(np.nanmean(out), day.depth_factor)
                        daily_outputs[k].loc['total',name] = copy.deepcopy(total)

                        s = np.nansum(daily_outputs['swi_vol'][name].values)
                        daily_outputs['swi_vol'].loc['total',name] = copy.deepcopy(s)

                    if k == 'precip_z':
                        nanmask = mask == 0
                        out = np.multiply(precip,mask)
                        out[nanmask] = np.nan
                        total = np.multiply(np.nanmean(out), day.depth_factor)
                        daily_outputs[k].loc['total',name] = copy.deepcopy(total)

                        s = np.nansum(daily_outputs['precip_vol'][name].values)
                        daily_outputs['precip_vol'].loc['total',name] = copy.deepcopy(s)

                    if k == 'rain_z':
                        nanmask = mask == 0
                        out = np.multiply(rain,mask)
                        out[nanmask] = np.nan

                        total = np.multiply(np.nanmean(out), day.depth_factor)
                        daily_outputs[k].loc['total',name] = copy.deepcopy(total)

                # Send daily results to database
                df = daily_outputs[k].copy()
                df = df.round(decimals = 3)

            results[out_date] = daily_outputs[value]

        day.basin = basin
        day.value = value
        day.results = results
        day.outputs = outputs
        day.daily_outputs = daily_outputs
        day.path = path
        day.nc_path = nc_path
        day.masks = masks
        day.plotorder = plotorder
        day.barcolors = barcolors
        day.edges = edges
        day.dplcs = 1
        day.xlims = xlims
        day.name_append = ''

        return day
