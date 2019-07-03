
from snowav.utils.OutputReader import iSnobalReader
import numpy as np
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
from snowav.utils.utilities import masks
import snowav.utils.get_topo_stats as ts
from snowav.utils.wyhr import handle_year_stradling, calculate_date_from_wyhr

def process(nc_path, topo_path, value):
    '''

    This version of process() is for the command line snowav_day, and
    processes either a single snow.nc file, or two snow.nc files to
    calculate a simple SWE volume difference. This could be expanded
    to include other output variables.

    '''

    print('Using topo {}...'.format(topo_path))

    out = masks(topo_path, plotorder = None, plotlabels = None)
    dem = out['dem']*3.28
    bmask = out['masks']
    nrows = out['nrows']
    ncols = out['ncols']
    plotorder = out['plotorder']
    labels = out['labels']
    depth_factor = 0.03937
    depthlbl = 'in'
    vollbl = 'TAF'
    elevlbl = 'ft'
    filetype = 'netcdf'
    snowbands = [0,1,2]
    embands = [6,7,8,9]
    cclimit = -5*1000*1000
    show = False
    name_append = ''
    edges = np.arange(4000,14000,1000)

    if type(nc_path) != list:
        nc_path = [nc_path]
    else:
        nc_path = nc_path

    topo = ts.get_topo_stats(topo_path, filetype = 'netcdf')
    pixel = int(topo['dv'])

    ncf = nc.Dataset(nc_path[0])
    t = nc.num2date(ncf.variables['time'][0],ncf.variables['time'].units)
    ncf.close()
    wy = handle_year_stradling(t) + 1

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

    results = {}

    conversion_factor = ((pixel**2) * 0.000000810713194*0.001)
    ixd = np.digitize(dem,edges)

    outputs = {'swi_z':[],'evap_z':[],'snowmelt':[],'swe_z':[],'depth':[],
               'dates':[],'time':[],'density':[],'coldcont':[]}

    rundirs_dict = {}
    for ncp in nc_path:
        print('Loading {}...'.format(ncp))

        path = ncp.split('snow.nc')[0]
        output = iSnobalReader(path, filetype, snowbands = snowbands,
                               embands = embands, wy = wy)

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

    date = outputs['dates'][0]

    # Daily, by elevation
    dz = pd.DataFrame(0, index = edges, columns = bmask.keys())

    if len(outputs['dates']) > 1:
        if outputs['dates'][0] == outputs['dates'][1]:
            outputs['dates'][1] = outputs['dates'][1] + timedelta(hours=1/60)

    print('Generating results for {}...'.format(value))

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

        # Get a snow-free mask ready
        swe = copy.deepcopy(outputs['swe_z'][iters])
        cold = copy.deepcopy(outputs['coldcont'][iters])

        # Loop over outputs (depths are copied, volumes are calculated)
        for k in vars.keys():

            # Mask by subbasin
            for name in bmask:
                mask = copy.deepcopy(bmask[name]['mask'])
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
                                                    depth_factor) )
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
                                                    depth_factor) )
                        daily_outputs['swi_vol'].loc[b,name] = (
                                        np.multiply(r,conversion_factor) )

                    elif k in ['evap_z','depth']:
                        # masked out by pixels with snow -> [ix]
                        r = np.nanmean(mask_out[ind])

                        if k == 'depth':
                            daily_outputs[k].loc[b,name] = (
                                        np.multiply(r,depth_factor*1000) )

                        else:
                            daily_outputs[k].loc[b,name] = (
                                            np.multiply(r,depth_factor) )

                    elif k in ['coldcont','density']:
                        # masked out by pixels with snow -> [ix]
                        daily_outputs[k].loc[b,name] = np.nanmean(mask_out[ind][ix])

                # Add basin total mean/volume field once all elevations are done
                isub = swe_mask_sub > 0

                if k == 'swe_z':
                    # calc swe_z derived values
                    out = np.multiply(outputs[k][iters],mask)
                    total = np.multiply(np.nanmean(out), depth_factor)

                    s = np.nansum(daily_outputs['swe_vol'][name].values)
                    s1 = np.nansum(daily_outputs['swe_avail'][name].values)
                    s2 = np.nansum(daily_outputs['swe_unavail'][name].values)
                    # daily_outputs['swe_vol'].loc['total',name] = copy.deepcopy(s)
                    # daily_outputs['swe_avail'].loc['total',name] = copy.deepcopy(s1)
                    # daily_outputs['swe_unavail'].loc['total',name] = copy.deepcopy(s2)

                    # daily_outputs[k].loc['total',name] = copy.deepcopy(total)

                if k == 'depth':
                    # Mask by snow-free
                    out = np.multiply(outputs[k][iters],mask)
                    total = np.multiply(np.nanmean(out), depth_factor*1000)
                    daily_outputs[k].loc['total',name] = copy.deepcopy(total)

                if k == 'density':
                    # Mask by snow-free
                    out = np.multiply(outputs[k][iters],mask)
                    total = np.nanmean(out[isub])
                    daily_outputs[k].loc['total',name] = copy.deepcopy(total)

                if k == 'swi_z':
                    # Do not mask by snow free
                    out = np.multiply(outputs[k][iters],mask)
                    total = np.multiply(np.nanmean(out), depth_factor)
                    daily_outputs[k].loc['total',name] = copy.deepcopy(total)

                    s = np.nansum(daily_outputs['swi_vol'][name].values)
                    daily_outputs['swi_vol'].loc['total',name] = copy.deepcopy(s)

            df = daily_outputs[k].copy()
            df = df.round(decimals = 3)

    # results[out_date] = daily_outputs[value]
    results['df'] = df.fillna(0)
    results['outputs'] = outputs
    results['masks'] = bmask
    results['plotorder'] = plotorder
    # results['lims'] = lims
    results['edges'] = edges
    results['labels'] = labels

    return results
