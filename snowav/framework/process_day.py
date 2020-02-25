
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

def process(nc_path, topo_path, value, log):
    '''
    Process single day snow.nc files from command line call.

    Currently only processing swe_z and swe_vol values.

    Args
    ------
    nc_path : str
    topo_path : str
    value : str
    log : object

    '''

    log.info('Using topo {}...'.format(topo_path))

    out = masks(topo_path, False, plotorder = None, plotlabels = None)
    dem = out['dem']*3.28
    bmask = out['masks']
    plotorder = out['plotorder']
    labels = out['labels']
    depth_factor = 0.03937
    depthlbl = 'in'
    vollbl = 'TAF'
    elevlbl = 'ft'
    snowbands = [0,1,2]
    embands = [6,7,8,9]
    cclimit = -5*1000*1000
    show = False
    edges = np.arange(4000,15000,1000)
    def_edges = np.arange(3000,14000,1000)
    ixd = np.digitize(dem,edges)

    if type(nc_path) != list:
        nc_path = [nc_path]
    else:
        nc_path = nc_path

    topo = ts.get_topo_stats(topo_path)
    pixel = int(topo['dv'])
    conversion_factor = ((pixel**2) * 0.000000810713194*0.001)

    ncf = nc.Dataset(nc_path[0])
    t = nc.num2date(ncf.variables['time'][0],ncf.variables['time'].units)
    ncf.close()
    wy = handle_year_stradling(t) + 1

    vars = OrderedDict([('coldcont','cold content'),
                 ('swe_z','snow water equivalent depth'),
                 ('swe_vol','snow water equivalent volume'),
                 ('swe_avail','snow water equivalent available for melt'),
                 ('swe_unavail','snow water equivalent unavailable for melt')])

    results = {}
    results['df'] = {}
    outputs = {'swi_z':[],'evap_z':[],'snowmelt':[],'swe_z':[],'depth':[],
               'dates':[],'time':[],'density':[],'coldcont':[], 'path':[]}

    rundirs_dict = {}
    for ncp in nc_path:
        log.info('Loading {}...'.format(ncp))

        path = ncp.split('snow.nc')[0]
        output = iSnobalReader(path, snowbands = snowbands, embands = embands,
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
        outputs['path'] = ncp

    dz = pd.DataFrame(0, index = def_edges, columns = bmask.keys())

    if len(outputs['dates']) > 1:
        if outputs['dates'][0] == outputs['dates'][1]:
            outputs['dates'][1] = outputs['dates'][1] + timedelta(hours=1/60)

    log.info('Generating results for {}...'.format(value))

    for iters, out_date in enumerate(outputs['dates']):

        # Initialize output dataframes for every day
        daily_outputs = {'swe_vol':dz.copy(),
                         'swe_avail':dz.copy(),
                         'swe_unavail':dz.copy(),
                         'swe_z':dz.copy(),
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
                elevbin = ixd*mask

                # If the key is something we've already calculated (i.e., not
                # volume), get the masked subbasin output, otherwise it will be
                # calculated below
                if k in outputs.keys():
                    mask_out = outputs[k][iters]*mask

                # make a subbasin, by elevation, snow-free mask
                swe_mask_sub = copy.deepcopy(swe)*mask
                cold_sub = cold*mask

                # Now mask the subbasins by elevation band
                for n in np.arange(0,len(def_edges)):
                    ind = elevbin == n
                    b = def_edges[n]

                    # index for pixels with SWE, in subbasin, in elevation band,
                    # needed for mean values
                    ix = swe_mask_sub[ind] > 0

                    # Assign values
                    if k == 'swe_z':
                        mask_out = outputs['swe_z'][iters]*mask
                        ccb = cold_sub[ind]
                        cind = ccb > cclimit

                        # masked out by pixels with snow -> [ix]
                        swe_bin = mask_out[ind]

                        if swe_bin.size:
                            r = np.nansum(swe_bin[cind])
                            ru = np.nansum(swe_bin[~cind])

                            daily_outputs['swe_unavail'].loc[b,name] = ru*conversion_factor
                            daily_outputs['swe_avail'].loc[b,name] = r*conversion_factor
                            daily_outputs['swe_z'].loc[b,name] = np.nanmean(swe_bin)*depth_factor
                            daily_outputs['swe_vol'].loc[b,name] = np.nansum(swe_bin)*conversion_factor

                        else:
                            daily_outputs['swe_unavail'].loc[b,name] = np.nan
                            daily_outputs['swe_avail'].loc[b,name] = np.nan
                            daily_outputs['swe_z'].loc[b,name] = np.nan
                            daily_outputs['swe_vol'].loc[b,name] = np.nan


                # Add basin total mean/volume field once all elevations are done
                isub = swe_mask_sub > 0

                if k == 'swe_z':
                    # calc swe_z derived values
                    out = outputs[k][iters]*mask
                    total = np.nanmean(out)*depth_factor

                    s = np.nansum(daily_outputs['swe_vol'][name].values)
                    s1 = np.nansum(daily_outputs['swe_avail'][name].values)
                    s2 = np.nansum(daily_outputs['swe_unavail'][name].values)


            if k == 'swe_vol':
                df = daily_outputs[k].copy()
                df = df.round(decimals = 3)

        results['df'][nc_path[iters]] = df.fillna(0)

    results['outputs'] = outputs
    results['masks'] = bmask
    results['plotorder'] = plotorder
    results['edges'] = def_edges
    results['labels'] = labels

    return results
