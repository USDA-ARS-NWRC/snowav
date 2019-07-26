
import numpy as np
import os
import copy
import pandas as pd
from datetime import datetime, timedelta
import logging
import coloredlogs
import netCDF4 as nc
from snowav.database.package_results import package
from snowav.database.database import delete
from snowav.utils.utilities import precip
import warnings

def process(args):
    '''
    Process AWSM model results made in awsm_daily format, remove any existing
    database results with the same run_name and date, and place results on
    database. See README.md and CoreConfig.ini for more information.

    Args
    ------
    args : dict

    includes the following fields:
        outputs : dict
            from outputs()
        ixs : int
            starting index
        ixe : int
            ending index
        rundirs_dict : dict
        edges : array
        masks : dict
            from masks()
        nrows : int
        ncols : int
        plotorder : list
        connector : str
        basins : dict
        run_name : str
        conversion_factor : float
        depth_factor : float
        wy : int
        vars : dict
            from parse()
        ixd : array
            from parse()
        run_id : int
        vollbl : str
        depthlbl : str
        vid : int
            variable id
        connector : str

    Returns
    -------
    precip_total : array
    rain_total : array
    log : list
    density : array
        currently the last day of processing
    flag : bool

    '''

    pd.options.mode.chained_assignment = None
    warnings.filterwarnings('ignore')

    log = []
    cclimit = -5*1000*1000
    t = 0
    outputs = args['outputs']
    ixs = args['ixs']
    ixe = args['ixe']
    rundirs_dict = args['rundirs_dict']
    edges = args['edges']
    masks = args['masks']
    ixd = args['ixd']
    plotorder = args['plotorder']
    run_name = args['run_name']
    fc = args['conversion_factor']
    fd = args['depth_factor']
    lbls = {}
    lbls['depthlbl'] = args['depthlbl']
    lbls['vollbl'] = args['vollbl']
    print('process, ', masks.keys())
    dz = pd.DataFrame(0, index = edges, columns = masks.keys())
    precip_total = np.zeros((args['nrows'],args['ncols']))
    rain_total = np.zeros((args['nrows'],args['ncols']))

    for iters, out_date in enumerate(outputs['dates']):

        # delete anything on the database with the same run_name and date
        for bid in plotorder:
            out = delete(args['connector'], args['basins'], out_date, out_date,
                         bid, args['run_name'])
            log.append(out[0])

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

        density = {}
        for name in masks:
            density[name] = {}

        # Get daily rain from hourly input data
        flag, path, pre, rain = precip(rundirs_dict[outputs['time'][iters]],
                                       '/smrfOutputs/precip.nc')

        if flag:
            precip_total = precip_total + pre
            rain_total = rain_total + rain

        else:
            log.append(' WARNING! Expected to find {} but it is not a valid '
                       'file, precip will not be calculated or put on the '
                       'database, no precip figures will be made'.format(path))

        # Get a snow-free mask ready
        swe = copy.deepcopy(outputs['swe_z'][iters])
        cold = copy.deepcopy(outputs['coldcont'][iters])

        # Loop over outputs (depths are copied, volumes are calculated)
        for k in args['vars'].keys():

            # Mask by subbasin
            for name in masks:
                mask = copy.deepcopy(masks[name]['mask'])
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
                for n in np.arange(0,len(edges)):
                    ind = elevbin == n
                    b = edges[n]

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
                            daily_outputs['swe_unavail'].loc[b,name] = ru*fc
                            daily_outputs['swe_avail'].loc[b,name] = r*fc
                            daily_outputs['swe_z'].loc[b,name] = np.nanmean(swe_bin)*fd
                            daily_outputs['swe_vol'].loc[b,name] = np.nansum(swe_bin)*fc

                        else:
                            daily_outputs['swe_unavail'].loc[b,name] = np.nan
                            daily_outputs['swe_avail'].loc[b,name] = np.nan
                            daily_outputs['swe_z'].loc[b,name] = np.nan
                            daily_outputs['swe_vol'].loc[b,name] = np.nan

                    elif k == 'swi_z':
                        # Not masked out by pixels with snow
                        mask_out = np.multiply(outputs['swi_z'][iters],mask)
                        r = np.nansum(mask_out[ind])
                        daily_outputs['swi_z'].loc[b,name] = np.nanmean(mask_out[ind])*fd
                        daily_outputs['swi_vol'].loc[b,name] = r*fc

                    elif k in ['evap_z','depth']:
                        # masked out by pixels with snow -> [ix]
                        r = np.nanmean(mask_out[ind])

                        if k == 'depth':
                            daily_outputs[k].loc[b,name] = r*fd*1000

                        else:
                            daily_outputs[k].loc[b,name] = r*fd

                    elif k in ['coldcont','density']:
                        # masked out by pixels with snow -> [ix]
                        daily_outputs[k].loc[b,name] = np.nanmean(mask_out[ind][ix])
                        density[name][edges[n]] = mask_out[ind][ix]

                    elif k == 'precip_z' and flag:
                        # not masked out by pixels with snow
                        nanmask = mask == 0
                        mask_out = pre*mask
                        mask_out[nanmask] = np.nan
                        r = np.nanmean(mask_out[ind])
                        rv = np.nansum(mask_out[ind])

                        daily_outputs[k].loc[b,name] = r*fd
                        daily_outputs['precip_vol'].loc[b,name] = rv*fc

                    elif k == 'rain_z' and flag:
                        # not masked out by pixels with snow
                        nanmask = mask == 0
                        mask_out = rain*mask
                        mask_out[nanmask] = np.nan

                        r = np.nanmean(mask_out[ind])
                        daily_outputs[k].loc[b,name] = r*fd

                # Add basin total mean/volume field once all elevations are done
                isub = swe_mask_sub > 0

                if k in ['evap_z','coldcont']:
                    # Mask by snow-free
                    out = outputs[k][iters]*mask
                    total = np.nanmean(out)*fd
                    daily_outputs[k].loc['total',name] = copy.deepcopy(total)

                if k == 'swe_z':
                    # calc swe_z derived values
                    out = np.multiply(outputs[k][iters],mask)
                    total = np.multiply(np.nanmean(out), fd)

                    s = np.nansum(daily_outputs['swe_vol'][name].values)
                    s1 = np.nansum(daily_outputs['swe_avail'][name].values)
                    s2 = np.nansum(daily_outputs['swe_unavail'][name].values)
                    daily_outputs['swe_vol'].loc['total',name] = copy.deepcopy(s)
                    daily_outputs['swe_avail'].loc['total',name] = copy.deepcopy(s1)
                    daily_outputs['swe_unavail'].loc['total',name] = copy.deepcopy(s2)

                    daily_outputs[k].loc['total',name] = copy.deepcopy(total)

                if k == 'depth':
                    # Mask by snow-free
                    out = outputs[k][iters]*mask
                    total = np.nanmean(out)*fd*1000
                    daily_outputs[k].loc['total',name] = copy.deepcopy(total)

                if k == 'density':
                    # Mask by snow-free
                    out = outputs[k][iters]*mask
                    total = np.nanmean(out[isub])
                    daily_outputs[k].loc['total',name] = copy.deepcopy(total)

                if k == 'swi_z':
                    # Do not mask by snow free
                    out = outputs[k][iters]*mask
                    total = np.nanmean(out)*fd
                    daily_outputs[k].loc['total',name] = copy.deepcopy(total)

                    s = np.nansum(daily_outputs['swi_vol'][name].values)
                    daily_outputs['swi_vol'].loc['total',name] = copy.deepcopy(s)

                if k == 'precip_z' and flag:
                    nanmask = mask == 0
                    out = pre*mask
                    out[nanmask] = np.nan
                    total = np.nanmean(out)*fd
                    daily_outputs[k].loc['total',name] = copy.deepcopy(total)

                    s = np.nansum(daily_outputs['precip_vol'][name].values)
                    daily_outputs['precip_vol'].loc['total',name] = copy.deepcopy(s)

                if k == 'rain_z' and flag:
                    nanmask = mask == 0
                    out = rain*mask
                    out[nanmask] = np.nan

                    total = np.nanmean(out)*fd
                    daily_outputs[k].loc['total',name] = copy.deepcopy(total)

            # Send daily results to database
            df = daily_outputs[k].copy()
            df = df.round(decimals = 3)
            package(args['connector'], lbls, args['basins'], df, args['run_id'],
                    args['vid'], k, out_date, run_name)

        hr = int(outputs['time'][iters])
        log.append(' processed: {}, elapsed hours: {}, date: {}'.format(
                   rundirs_dict[outputs['time'][iters]],
                   str(hr - t), out_date.date().strftime("%Y-%-m-%-d")))

        # hours between outputs
        t = int(hr)

    return flag, log, precip_total, rain_total, density
