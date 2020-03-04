
import coloredlogs
import copy
from datetime import datetime, timedelta
import logging
import netCDF4 as nc
import numpy as np
import pandas as pd
import os
from sys import exit
import warnings

from snowav.database.package_results import package
from snowav.database.database import delete, query
from snowav.utils.utilities import calculate, sum_precip, snow_line, \
    input_summary
from tablizer.defaults import Units
from tablizer.tablizer import get_existing_records

class Process(object):
    """ Process AWSM model results made in awsm_daily format.

    Args
    ------
    cfg: config object

    """

    def __init__(self, cfg):

        pd.options.mode.chained_assignment = None
        warnings.filterwarnings('ignore')

        elapsed_hours = 0
        precip_path = '/smrfOutputs/precip.nc'
        variables = cfg.variables.variables
        lbls = {'depthlbl': cfg.depthlbl,
                'vollbl': cfg.vollbl,
                'elevlbl': cfg.elevlbl}

        if os.path.splitext(cfg.connector)[1] == '.db':
            db = 'sqlite'
        else:
            db = 'sql'

        # images for figures
        precip_total = np.zeros((cfg.nrows,cfg.ncols))
        rain_total = np.zeros((cfg.nrows,cfg.ncols))

        # Check that topo and outputs are the same dimensions
        if cfg.outputs['swe_z'][0].shape != cfg.dem.shape:
            logging.error(" topo.nc dimensions: {} do not match output "
                "dimensions: {}".format(cfg.dem.shape,
                                        cfg.outputs['swe_z'][0].shape))
            exit()

        # process each date
        for iters, out_date in enumerate(cfg.outputs['dates']):
            wy_hour = int(cfg.outputs['time'][iters])

            # assume to start that we are processing everything
            proc_list = copy.deepcopy(cfg.variables.variables)
            pass_flag = False

            # delete anything on the database with the same run_name and date
            for bid in cfg.plotorder:
                results = query(cfg.connector, out_date, out_date,
                                cfg.run_name, cfg.basins, bid, 'swe_vol')

                if not results.empty:
                    if not cfg.db_overwrite:
                        # If database records exist we can skip all outputs
                        # but density on the last day, if making the figure
                        pass_flag = True
                        flag = False

                        if bid == cfg.plotorder[0]:
                            logging.info(' Skipping outputs for {}, {}, '
                                         'database records exist...'.format(bid,
                                         out_date.strftime("%Y-%-m-%-d %H:00")))
                        else:
                            logging.debug(' Skipping outputs for {}, {}, '
                                         'database records exist...'.format(bid,
                                         out_date.strftime("%Y-%-m-%-d %H:00")))

                        if (out_date == cfg.outputs['dates'][-1] and
                        cfg.density_flag):
                            proc_list = ['density']
                        else:
                            proc_list = []

                    else:
                        out = delete(cfg.connector, cfg.basins, out_date,
                                     out_date, bid, cfg.run_name)
                        if out != []:
                            if bid == cfg.plotorder[0]:
                                logging.info(out[0])
                            else:
                                logging.debug(out[0])

            density = {}
            for name in cfg.masks:
                density[name] = {}

            # Get daily rain from hourly input data
            # Only run this processing if database records don't exist, or if
            # we are making the precip figure
            # if not pass_flag or process.precip_depth_flag:
            #     logging.debug(' Processing precip for {}'.format(
            #                  out_date.strftime("%Y-%-m-%-d %H:00")))
            #     flag, path, pre, rain = precip(cfg.rundirs_dict[wy_hour],
            #                                    precip_path)
            #
            #     if flag:
            #         precip_total = precip_total + pre
            #         rain_total = rain_total + rain
            #     else:
            #         logging.warning(' Expected to find {} but it is not a valid '
            #                    'file, precip will not be calculated or put on '
            #                    'the database, no precip figures will be '
            #                    'made'.format(path))



            if cfg.inputs_flag:
                mask_list = []
                for name in cfg.masks:
                    mask_list.append(copy.deepcopy(cfg.masks[name]['mask']))

                sf = cfg.rundirs_dict[wy_hour].replace('runs','data')
                sf = sf.replace('run','data') + '/smrfOutputs/'
                inputs = os.listdir(sf)

                df = get_existing_records(cfg.connector, db)
                df = df.set_index('date_time')
                df.sort_index(inplace=True)

                for i, input in enumerate(inputs):
                    input_path = os.path.join(sf,input)
                    variable = os.path.splitext(input)[0]

                    if variable in cfg.inputs_variables:
                        for basin in cfg.inputs_basins:

                            basin_id = int(cfg.basins[basin]['basin_id'])

                            dfs = df[((df['variable'] == variable) &
                                     (df['run_name'] == cfg.run_name) &
                                     (df['basin_id'] == basin_id))]
                            idx = dfs.index

                            if (not idx.contains(out_date)
                            or (idx.contains(out_date) and cfg.db_overwrite)):
                                if basin == cfg.inputs_basins[0] and i == 0:
                                    logging.info(' Processing {} inputs, '
                                                 '{} '.format(basin,
                                         out_date.strftime("%Y-%-m-%-d %H:00")))
                                else:
                                    logging.debug(' Processing {}, '
                                         '{}, {} '.format(variable, basin,
                                         out_date.strftime("%Y-%-m-%-d %H:00")))

                                input_summary(input_path,
                                              variable,
                                              cfg.inputs_methods,
                                              cfg.inputs_percentiles,
                                              db,
                                              cfg.connector,
                                              cfg.run_name,
                                              basin_id,
                                              cfg.run_id,
                                              cfg.masks[basin]['mask'])

                                # handle precip different because we need
                                # summed images for figures
                                if variable == 'precip':
                                    ps_path = input_path.replace('precip.nc',
                                                            'percent_snow.nc')
                                    sum_precip(input_path,
                                               ps_path,
                                               precip_total,
                                               rain_total)

                            if idx.contains(out_date) and not cfg.db_overwrite:
                                if ((basin == cfg.inputs_basins[0]) and
                                    (variable == cfg.inputs_variables[0])):
                                    logging.info(' Skipping inputs for {}, {}, '
                                         '{}, database records exist...'
                                         ''.format(variable, basin,
                                         out_date.strftime("%Y-%-m-%-d %H:00")))
                                else:
                                    logging.debug(' Skipping inputs for {}, {}, '
                                         '{}, database records exist...'
                                         ''.format(variable, basin,
                                         out_date.strftime("%Y-%-m-%-d %H:00")))

            swe = copy.deepcopy(cfg.outputs['swe_z'][iters])
            cold = copy.deepcopy(cfg.outputs['coldcont'][iters])

            # Loop over outputs (depths are copied, volumes are calculated)
            for k in proc_list:

                # Mask by subbasin
                for name in cfg.masks:
                    if name == cfg.plotorder[0] and k == list(proc_list.keys())[0]:
                        logging.info(' Processing {}, {}'.format(name,
                                     cfg.rundirs_dict[wy_hour].split('/')[-1]))
                    else:
                        logging.debug(' Processing {}, {}, {}'.format(k, name,
                                     cfg.rundirs_dict[wy_hour].split('/')[-1]))

                    mask = copy.deepcopy(cfg.masks[name]['mask'])
                    elevbin = cfg.ixd*mask

                    snow_mask = swe > 0
                    avail_mask = cold > cfg.cclimit
                    unavail_mask = cold <= cfg.cclimit

                    if k in cfg.variables.awsm_variables:
                        o = copy.deepcopy(cfg.outputs[k][iters])

                    # elevation band
                    for n in np.arange(0,len(cfg.edges)):
                        b = cfg.edges[n]
                        elev_mask = elevbin == n

                        # list of masks
                        be = [mask, elev_mask]
                        ba = [mask, avail_mask]
                        bu = [mask, unavail_mask]
                        bes = [mask, elev_mask, snow_mask]
                        besa = [mask, elev_mask, snow_mask, avail_mask]
                        besu = [mask, elev_mask, snow_mask, unavail_mask]

                        if k == 'swe_z':
                            variables['swe_unavail']['df'].loc[b,name] = calculate(o, cfg.pixel, besu, 'sum','volume', cfg.units, cfg.dplcs)
                            variables['swe_avail']['df'].loc[b,name] = calculate(o, cfg.pixel, besa, 'sum','volume', cfg.units, cfg.dplcs)
                            variables['swe_vol']['df'].loc[b,name] = calculate(o, cfg.pixel, be, 'sum','volume', cfg.units, cfg.dplcs)
                            variables[k]['df'].loc[b,name] = calculate(o, cfg.pixel, be, 'mean','depth', cfg.units, cfg.dplcs)

                        if k == 'swi_z':
                            variables[k]['df'].loc[b,name] = calculate(o,cfg.pixel,be,'mean','depth', cfg.units, cfg.dplcs)
                            variables['swi_vol']['df'].loc[b,name] = calculate(o,cfg.pixel,be,'sum','volume', cfg.units, cfg.dplcs)

                        if k in cfg.variables.process_depth_units:

                            # iSnobal depth units are m
                            if k == 'depth':
                                type = 'snow_depth'
                            else:
                                type = variables[k]['unit_type']

                            variables[k]['df'].loc[b,name] = calculate(o, cfg.pixel, bes, variables[k]['calculate'], type, cfg.units, cfg.dplcs)

                            if out_date == cfg.outputs['dates'][-1] and k == 'density':
                                od = copy.deepcopy(o)
                                ml = [mask, elev_mask, snow_mask]
                                for m in ml:
                                    m = m.astype('float')
                                    m[m < 1] = np.nan
                                    od = od * m

                                density[name][cfg.edges[n]] = copy.deepcopy(od)

                        if k in ['precip', 'rain'] :
                            variables[k]['df'].loc[b,name] = calculate(pre, cfg.pixel, be, variables[k]['calculate'],variables[k]['unit_type'], cfg.units, cfg.dplcs)
                            if k == 'precip':
                                variables['precip_vol']['df'].loc[b,name] = calculate(pre, cfg.pixel, be, 'sum', 'volume', cfg.units, cfg.dplcs)

                    if k in cfg.variables.process_depth_units:

                        if k == 'depth':
                            type = 'snow_depth'
                        else:
                            type = variables[k]['unit_type']
                        variables[k]['df'].loc['total',name] = calculate(o, cfg.pixel, [mask, snow_mask], variables[k]['calculate'],type, cfg.units, cfg.dplcs)

                    if k == 'swe_z':
                        variables['swe_vol']['df'].loc['total',name] = calculate(o, cfg.pixel, mask, 'sum','volume', cfg.units, cfg.dplcs)
                        variables['swe_avail']['df'].loc['total',name] = calculate(o, cfg.pixel, ba, 'sum','volume', cfg.units, cfg.dplcs)
                        variables['swe_unavail']['df'].loc['total',name] = calculate(o, cfg.pixel, bu, 'sum','volume', cfg.units, cfg.dplcs)
                        variables[k]['df'].loc['total',name] = calculate(o, cfg.pixel, mask, 'mean','depth', cfg.units, cfg.dplcs)
                        variables['snow_line']['df'].loc['total',name] = snow_line(o, cfg.dem, mask, cfg.diagnostics_flag)

                    if k == 'swi_z':
                        variables[k]['df'].loc['total',name] = calculate(o,cfg.pixel,mask,'mean','depth', cfg.units, cfg.dplcs)
                        variables['swi_vol']['df'].loc['total',name] = calculate(o,cfg.pixel,mask,'sum','volume', cfg.units, cfg.dplcs)

                    if k in ['precip', 'rain']:
                        variables[k]['df'].loc['total',name] = calculate(pre, cfg.pixel, mask, 'mean', 'depth', cfg.units, cfg.dplcs)
                        if k == 'precip':
                            variables['precip_vol']['df'].loc['total',name] = calculate(pre, cfg.pixel, mask, 'sum', 'volume', cfg.units, cfg.dplcs)

                df = copy.deepcopy(variables[k]['df'])
                df = df.round(decimals = cfg.dplcs)

                if not pass_flag:
                    package(cfg.connector, lbls, cfg.basins, df, cfg.run_id,
                            cfg.vid, k, out_date, cfg.run_name)

                    # snow line
                    if k == 'swe_z':
                        package(cfg.connector, lbls, cfg.basins, df,
                                cfg.run_id, cfg.vid, 'snow_line', out_date,
                                cfg.run_name)

            stamp = datetime.now().strftime("%Y-%-m-%-d %H:%M:%S")
            if not pass_flag:
                logging.info(' Completed {}, date: {}, at: {} UTC'.format(
                     cfg.rundirs_dict[wy_hour].split('/')[-1],
                     str(wy_hour - elapsed_hours),
                     out_date.date().strftime("%Y-%-m-%-d"),
                     stamp))

                if iters != 0 and int(wy_hour - elapsed_hours) != 24:
                    logging.warning(' Elapsed hours = {}, there may be a '
                                    'missing output '
                                    'day'.format(str(wy_hour-elapsed_hours)))

            elapsed_hours = int(wy_hour)

        self.density = density
        self.precip_total = precip_total
        self.rain_total = rain_total
        self.precip_flag = flag
