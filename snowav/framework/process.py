from copy import deepcopy
from datetime import datetime
import logging
import numpy as np
import os
import pandas as pd
from sys import exit
import warnings

from snowav.database.database import package
from snowav.database.tables import Results, RunMetadata, Inputs
from snowav.utils.utilities import calculate, sum_precip, snow_line, \
    input_summary


class Process(object):
    """ Process AWSM model results made in awsm_daily format.

    Args
    ------
    cfg: config class
    db: database class
    """

    def __init__(self, cfg, db):

        pd.options.mode.chained_assignment = None
        warnings.filterwarnings('ignore')

        elapsed_hours = 0
        variables = cfg.variables.variables

        if os.path.splitext(cfg.connector)[1] == '.db':
            dbs = 'sqlite'
        else:
            dbs = 'sql'

        # images for figures
        precip_total = np.zeros((cfg.nrows, cfg.ncols))
        rain_total = np.zeros((cfg.nrows, cfg.ncols))

        # Check that topo and outputs are the same dimensions
        if cfg.outputs['swe_z'][0].shape != cfg.dem.shape:
            logging.error(" topo.nc dimensions: {} do not match output "
                          "dimensions: "
                          "{}".format(cfg.dem.shape,
                                      cfg.outputs['swe_z'][0].shape))
            exit()

        # process each date
        for iters, out_date in enumerate(cfg.outputs['dates']):
            wy_hour = int(cfg.outputs['time'][iters])
            dir_str = cfg.rundirs_dict[wy_hour].split('/')[-1]
            odate_str = out_date.strftime("%Y-%-m-%-d %H:00")

            # assume to start that we are processing everything
            proc_list = deepcopy(cfg.variables.snowav_results_variables)
            pass_flag = False

            # delete anything on the database with the same run_name and date
            for bid in cfg.plotorder:
                params = {Results: [('date_time', 'eq', out_date),
                                    ('basin_id', 'eq', cfg.basins[bid]['basin_id']),
                                    ('run_id', 'ge', 0)
                                    ],
                          RunMetadata: [('run_name', 'eq', cfg.run_name)]}
                results = db.query(params)

                if len(results['run_id'].unique()) > 1:
                    cfg._logger.warn(" Multiple database 'run_id' values "
                                     "returned for single query")

                if not results.empty:
                    if not cfg.db_overwrite:
                        # If database records exist we can skip all outputs
                        # but density on the last day, if making the figure
                        pass_flag = True

                        if bid == cfg.plotorder[0]:
                            logging.info(' Skipping outputs for {}, {}, '
                                         'database records exist'
                                         '...'.format(bid, odate_str))
                        else:
                            logging.debug(' Skipping outputs for {}, {}, '
                                          'database records exist'
                                          '...'.format(bid, odate_str))

                        if (out_date == cfg.outputs['dates'][-1] and
                                cfg.density_flag):
                            proc_list = ['density']
                        else:
                            proc_list = []

                    else:
                        del_params = {
                            Results: [
                                ('date_time', 'eq', out_date),
                                ('basin_id', 'eq', cfg.basins[bid]['basin_id']),
                                ('run_id', 'in', [int(x) for x in results['run_id'].unique()])
                            ]
                        }

                        db.delete(del_params)
                        logging.debug(' Deleting database records for '
                                      '{}, {}'.format(bid, odate_str))

                        # query again, and if we have deleted everything from
                        # that run, also delete the metadata
                        results = db.query(params)
                        del_params = {
                            RunMetadata: [
                                ('run_id', 'in', [int(x) for x in results['run_id'].unique()])
                            ]
                        }

                        if results.empty:
                            db.delete(del_params)

            density = {}
            for name in cfg.masks:
                density[name] = {}

            # if database records don't exist, or we are making the
            # precip_depth figure we need to process precip
            if not pass_flag or cfg.precip_depth_flag:
                logging.info(' Processing precip {}'.format(
                    dir_str))
                run_dir = cfg.rundirs_dict[cfg.outputs['time'][iters]]
                precip_path = os.path.join(
                    run_dir.replace('/runs/run', '/data/data'),
                    'smrfOutputs/precip.nc')

                percent_snow_path = precip_path.replace('precip.nc',
                                                        'percent_snow.nc')

                if (os.path.isfile(precip_path) and
                        os.path.isfile(percent_snow_path)):
                    precip, rain = sum_precip(precip_path,
                                              percent_snow_path)
                    precip_total = precip_total + precip
                    rain_total = rain_total + rain
                else:
                    logging.warning(' One or both of precip.nc, '
                                    'percent_snow.nc does not exist')

            if cfg.inputs_flag:
                mask_list = []
                for name in cfg.masks:
                    mask_list.append(deepcopy(cfg.masks[name]['mask']))

                # sf = cfg.rundirs_dict[wy_hour].replace('runs', 'data')
                # sf = sf.replace('run', 'data') + '/smrfOutputs/'
                sf = cfg.rundirs_dict[wy_hour].replace('/runs/run', '/data/data') + '/smrfOutputs/'

                params = {Inputs: [
                    ('date_time', 'eq', out_date),
                    ('basin_id', 'eq', cfg.basins[bid]['basin_id']),
                    ('run_name', 'eq', cfg.run_name),
                    ('variable', 'in', cfg.variables.snowav_inputs_variables)
                ],
                }

                df = db.query(params)
                idx = df.index

                for i, input in enumerate(cfg.variables.snowav_inputs_variables):
                    input_path = os.path.join(sf, input + '.nc')

                    for basin in cfg.inputs_basins:
                        basin_id = int(cfg.basins[basin]['basin_id'])

                        if (not idx.contains(out_date) or
                                (idx.contains(out_date) and cfg.db_overwrite)):
                            if basin == cfg.inputs_basins[0] and i == 0:
                                logging.info(' Processing {} inputs, '
                                             '{} '.format(basin, odate_str))
                            else:
                                logging.debug(' Processing {}, {}, '
                                              '{}'.format(input, basin,
                                                          odate_str))

                            mask = cfg.masks[basin]['mask']

                            input_summary(input_path,
                                          input,
                                          cfg.inputs_methods,
                                          cfg.inputs_percentiles,
                                          dbs,
                                          cfg.connector,
                                          cfg.run_name,
                                          basin_id,
                                          cfg.run_id,
                                          mask)

                        if idx.contains(out_date) and not cfg.db_overwrite:
                            if ((basin == cfg.inputs_basins[0]) and
                                    (input == cfg.inputs_variables[0])):
                                logging.info(' Skipping inputs for {}, {}, {},'
                                             ' database records exist...'
                                             ''.format(input,
                                                       basin,
                                                       odate_str))
                            else:
                                logging.debug(' Skipping inputs for {}, {}, '
                                              '{}, database records exist...'
                                              ''.format(input,
                                                        basin,
                                                        odate_str))

            swe = deepcopy(cfg.outputs['swe_z'][iters])
            cold = deepcopy(cfg.outputs['coldcont'][iters])

            # Loop over outputs (depths are copied, volumes are calculated)
            for k in proc_list:

                # Mask by subbasin
                for name in cfg.masks:
                    if name == cfg.plotorder[0] and k == proc_list[0]:
                        logging.info(' Processing {}, {}'.format(name,
                                                                 dir_str))
                    else:
                        logging.debug(' Processing {}, {}, {}'.format(k, name,
                                                                      dir_str))

                    mask = deepcopy(cfg.masks[name]['mask'])
                    elevbin = cfg.ixd * mask

                    snow_mask = swe > 0
                    avail_mask = cold > cfg.cclimit
                    unavail_mask = cold <= cfg.cclimit

                    if k in cfg.variables.awsm_variables:
                        o = deepcopy(cfg.outputs[k][iters])

                    # elevation band
                    for n in np.arange(0, len(cfg.edges)):
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
                            variables['swe_unavail']['df'].loc[b, name] = \
                                calculate(o, cfg.pixel, besu, 'sum', 'volume',
                                          cfg.units, cfg.dplcs)
                            variables['swe_avail']['df'].loc[b, name] = \
                                calculate(o, cfg.pixel, besa, 'sum', 'volume',
                                          cfg.units, cfg.dplcs)
                            variables['swe_vol']['df'].loc[b, name] = \
                                calculate(o, cfg.pixel, be, 'sum', 'volume',
                                          cfg.units, cfg.dplcs)
                            variables[k]['df'].loc[b, name] = \
                                calculate(o, cfg.pixel, be, 'mean', 'depth',
                                          cfg.units, cfg.dplcs)

                        if k == 'swi_z':
                            variables[k]['df'].loc[b, name] = \
                                calculate(o, cfg.pixel, be, 'mean', 'depth',
                                          cfg.units, cfg.dplcs)
                            variables['swi_vol']['df'].loc[b, name] = \
                                calculate(o, cfg.pixel, be, 'sum', 'volume',
                                          cfg.units, cfg.dplcs)

                        if k in cfg.variables.process_depth_units:

                            # iSnobal depth units are m
                            if k == 'depth':
                                type = 'snow_depth'
                            else:
                                type = variables[k]['unit_type']

                            variables[k]['df'].loc[b, name] = \
                                calculate(o, cfg.pixel, bes,
                                          variables[k]['calculate'], type,
                                          cfg.units, cfg.dplcs)

                            if (out_date == cfg.outputs['dates'][-1] and
                                    k == 'density'):
                                od = deepcopy(o)
                                ml = [mask, elev_mask, snow_mask]
                                for m in ml:
                                    m = m.astype('float')
                                    m[m < 1] = np.nan
                                    od = od * m

                                density[name][cfg.edges[n]] = deepcopy(od)

                        if k == 'precip_z':
                            variables['precip_z']['df'].loc[b, name] = \
                                calculate(precip, cfg.pixel, be,
                                          variables[k]['calculate'],
                                          variables[k]['unit_type'], cfg.units,
                                          cfg.dplcs)
                            variables['rain_z']['df'].loc[b, name] = \
                                calculate(rain, cfg.pixel, be,
                                          variables[k]['calculate'],
                                          variables[k]['unit_type'], cfg.units,
                                          cfg.dplcs)
                            variables['precip_vol']['df'].loc[b, name] = \
                                calculate(precip, cfg.pixel, be, 'sum',
                                          'volume', cfg.units, cfg.dplcs)

                    if k in cfg.variables.process_depth_units:

                        if k == 'depth':
                            type = 'snow_depth'
                        else:
                            type = variables[k]['unit_type']
                        variables[k]['df'].loc['total', name] = \
                            calculate(o, cfg.pixel, [mask, snow_mask],
                                      variables[k]['calculate'], type,
                                      cfg.units, cfg.dplcs)

                    if k == 'swe_z':
                        variables['swe_vol']['df'].loc['total', name] = \
                            calculate(o, cfg.pixel, mask, 'sum', 'volume',
                                      cfg.units, cfg.dplcs)
                        variables['swe_avail']['df'].loc['total', name] = \
                            calculate(o, cfg.pixel, ba, 'sum', 'volume',
                                      cfg.units, cfg.dplcs)
                        variables['swe_unavail']['df'].loc['total', name] = \
                            calculate(o, cfg.pixel, bu, 'sum', 'volume',
                                      cfg.units, cfg.dplcs)
                        variables[k]['df'].loc['total', name] = \
                            calculate(o, cfg.pixel, mask, 'mean', 'depth',
                                      cfg.units, cfg.dplcs)

                        if 'snow_line' in proc_list:
                            variables['snow_line']['df'].loc['total', name] = \
                                snow_line(o, cfg.dem, mask,
                                          cfg.diagnostics_flag)

                    if k == 'swi_z':
                        variables[k]['df'].loc['total', name] = \
                            calculate(o, cfg.pixel, mask, 'mean', 'depth',
                                      cfg.units, cfg.dplcs)
                        variables['swi_vol']['df'].loc['total', name] = \
                            calculate(o, cfg.pixel, mask, 'sum', 'volume',
                                      cfg.units, cfg.dplcs)

                    if k == 'precip_z':
                        variables['precip_z']['df'].loc['total', name] = \
                            calculate(precip, cfg.pixel, mask, 'mean', 'depth',
                                      cfg.units, cfg.dplcs)
                        variables['precip_vol']['df'].loc['total', name] = \
                            calculate(precip, cfg.pixel, mask, 'sum', 'volume',
                                      cfg.units, cfg.dplcs)
                        variables['rain_z']['df'].loc['total', name] = \
                            calculate(rain, cfg.pixel, mask, 'mean', 'depth',
                                      cfg.units, cfg.dplcs)

                df = deepcopy(variables[k]['df'])
                df = df.round(decimals=cfg.dplcs)

                if k == 'snow_line':
                    df = df[df.index == 'total']

                if not pass_flag:
                    package(cfg.connector, cfg.basins, df, cfg.run_id, cfg.vid,
                            k, out_date)

            stamp = datetime.now().strftime("%Y-%-m-%-d %H:%M:%S")
            if not pass_flag:
                logging.info(' Completed {} at {} UTC'.format(dir_str, stamp))

                if iters != 0 and int(wy_hour - elapsed_hours) != 24:
                    logging.warning(' Elapsed hours: {}, there may be a '
                                    'missing output '
                                    'day'.format(str(wy_hour - elapsed_hours)))

            elapsed_hours = int(wy_hour)

        self.density = density
        self.precip_total = precip_total
        self.rain_total = rain_total
