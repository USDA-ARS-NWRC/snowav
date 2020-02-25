
from copy import deepcopy
import os
import numpy as np
from datetime import datetime, timedelta
from snowav.utils.wyhr import calculate_wyhr_from_date
from snowav.utils.OutputReader import iSnobalReader
import netCDF4 as nc

def outputs(run_dirs, wy, properties, start_date = None,
            end_date = None, flight_dates = None, loglevel = None):
    '''
    This uses start_date and end_date to load the snow.nc and em.nc of interest
    within a report period to the outputs format that will be used in process().
    Also returns a clipped run_dirs that only contains paths with the specified
    date range. If start_date and end_date are not supplied, no run_dirs will
    be clipped.

    Note: this is assuming awsm_daily output folders and dates in the snow.nc
    file.

    Args
    -----
    run_dirs : list
        list of run directories
    start_date : datetime
        report period start date (optional)
    end_date : datetime
        report period end date (optional)
    wy : int
        water year
    flight_dates : array
        (optional)

    Returns
    ------
    results : dict

    '''

    dirs = deepcopy(run_dirs)
    snowbands = []
    embands = []
    log = []
    rdict = {}
    outputs = {'dates': [], 'time': []}

    bands_map = {'snow':{'depth': 0,
                         'density': 1,
                         'swe_z': 2,
                         'lwc': 3,
                         'temp_surface': 4,
                         'temp_lower': 5,
                         'temp_bulk': 6,
                         'depth_lower_layer': 7,
                         'h20_sat': 8},
                 'em':{'R_n': 0,
                       'H': 1,
                       'L_v_E': 2,
                       'G': 3,
                       'M': 4,
                       'delta_Q': 5,
                       'evap_z': 6,
                       'melt': 7,
                       'swi_z': 8,
                       'coldcont': 9}}

    for band in properties:
        if band in [*bands_map['snow'].keys()]:
            snowbands.append(bands_map['snow'][band])
        if band in [*bands_map['em'].keys()]:
            embands.append(bands_map['em'][band])
        outputs[band] = []

    start = deepcopy(start_date)
    end = deepcopy(end_date)

    # Run this with standard processing, and forecast processing
    if flight_dates is None:

        for path in dirs:
            snowfile = os.path.join(path, 'snow.nc')

            if loglevel == 'DEBUG':
                log.append(' Reading date: {}'.format(snowfile))

            # Consider making this a warning, with an else: .remove(path)
            # to catch other files that are in these directories
            if not os.path.isfile(snowfile):
                log.append(' {} not a valid file'.format(snowfile))
                print(' {} not a valid file, snowav may '
                    'error...'.format(snowfile))
                run_dirs.remove(path)
                results = {'outputs': outputs,
                           'dirs': dirs,
                           'run_dirs': run_dirs,
                           'rdict': rdict,
                           'log': log}
                return results

            ncf = nc.Dataset(snowfile)

            # Catch 'empty' snow.nc and em.nc file from certain awsm crash
            # scenarios in awsm<=0.10.0
            if 'specific_mass' not in ncf.variables:
                log.append(' No "specific_mass" variable in {}, this may be the result '
                    'of awsm crashing without writing variables to file, '
                    'consider deleting and re-running awsm'.format(snowfile))
                raise Exception(' No "specific_mass" variable in {}'.format(snowfile))

            ta = nc.num2date(ncf.variables['time'][:],ncf.variables['time'].units)
            ncf.close()

            ta = np.sort(ta)

            if start_date is None:
                start = deepcopy(ta[0])

            if end_date is None:
                end = deepcopy(ta[-1])

            for idx,t in enumerate(ta):

                # Only load the rundirs that we need
                if (t.date() >= start.date()) and (t.date() <= end.date()):

                    log.append(' Loading: {}'.format(snowfile))

                    st_hr = calculate_wyhr_from_date(start)
                    en_hr = calculate_wyhr_from_date(end)

                    output = iSnobalReader(path, snowbands = snowbands,
                                           embands = embands, wy = wy,
                                           time_start = st_hr, time_end = en_hr)

                    # Make a dict for wyhr-rundir lookup
                    for ot in output.time:
                        rdict[int(ot)] = path

                    # Apply snow bands to outputs
                    for band in properties:
                        if band in [*bands_map['snow'].keys()]:
                            s = bands_map['snow'][band]
                            outputs[band].append(output.snow_data[s][idx,:,:])
                        if band in [*bands_map['em'].keys()]:
                            e = bands_map['em'][band]
                            outputs[band].append(output.em_data[e][idx,:,:])

                    outputs['dates'].append(output.dates[idx])
                    outputs['time'].append(output.time[idx])

                else:
                    run_dirs.remove(path)

    # Run this when flight updates are present to make custom outputs
    else:

        for path in dirs:
            snowfile = os.path.join(path, 'snow.nc')

            # If the run_dirs isn't empty use it, otherwise remove
            if not os.path.isfile(snowfile):
                raise Exception('{} not a valid file'.format(snowfile))

            ncf = nc.Dataset(snowfile)
            ta = nc.num2date(ncf.variables['time'][:],ncf.variables['time'].units)
            ncf.close()

            for idx,t in enumerate(ta):
                if (t.date() in [x.date() for x in flight_dates]):

                    output = iSnobalReader(path, snowbands=[0,1,2], wy=wy)

                    for ot in output.time:
                        rdict[int(ot)] = path

                    outputs['swe_z'].append(output.snow_data[2][idx,:,:])
                    outputs['depth'].append(output.snow_data[0][idx,:,:])
                    outputs['density'].append(output.snow_data[1][idx,:,:])
                    outputs['dates'].append(output.dates[idx])
                    outputs['time'].append(output.time[idx])

    results = {'outputs': outputs,
               'dirs': dirs,
               'run_dirs': run_dirs,
               'rdict': rdict,
               'log': log}

    return results
