
import copy
import os
import numpy as np
from datetime import datetime, timedelta
from snowav.utils.wyhr import calculate_wyhr_from_date
from snowav.utils.OutputReader import iSnobalReader
import netCDF4 as nc

def outputs(run_dirs = None, start_date = None, end_date = None,
            filetype = None, wy = None, flight_dates = None, loglevel = None):
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
    filetype : str
        currently only supporting 'netcdf'
    wy : int
        water year
    flight_dates : array
        (optional)

    Returns
    ------
    outputs : dict
        dictionary of snow.nc and em.nc outputs within time period
    dirs : list
        all dirs within run_dir
    run_dirs : list
        modified run_dirs, with paths outside of start_date, end_date removed
    rdict : dict
        process() lookup

    '''

    log = []
    rdict = {}
    dirs = copy.deepcopy(run_dirs)
    outputs = {'swi_z':[], 'evap_z':[], 'snowmelt':[], 'swe_z':[],'depth':[],
               'dates':[], 'time':[], 'density':[], 'coldcont':[] }

    start = copy.deepcopy(start_date)
    end = copy.deepcopy(end_date)

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
                return [], [], [], [], log

            ncf = nc.Dataset(snowfile)
            ta = nc.num2date(ncf.variables['time'][:],ncf.variables['time'].units)
            ncf.close()

            ta = np.sort(ta)

            if start_date is None:
                start = copy.deepcopy(ta[0])

            if end_date is None:
                end = copy.deepcopy(ta[-1])

            for idx,t in enumerate(ta):

                # Only load the rundirs that we need
                if (t.date() >= start.date()) and (t.date() <= end.date()):

                    log.append(' Loading: {}'.format(snowfile))

                    st_hr = calculate_wyhr_from_date(start)
                    en_hr = calculate_wyhr_from_date(end)

                    output = iSnobalReader(path, filetype, snowbands = [0,1,2],
                                           embands = [6,7,8,9], wy = wy,
                                           time_start = st_hr, time_end = en_hr)

                    # Make a dict for wyhr-rundir lookup
                    for ot in output.time:
                        rdict[int(ot)] = path

                    outputs['swi_z'].append(output.em_data[8][idx,:,:])
                    outputs['snowmelt'].append(output.em_data[7][idx,:,:])
                    outputs['evap_z'].append(output.em_data[6][idx,:,:])
                    outputs['coldcont'].append(output.em_data[9][idx,:,:])
                    outputs['swe_z'].append(output.snow_data[2][idx,:,:])
                    outputs['depth'].append(output.snow_data[0][idx,:,:])
                    outputs['density'].append(output.snow_data[1][idx,:,:])
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

                    output = iSnobalReader(path, filetype, snowbands=[0,1,2], wy=wy)

                    for ot in output.time:
                        rdict[int(ot)] = path

                    outputs['swe_z'].append(output.snow_data[2][idx,:,:])
                    outputs['depth'].append(output.snow_data[0][idx,:,:])
                    outputs['density'].append(output.snow_data[1][idx,:,:])
                    outputs['dates'].append(output.dates[idx])
                    outputs['time'].append(output.time[idx])

    return outputs, dirs, run_dirs, rdict, log
