
import copy
import os
import numpy as np
from datetime import datetime, timedelta
from snowav.utils.wyhr import calculate_wyhr_from_date
from snowav.utils.OutputReader import iSnobalReader
import netCDF4 as nc

def outputs(run_dirs = None, start_date = None, end_date = None,
            filetype = None, wy = None, flight_dates = None,
            pre_flight_dates = None):
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
    pre_flight_dates : array
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

            # Consider making this a warning, with an else: .remove(path)
            # to catch other files that are in these directories
            if not os.path.isfile(snowfile):
                raise Exception('{} not a valid file'.format(snowfile))

            ncf = nc.Dataset(snowfile)
            t = nc.num2date(ncf.variables['time'][0],ncf.variables['time'].units)
            ncf.close()

            if start_date is None:
                start = copy.deepcopy(t)

            if end_date is None:
                end = copy.deepcopy(t)

            # Only load the rundirs that we need
            if (t.date() >= start.date()) and (t.date() <= end.date()):

                st_hr = calculate_wyhr_from_date(start)
                en_hr = calculate_wyhr_from_date(end)

                output = iSnobalReader(path, filetype, snowbands = [0,1,2],
                                       embands = [6,7,8,9], wy = wy,
                                       time_start = st_hr, time_end = en_hr)

                outputs['dates'] = np.append(outputs['dates'],output.dates)
                outputs['time'] = np.append(outputs['time'],output.time)

                # Make a dict for wyhr-rundir lookup
                for t in output.time:
                    rdict[int(t)] = path

                for n in range(0,len(output.em_data[8])):
                    outputs['swi_z'].append(output.em_data[8][n,:,:])
                    outputs['snowmelt'].append(output.em_data[7][n,:,:])
                    outputs['evap_z'].append(output.em_data[6][n,:,:])
                    outputs['coldcont'].append(output.em_data[9][n,:,:])
                    outputs['swe_z'].append(output.snow_data[2][n,:,:])
                    outputs['depth'].append(output.snow_data[0][n,:,:])
                    outputs['density'].append(output.snow_data[1][n,:,:])

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
            t = nc.num2date(ncf.variables['time'][0],ncf.variables['time'].units)
            ncf.close()

            # Only load the rundirs that we need
            if ((t.date() in [x.date() for x in flight_dates]) and
               (t.date() <= end_date.date())):

                output = iSnobalReader(path, filetype, snowbands=[0,1,2], wy=wy)
                outputs['dates'] = np.append(outputs['dates'],output.dates)
                outputs['time'] = np.append(outputs['time'],output.time)

                for t in output.time:
                    rdict[int(t)] = path

                for n in range(0,len(output.snow_data[0])):
                    outputs['swe_z'].append(output.snow_data[2][n,:,:])
                    outputs['depth'].append(output.snow_data[0][n,:,:])
                    outputs['density'].append(output.snow_data[1][n,:,:])

            # Always get the previous snow.nc to do flight difference
            if ((t.date() in [x.date()-timedelta(hours = 24) for x in pre_flight_dates]) and
                (t.date() <= end_date.date())):

                output = iSnobalReader(path, filetype, snowbands=[0,1,2], wy=wy)
                outputs['dates'] = np.append(outputs['dates'],output.dates)
                outputs['time'] = np.append(outputs['time'],output.time)

                for n in range(0,len(output.snow_data[0])):
                    outputs['swe_z'].append(output.snow_data[2][n,:,:])
                    outputs['depth'].append(output.snow_data[0][n,:,:])
                    outputs['density'].append(output.snow_data[1][n,:,:])

            else:
                run_dirs.remove(path)

    return outputs, dirs, run_dirs, rdict
