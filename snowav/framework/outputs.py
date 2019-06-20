
import copy
import os
import numpy as np
from datetime import datetime
from snowav.utils.wyhr import calculate_wyhr_from_date
from snowav.utils.OutputReader import iSnobalReader
import netCDF4 as nc

def outputs(run_dirs = None, start_date = None, end_date = None,
            filetype = None, wy = None):
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

    Returns
    ------
    outputs : dict
        dictionary of snow.nc and em.nc outputs within time period
    run_dirs : list
        modified run_dirs, with paths outside of start_date, end_date removed
    rdict : dict
        process() lookup

    '''

    outputs = {'swi_z':[], 'evap_z':[], 'snowmelt':[], 'swe_z':[],'depth':[],
               'dates':[], 'time':[], 'density':[], 'coldcont':[] }

    rdict = {}

    dirs = copy.deepcopy(run_dirs)

    for path in dirs:

        snowfile = os.path.join(path, 'snow.nc')

        # If the run_dirs isn't empty use it, otherwise remove
        if not os.path.isfile(snowfile):
            raise Exception('{} not a valid file'.format(snowfile))

        ncf = nc.Dataset(snowfile)
        t = nc.num2date(ncf.variables['time'][0],ncf.variables['time'].units)
        ncf.close()

        # Only load the rundirs that we need
        if (t.date() >= start_date.date()) and (t.date() <= end_date.date()):

            # pass the reader start and end times
            if (start_date is not None) and (end_date is not None):
                st_hr = calculate_wyhr_from_date(start_date)
                en_hr = calculate_wyhr_from_date(end_date)

            else:
                st_hr = calculate_wyhr_from_date(t)
                en_hr = calculate_wyhr_from_date(t)

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

            # Everything but 'dates' gets clipped in the reader
            outputs['dates'] = np.asarray(([d for (d, remove) in
                                    zip(outputs['dates'],
                                    (outputs['dates'] > end_date))
                                    if not remove]))
            outputs['dates'] = np.asarray(([d for (d, remove) in
                                    zip(outputs['dates'],
                                    (outputs['dates'] < start_date))
                                    if not remove]))

        else:
            run_dirs.remove(path)

    return outputs, run_dirs, rdict
