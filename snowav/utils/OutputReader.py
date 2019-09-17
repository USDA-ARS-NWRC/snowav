import pandas as pd
import numpy as np
from spatialnc import ipw
from snowav.utils.wyhr import calculate_wyhr_from_date
import os
import glob
import pytz
from netCDF4 import Dataset
import netCDF4 as nc
import warnings

class iSnobalReader():
    def __init__(self, outputdir, filetype='netcdf', timesteps=None, embands=None,
                snowbands=None, mask=None, time_start=None, time_end=None,
                wy=None):
        """
        Inputs:
        outputdir - abosulte path to location of outputs
        filetype - either ipw or netcdf
        timestep - list of time steps to return, default all
        embands - list of bands numbers to grab, default all
        snowbands - list of bands numbers to grab, default all
        time_start - water year hour to start grabbing files
        time_end - water year hour to end grabbing files
        wy - water year for finding dates from ipw file wyhrs
        """

        # parse innputs
        self.outputdir = outputdir
        self.filetype = filetype
        self.timesteps = timesteps
        self.mask = mask
        self.time_start = time_start
        self.time_end = time_end
        self.wy = wy

        # list of band numbers possible
        emnums = range(10)
        snownums = range(9)
        self.embands = embands
        self.snowbands = snowbands
        # if no bands specified, get them all
        if self.embands is None:
            self.embands = emnums
        if self.snowbands is None:
            self.snowbands = snownums

        emnames = ['net_rad', 'sensible_heat', 'latent_heat', 'snow_soil',
                    'precip_advected', 'sum_EB', 'evaporation',
                    'snowmelt', 'SWI', 'cold_content'
                  ]
        snownames = ['thickness', 'snow_density', 'specific_mass', 'liquid_water',
                    'temp_surf', 'temp_lower', 'temp_snowcover',
                    'thickness_lower', 'water_saturation'
                    ]

        self.emdict = {}
        self.snowdict = {}
        for idn, nm in enumerate(emnames):
            self.emdict[emnums[idn]] = nm
        for idn, nm in enumerate(snownames):
            self.snowdict[snownums[idn]] = nm

        # read in and store the data
        self.snow_data, self.em_data = self.read_isnobal_outputs()

    def read_isnobal_outputs(self):
        """
        Inputs:
        self

        Returns:
        snow_data - dictionary of time series spatial data np arrays as well as a
                    time series array
        em_date - dictionary of time series spatial data np arrays as well as a
                    time series array
        """

        # check for errors in time inputs
        if self.timesteps is not None and (self.time_start is not None or self.time_end is not None):
            print('timesteps: {}, time_start: {}, time_end: {}'.format(self.timesteps, self.time_start, self.time_end))
            raise ValueError('Incorrect combination of time inputs for reading IPW files')
        if (self.time_start is None and self.time_end is not None) or (self.time_start is None and self.time_end is not None):
            print('time_start: {}, time_end: {}'.format(self.time_start, self.time_end))
            raise ValueError('Incorrect combination of time inputs for reading IPW files')

        # if we're dealing with ipw outputs
        if self.filetype == 'ipw':
            if self.wy is None:
                raise IOError('WY not included in iSnobalReader call with type ipw')
            # set start date for wyhrs from ipw files to allow datetime conversion
            self.start_ipw = pd.to_datetime('10-01-{} 00:00'.format(self.wy-1))
            # get the data
            snow_data, em_data = self.read_isnobal_ipw()

        # if we're dealing with netcdf output files
        elif self.filetype == 'netcdf':
            snow_data, em_data = self.read_isnobal_netcdf()

        # if not a valid option
        else:
            raise ValueError('specified filetype incorrect for readings iSnobal outputs')

        return snow_data, em_data

    def read_isnobal_ipw(self):
        """
        Read in specific timesteps and bands from ipw files output from isnobal.

        Input:
        self

        Returns:
        snow_data - dictionary of time series spatial data np arrays as well as a
                    time series array
        em_date - dictionary of time series spatial data np arrays as well as a
                    time series array
        """
        # path to snow and em files
        pathsnow = os.path.join(self.outputdir, "snow.*".format(self.outputdir))
        pathem = os.path.join(self.outputdir, "em.*".format(self.outputdir))

        # get list of isnobal outputs and sort by time step
        snowfiles = sorted(glob.glob(pathsnow), key=os.path.getmtime)
        snowfiles.sort(key=lambda f: os.path.basename(f).split('snow.')[1])

        emfiles = sorted(glob.glob(pathem), key=os.path.getmtime)
        emfiles.sort(key=lambda f: os.path.basename(f).split('em.')[1])

        # create list of timesteps if none given
        if self.timesteps is None and self.time_start is None and self.time_end is None:
            timesteps = []
            for nm in snowfiles:
                ts = os.path.basename(nm).split('snow.')[1]
                ts = int(ts)
                timesteps.append(ts)

        # if we are grabbing files between 2 times
        if self.time_start is not None and self.time_end is not None and self.timesteps is None:
            timesteps = []
            for nm in snowfiles:
                ts = os.path.basename(nm).split('snow.')[1]
                ts = int(ts)
                if ts >= self.time_start and ts <= self.time_end:
                    timesteps.append(ts)

        # empty dictionaries to store numpy arrays of each variable
        snow_data = {}
        em_data = {}
        # read in one file to get size
        size_var = ipw.IPW(emfiles[0]).bands[0].data.shape
        if self.mask is None:
            mask = np.ones(size_var)
        else:
            mask = self.mask

        # initilize numpy arrays
        for ebd in self.embands:
            em_data[ebd] = np.zeros((len(timesteps), size_var[0], size_var[1]))
        for sbd in self.snowbands:
            snow_data[sbd] = np.zeros((len(timesteps), size_var[0], size_var[1]))

        # assign numpy array of timesteps
        self.time = np.array(timesteps)
        # get dates for wyhr from ipw files
        ipw_dates = [self.start_ipw + pd.to_timedelta(dts, unit='h') for dts in self.time]
        self.dates = np.array(ipw_dates)

        for idt, ts in enumerate(timesteps):
            # get correct file
            fps = os.path.join(self.outputdir, 'snow.%04i'%(ts))
            fpe = os.path.join(self.outputdir, 'em.%04i'%(ts))
            # read in ipw
            i_snow = ipw.IPW(fps)
            i_em = ipw.IPW(fpe)

            # assign each correct band
            for sbd in self.snowbands:
                tmp_snow_data = i_snow.bands[sbd].data
                if sbd in [4, 5, 6]:
                    tmp_snow_data[tmp_snow_data == -75.0] = np.nan
                #tmp_snow_data[np.isnan(tmp_snow_data)] = 0.0
                snow_data[sbd][idt,:] = tmp_snow_data*mask

            for ebd in self.embands:
                tmp_em_data = i_em.bands[ebd].data
                #tmp_em_data[np.isnan(tmp_em_data)] = 0.0
                em_data[ebd][idt,:] = tmp_em_data*mask


        return snow_data, em_data

    def read_isnobal_netcdf(self):
        """
        Read in specific timesteps and bands from ipw files output from isnobal.

        Input:
        self

        Returns:
        snow_data - dictionary of time series spatial data np arrays as well as a
                    time series array
        em_date - dictionary of time series spatial data np arrays as well as a
                    time series array
        """

        warnings.simplefilter("ignore")

        name_snow_fp = 'snow'
        name_em_fp = 'em'
        # path to snow and em files
        pathsnow = os.path.join(self.outputdir, name_snow_fp+'.nc')
        pathem = os.path.join(self.outputdir, name_em_fp+'.nc')
        # read datasets
        ds_snow = Dataset(pathsnow, 'r')
        ds_em = Dataset(pathem, 'r')

        # hack for different swi names
        if 'runoff' in ds_em.variables:
            self.emdict[8] = 'runoff'

        # get time array from dataset
        ny = len(ds_snow.variables['y'][:])
        nx = len(ds_snow.variables['x'][:])

        # find dates the right way
        nc_time = ds_snow.variables['time'][:]
        t_units = ds_snow.variables['time'].units
        nc_calendar = ds_snow.variables['time'].calendar
        nc_dates = nc.num2date(nc_time, t_units, nc_calendar)
        self.dates = nc_dates

        # calculate water year hours
        time = np.zeros(len(nc_dates))
        for ix, t in enumerate(nc_dates):
            # calculate new date
            dtt = t
            # convert to wyhr
            wyhr = calculate_wyhr_from_date(dtt)
            # store wyhr in time array
            time[ix] = wyhr

        # mask is ones if none specified
        if self.mask is None:
            mask = np.ones((ny,nx))
        else:
            mask = self.mask

        # empty dictionaries to store numpy arrays of each variable
        snow_data = {}
        em_data = {}

        # print('Processing datasets')
        # parse user input timesteps or use all
        if self.timesteps is None and self.time_start is None and self.time_end is None:
            timesteps = time
            itd = len(time)
            for ebd in self.embands:
                e_name = self.emdict[ebd]
                em_data[ebd] = ds_em.variables[e_name][:]*mask
                #em_data[ebd][em_data[ebd] <= -75.0] = np.nan
            for sbd in self.snowbands:
                s_name = self.snowdict[sbd]
                snow_data[sbd] = ds_snow.variables[s_name][:]*mask
                #snow_data[sbd][snow_data[sbd] <= -75.0] = np.nan

        # if given a start and stop time
        elif self.timesteps is None and self.time_start is not None and self.time_end is not None:
            # only time steps in desired timeframe
            timesteps = time[(time >= self.time_start) & (time <= self.time_end)]
            itd = timesteps

            # pull data out at specified times
            for ebd in self.embands:
                # find variable name
                e_name = self.emdict[ebd]

                em_data[ebd] = np.zeros((len(itd), ny, nx))
                for ide, it in enumerate(itd):
                    time_samp = np.where(time == it)[0]
                    tmp_em_data = ds_em.variables[e_name][time_samp,:]
                    #tmp_em_data[np.isnan(tmp_em_data)] = 0.0
                    em_data[ebd][ide,:] = tmp_em_data*mask

            for sbd in self.snowbands:
                # find variable name
                s_name = self.snowdict[sbd]
                snow_data[sbd] = np.zeros((len(itd), ny, nx))
                for ids, it in enumerate(itd):
                    time_samp = np.where(time == it)[0]
                    tmp_snow_data = ds_snow.variables[s_name][time_samp,:]
                    if sbd in [4, 5, 6]:
                        tmp_snow_data[tmp_snow_data == -75.0] = np.nan
                    #tmp_snow_data[np.isnan(tmp_snow_data)] = 0.0

                    snow_data[sbd][ids,:] = tmp_snow_data*mask

        # if given a list of timesteps
        elif self.timesteps is not None:
            timesteps = np.array(self.timesteps)
            # find indices that have matching times
            itd = []
            for it, t in enumerate(time):
                if t in timesteps:
                    itd.append(it)
            if len(itd) == 0:
                raise ValueError('No matching timesteps to those given')

            itd = np.array(itd)

            # pull data out at specified times
            for ebd in self.embands:
                # find variable name
                e_name = self.emdict[ebd]

                em_data[ebd] = np.zeros((len(itd), ny, nx))
                for ide, it in enumerate(itd):
                    em_data[ebd][ide,:] = ds_em.variables[e_name][it,:]*mask
                    #em_data[ebd][ide,:][em_data[ebd][ide,:] <= -75.0] = np.nan

            for sbd in self.snowbands:
                # find variable name
                s_name = self.snowdict[sbd]
                snow_data[sbd] = np.zeros((len(itd), ny, nx))
                for ids, it in enumerate(itd):
                    snow_data[sbd][ids,:] = ds_snow.variables[s_name][it,:]*mask
                    #snow_data[sbd][ids,:][snow_data[sbd][ids,:] <= -75.0] = np.nan

        # store timesteps
        self.time = timesteps
        # close datasets
        ds_snow.close()
        ds_em.close()

        return snow_data, em_data
