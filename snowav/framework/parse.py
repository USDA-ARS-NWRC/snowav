
import numpy as np
from shutil import copyfile
import os
import copy
import datetime
import snowav.utils.wyhr_to_datetime as wy
import snowav.utils.get_topo_stats as ts
from snowav.utils.OutputReader import iSnobalReader
import logging
import coloredlogs
import netCDF4 as nc
from dateutil.relativedelta import relativedelta
# from sqlalchemy import create_engine
# from snowav.database.tables import Base, RunMetadata, Watershed, Basin, Results,VariableUnits, Watersheds, Basins
from sqlalchemy.ext.declarative import declarative_base
from sqlalchemy.orm import sessionmaker
from collections import OrderedDict
from snowav import database
from datetime import timedelta
from shutil import copyfile
from sys import exit


def parse(self):
    '''
    Parse options from reading in config file.

    '''

    ####################################################
    #             snowav system                        #
    ####################################################
    if self.save_path is None:
        self.save_path = self.snowav_path + '/snowav/data/'
        self._logger.info('No save_path specified, using ' +
                            '{}'.format(self.save_path))

    ####################################################
    #           outputs                                #
    ####################################################
    if self.flt_start_date is not None:
        self.flt_flag = True

        if not isinstance(self.flt_start_date, datetime.date):
            self.flt_start_date = self.flt_start_date.to_pydatetime()
            self.flt_end_date = self.flt_end_date.to_pydatetime()

    else:
        self.flt_flag = False

    ####################################################
    #         results
    ####################################################
    if type(self.write_csv) != list and self.write_csv != None:
        self.write_csv = [self.write_csv]

    if self.write_csv != None:
        self.write_csv_flag = True
    else:
        self.write_csv_flag = False

    self.barcolors = ['xkcd:cobalt',
                      'xkcd:mustard green',
                      'xkcd:lichen',
                      'xkcd:pale green',
                      'xkcd:blue green',
                      'xkcd:bluish purple',
                      'xkcd:lightish purple',
                      'xkcd:deep magenta']

    # These are used in process() - and order matters!
    # Also used to define database table VariableUnits
    self.vars = OrderedDict([('coldcont','cold content'),
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

    # Establish database connection, create if necessary
    database.database.connect(self)

    self.plotorder = []
    maskpaths = []

    ncf = nc.Dataset(self.dempath, 'r')
    self.dem = ncf.variables['dem'][:]
    self.mask = ncf.variables['mask'][:]
    self.nrows = len(self.dem[:,0])
    self.ncols = len(self.dem[0,:])

    blank = np.zeros((self.nrows,self.ncols))
    self.masks = dict()

    for lbl in self.plotorder:

        # Exceptions
        if lbl == 'Cherry Creek':
            nclbl = 'Cherry'
        else:
            nclbl = lbl

        if lbl != self.plotorder[0]:
            self.masks[lbl] = {'border': blank,
                               'mask': ncf[nclbl + ' mask'][:],
                               'label': lbl}

        else:
            self.masks[lbl] = {'border': blank,
                               'mask': ncf['mask'][:],
                               'label': nclbl}

    ncf.close()

    fp = os.path.join(self.run_dirs[0],'snow.nc')
    topo = ts.get_topo_stats(fp, filetype = self.filetype)
    self.snow_x = topo['x']
    self.snow_y = topo['y']
    self.pixel = int(topo['dv'])
    self.edges = np.arange(self.elev_bins[0],
                           self.elev_bins[1]+self.elev_bins[2],
                           self.elev_bins[2])

    if self.basin == 'LAKES':
        self.imgx = (1250,1475)
        self.imgy = (450,225)

    # Right now this is a placeholder, could edit by basin...
    self.xlims = (0,len(self.edges))

    # Conversion factors and labels
    # Note! If new units are introduced, may need a second look at figure
    # labels, database fields, and writing out variable csv...
    if self.units == 'TAF':
        # KAF, inches, ft
        self.conversion_factor = ((self.pixel**2)
                                 * 0.000000810713194*0.001)
        self.depth_factor = 0.03937
        self.dem = self.dem * 3.28
        self.ixd = np.digitize(self.dem,self.edges)
        self.depthlbl = 'in'
        self.vollbl = self.units
        self.elevlbl = 'ft'

    if self.units == 'SI':
        # m^3, m, mm
        self.conversion_factor = ((self.pixel**2)
                                  * 0.000000810713194*1233.48/1e9)
        self.depth_factor = 1
        self.ixd = np.digitize(self.dem,self.edges)
        self.depthlbl = 'mm'
        self.vollbl = '$km^3$'
        self.elevlbl = 'm'

    # If the user isn't just trying to create figures from the database, load
    # iSnobal outputs and so forth
    if self.plot_flag is False:
        self.outputs = {'swi_z':[], 'evap_z':[], 'snowmelt':[], 'swe_z':[],
                        'depth':[], 'dates':[], 'time':[], 'density':[],
                        'coldcont':[] }

        self.rundirs_dict = {}

        # these were made with an old version of awsm and have the wrong date
        fdirs = ['brb/ops/wy2018/runs/run20171001_20180107/',
                 'brb/ops/wy2018/runs/run20180108_20180117/',
                 'brb/devel/wy2018/hrrr_comparison/run20171001_20180107/',
                 'brb/devel/wy2018/hrrr_comparison/run20180108_20180117/']

        lrdirs = copy.deepcopy(self.run_dirs)
        for rd in lrdirs:

            path = rd

            # If the run_dirs isn't empty use it, otherwise remove
            if (any(os.path.isfile(os.path.join(path, i)) for i in os.listdir(path))
               and not (os.path.isfile(path))):

                d = path.split('runs/run')[-1]
                folder_date = datetime.datetime(int(d[:4]),int(d[4:6]),int(d[6:]))

                # Only load the rundirs that we need
                # If we do this, stn_validate isn't right...
                # if ((self.start_date is not None) and
                #    (folder_date.date() >= self.start_date.date()) and
                #    (folder_date.date() <= self.end_date.date())):
                if (
                    (
                    (self.start_date is not None) and
                    (folder_date.date() <= self.end_date.date()) )
                    or (self.start_date is None)
                    ):

                    output = iSnobalReader(path,
                                           self.filetype,
                                           snowbands = [0,1,2],
                                           embands = [6,7,8,9],
                                           wy = self.wy)

                    if (fdirs[0] in rd) or (fdirs[1] in rd) or (fdirs[2] in rd):
                        self.outputs['dates'] = np.append(
                                self.outputs['dates'],output.dates-relativedelta(years=1) )
                    else:
                        self.outputs['dates'] = np.append(self.outputs['dates'],output.dates)

                    self.outputs['time'] = np.append(self.outputs['time'],output.time)

                    # Make a dict for wyhr-rundir lookup
                    for t in output.time:
                        self.rundirs_dict[int(t)] = rd

                    for n in range(0,len(output.em_data[8])):
                        self.outputs['swi_z'].append(output.em_data[8][n,:,:])
                        self.outputs['snowmelt'].append(output.em_data[7][n,:,:])
                        self.outputs['evap_z'].append(output.em_data[6][n,:,:])
                        self.outputs['coldcont'].append(output.em_data[9][n,:,:])
                        self.outputs['swe_z'].append(output.snow_data[2][n,:,:])
                        self.outputs['depth'].append(output.snow_data[0][n,:,:])
                        self.outputs['density'].append(output.snow_data[1][n,:,:])

                else:
                    self.run_dirs.remove(rd)

            else:
                self.run_dirs.remove(rd)

        # If no dates are specified, use first and last
        if (self.start_date is None) and (self.end_date is None):
            self.start_date = self.outputs['dates'][0]
            self.end_date = self.outputs['dates'][-1]

            self._logger.info('start_date and/or end_date not specified, '
                              'using: {} and {}'.format(self.start_date,
                                                        self.end_date))
            self.ixs = 0
            self.ixe = len(self.outputs['dates']) - 1

        # Otherwise, get closest dates and make indices
        else:
            s = min(self.outputs['dates'],key=lambda x: abs(x-self.start_date))
            e = min(self.outputs['dates'],key=lambda x: abs(x-self.end_date))
            self.ixs = np.where(self.outputs['dates'] == s)[0][0]
            self.ixe = np.where(self.outputs['dates'] == e)[0][0]

        if ((self.start_date.date() < self.outputs['dates'][0].date())
            or (self.end_date.date() > self.outputs['dates'][-1].date())):
            self._logger.info('ERROR! [Outputs] -> start_date or end_date ' +
                              'outside of range in [runs] -> run_dirs')
            # print(self.start_date.date(), self.outputs['dates'][0].date(),
            # self.end_date.date() , self.outputs['dates'][-1].date())
            exit()

        # Since model outputs at 23:00, step the figure and report dates to
        # show 00:00 the next day (unless start of water year)
        if self.start_date == datetime.datetime(self.wy-1,10,1,23,0,0):
            self.report_start = self.start_date
        else:
            self.report_start = self.start_date + timedelta(hours=1)

        # Copy the config file where figs will be saved
        # use name_append if only plotting figures from database and don't
        # have start_date, end_date
        extf = os.path.splitext(os.path.split(self.config_file)[1])
        ext_shr = (self.name_append +
                  '_'  +
                  self.start_date.date().strftime("%Y%m%d") +
                  '_' +
                  self.end_date.date().strftime("%Y%m%d") )
        self.figs_path = os.path.join(self.save_path, '%s/'%(ext_shr))

    # Otherwise, all we need to do is create the figs_path
    else:
        extf = os.path.splitext(os.path.split(self.config_file)[1])
        # ext_shr = self.name_append + '_'  + self.start_date.date().strftime("%Y%m%d") + '_' + self.end_date.date().strftime("%Y%m%d")
        ext_shr = self.name_append
        self.figs_path = os.path.join(self.save_path, '%s/'%(ext_shr))


    self.report_date = self.end_date + timedelta(hours=1)
    parts = self.report_name.split('.')
    self.report_name = ( parts[0]
                       + self.report_date.date().strftime("%Y%m%d")
                       + '.' + parts[1] )

    if not os.path.exists(self.figs_path):
        os.makedirs(self.figs_path)
