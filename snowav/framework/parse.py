
import logging
import numpy as np
import os
import copy
import snowav.utils.get_topo_stats as ts
from snowav.utils.utilities import masks
from snowav.framework.outputs import outputs
import coloredlogs
import netCDF4 as nc
from collections import OrderedDict
from snowav.database.database import connect
from datetime import timedelta, datetime
from inicheck.output import print_config_report, generate_config
from snowav.utils.wyhr import handle_year_stradling, calculate_date_from_wyhr

def parse(self, external_logger=None):
    '''
    Parse options from config file after read_config().

    '''

    self.barcolors = ['xkcd:cobalt',
                      'xkcd:mustard green',
                      'xkcd:lichen',
                      'xkcd:pale green',
                      'xkcd:blue green',
                      'xkcd:bluish purple',
                      'xkcd:lightish purple',
                      'xkcd:deep magenta']

    # These are used in process() and database table VariableUnits
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

    out = masks(self.dempath, self.db_convert, plotorder = self.plotorder,
                plotlabels = self.plotlabels)

    self.dem = out['dem']
    self.masks = out['masks']
    self.nrows = out['nrows']
    self.ncols = out['ncols']
    self.plotorder = out['plotorder']
    self.labels = out['labels']

    for log in out['logger']:
        self.tmp_log.append(log)

    # Establish database connection
    self.basins, cnx, out = connect(sqlite = self.sqlite, sql = self.mysql,
                               plotorder = self.plotorder, user = self.db_user,
                               password = self.db_password, host = self.db_host,
                               port = self.db_port, convert = self.db_convert,
                               add = self.add_basins)
    self.connector = cnx

    for log in out:
        self.tmp_log.append(log)

    if self.loglevel == 'DEBUG':
        for basin in self.basins:
            self.tmp_log.append(' {}: {}'.format(basin, self.basins[basin]))

        self.tmp_log.append(' Connection: {}'.format(self.connector))

    # Check snow.nc file location, get topo stats and water year
    sfile = os.path.join(self.run_dirs[0],'snow.nc')

    if os.path.isfile(sfile):
        topo = ts.get_topo_stats(sfile, filetype = self.filetype)
        self.snow_x = topo['x']
        self.snow_y = topo['y']
        self.pixel = int(topo['dv'])

        ncf = nc.Dataset(sfile)
        t = nc.num2date(ncf.variables['time'][0],ncf.variables['time'].units)
        ncf.close()
        self.wy = handle_year_stradling(t) + 1

    else:
        print('\nGiven config options, expecting to find:\n   {}\nto load topo '
              'stats but is not a valid file\nCheck config [run] options, see '
              'CoreConfig.ini for details\n'.format(sfile))
        raise Exception('{} not a valid file'.format(sfile))

    # make the bins
    edges = np.arange(self.elev_bins[0],
                      self.elev_bins[1]+self.elev_bins[2],
                      self.elev_bins[2])

    # use for definition
    self.edges = np.arange(self.elev_bins[0]-self.elev_bins[2],
                           self.elev_bins[1],
                           self.elev_bins[2])

    if self.units == 'TAF':
        self.conversion_factor = ((self.pixel**2)*0.000000810713194*0.001)
        self.depth_factor = 0.03937
        self.dem = self.dem * 3.28
        self.depthlbl = 'in'
        self.vollbl = self.units
        self.elevlbl = 'ft'

    self.ixd = np.digitize(self.dem,edges)
    self.xlims = (0,len(edges))

    if self.loglevel == 'DEBUG' and self.log_to_file is not True:
        print('Reading files in {}...'.format(self.run_dirs[0].split('runs')[0]))

    out, all_dirs, dirs, rdict, log = outputs(run_dirs = self.run_dirs,
                                              start_date = self.start_date,
                                              end_date = self.end_date,
                                              filetype = self.filetype,
                                              wy = self.wy,
                                              loglevel = self.loglevel)

    for l in log:
        self.tmp_log.append(l)

    if out['dates'] == []:
        raise Exception('Supplied [run] directory, start_date, and end_date '
              'give no valid snow files')

    self.outputs = out
    self.run_dirs = dirs
    self.all_dirs = all_dirs
    self.rundirs_dict = rdict
    self.all_dirs_flt = copy.deepcopy(all_dirs)

    if self.start_date is not None and self.end_date is None:
        self.end_date = self.outputs['dates'][-1]
        self.tmp_log.append(' Config options [run] end_date '
                            'not specified, assigning '
                            '{} and {}'.format(self.start_date,self.end_date))

        self.ixs = 0
        self.ixe = len(self.outputs['dates']) - 1

    # Otherwise, get closest dates and make indices
    else:
        self.start_date = self.outputs['dates'][0]
        self.end_date = self.outputs['dates'][-1]
        self.ixs = 0
        self.ixe = len(self.outputs['dates']) - 1

    if ((self.start_date.date() < self.outputs['dates'][0].date())
        or (self.end_date.date() > self.outputs['dates'][-1].date())):
        raise Exception('ERROR! Config option [run] start_date or end_date '
                        'outside of date range found in [run] directory')

    # Since model outputs at 23:00, step the figure and report dates to
    # show 00:00 the next day (unless start of water year)
    if self.start_date == datetime(self.wy-1,10,1,23,0,0):
        self.report_start = self.start_date

    else:
        self.report_start = self.start_date + timedelta(hours=1)

    # Copy the config file where figs will be saved
    # use directory if only plotting figures from database and don't
    # have start_date, end_date
    extf = os.path.splitext(os.path.split(self.config_file)[1])
    ext_shr = (self.directory +
              '_'  +
              self.start_date.date().strftime("%Y%m%d") +
              '_' +
              self.end_date.date().strftime("%Y%m%d") )
    self.figs_path = os.path.join(self.save_path, '{}/'.format(ext_shr))

    # get forecast outputs
    if self.forecast_flag:

        out, alldirs, dirs, rdict, log = outputs(run_dirs = self.for_run_dirs,
                                                 start_date = None,
                                                 end_date = None,
                                                 filetype = self.filetype,
                                                 wy = self.wy,
                                                 loglevel = self.loglevel)
        self.for_outputs = out
        self.for_run_dirs = dirs
        self.for_rundirs_dict = rdict
        self.for_ixs = 0
        self.for_ixe = len(self.for_outputs['swe_z']) - 1

    if self.flt_flag:

        file = self.update_file
        p = nc.Dataset(file, 'r')

        if self.update_numbers is None:
            times = p.variables['time'][:]
        else:
            if sum([x > len(p.variables['time']) for x in self.update_numbers]) > 0:
                self.tmp_log.append(' Value in [plots] update_numbers out of '
                                    'range, max is {}, flight update figs '
                                    'being set to False'.format(len(p.variables['time'])))
                times = []
                self.flt_flag = False

            else:
                times = p.variables['time'][self.update_numbers]

        p.close()

        flight_dates = []
        pre_flight_dates = []

        for time in times:
            wydate = calculate_date_from_wyhr(int(time), self.wy)
            pre_wydate = calculate_date_from_wyhr(int(time-24), self.wy)
            flight_dates = np.append(flight_dates,wydate)
            pre_flight_dates = np.append(pre_flight_dates,pre_wydate)

        if self.loglevel == 'DEBUG' and self.log_to_file is not True:
            print('Reading files in {} for flight updates...'.format(self.run_dirs[0].split('runs')[0]))

        out, x, dirs, rdict, log = outputs(run_dirs = self.all_dirs_flt,
                                           start_date = None,
                                           end_date = None,
                                           filetype = self.filetype,
                                           wy = self.wy,
                                           flight_dates = flight_dates,
                                           loglevel = self.loglevel)

        self.flight_outputs = out
        self.run_dirs_flt = dirs
        self.for_rundirs_dict = rdict
        self.flight_diff_dates = out['dates']

        out, x, dirs, rdict, log = outputs(run_dirs = self.all_dirs_flt,
                                           start_date = None,
                                           end_date = None,
                                           filetype = self.filetype,
                                           wy = self.wy,
                                           loglevel = self.loglevel,
                                           flight_dates = pre_flight_dates)

        self.pre_flight_outputs = out

        # If there are no flights in the period, set to false for the flight
        # difference figure and report
        if not self.run_dirs_flt:
            self.flt_flag = False
            self.tmp_log.append(' Config option [plots] update_file was '
                                'supplied, but no snow.nc files were found in '
                                '[run] directory that fit the date range, no '
                                'flight difference figure will be made')

    self.report_date = self.end_date + timedelta(hours=1)
    parts = self.report_name.split('.')
    self.report_name = (parts[0] + self.report_date.date().strftime("%Y%m%d") +
                       '.' + parts[1])

    if not os.path.exists(self.figs_path):
        os.makedirs(self.figs_path)

    config_copy = '{}{}{}'.format(self.figs_path, ext_shr, extf[1])
    generate_config(self.ucfg, config_copy)

    if external_logger == None:
        createLog(self)
    else:
        self._logger = external_logger

    if self.inputs_basins is None:
        self.inputs_basins = [self.plotorder[0]]

    # set up process() inputs for standard run
    self.pargs = {}
    self.pargs['outputs'] = self.outputs
    self.pargs['ixs'] = self.ixs
    self.pargs['ixe'] = self.ixe
    self.pargs['rundirs_dict'] = self.rundirs_dict
    self.pargs['edges'] = self.edges
    self.pargs['masks'] = self.masks
    self.pargs['nrows'] = self.nrows
    self.pargs['ncols'] = self.ncols
    self.pargs['plotorder'] = self.plotorder
    self.pargs['connector'] = self.connector
    self.pargs['basins'] = self.basins
    self.pargs['run_name'] = self.run_name
    self.pargs['db_overwrite'] = self.db_overwrite
    self.pargs['conversion_factor'] = self.conversion_factor
    self.pargs['depth_factor'] = self.depth_factor
    self.pargs['wy'] = self.wy
    self.pargs['vars'] = self.vars
    self.pargs['ixd'] = self.ixd
    self.pargs['vollbl'] = self.vollbl
    self.pargs['depthlbl'] = self.depthlbl
    self.pargs['elevlbl'] = self.elevlbl
    self.pargs['pixel'] = self.pixel
    self.pargs['dem'] = self.dem
    self.pargs['snow_limit'] = self.diag_limit
    self.pargs['inputs_flag'] = self.inputs_flag
    self.pargs['inputs_methods'] = self.inputs_methods
    self.pargs['inputs_basins'] = self.inputs_basins
    self.pargs['inputs_variables'] = self.inputs_variables
    self.pargs['inputs_percentiles'] = self.inputs_percentiles
    if self.mysql is not None:
        self.pargs['dbs'] = 'sql'
    else:
        self.pargs['dbs'] = 'sqlite'

def createLog(self):
    '''
    Create log file and print out saved logging statements.
    '''

    level_styles = {'info': {'color': 'white'},
                    'notice': {'color': 'magenta'},
                    'verbose': {'color': 'blue'},
                    'success': {'color': 'green', 'bold': True},
                    'spam': {'color': 'green', 'faint': True},
                    'critical': {'color': 'red', 'bold': True},
                    'error': {'color': 'red'},
                    'debug': {'color': 'green'},
                    'warning': {'color': 'yellow'}}

    field_styles =  {'hostname': {'color': 'magenta'},
                     'programname': {'color': 'cyan'},
                     'name': {'color': 'white'},
                     'levelname': {'color': 'white', 'bold': True},
                     'asctime': {'color': 'green'}}

    # start logging
    loglevel = self.loglevel
    numeric_level = getattr(logging, loglevel, None)

    if not isinstance(numeric_level, int):
        raise ValueError('Invalid log level: %s' % loglevel)

    # setup the logging
    logfile = None
    if self.log_to_file:
        logfile = os.path.join(self.figs_path, 'log_snowav.out')
        # let user know
        print('Logging to file: {}'.format(logfile))

    fmt = '%(levelname)s:%(module)s:%(message)s'

    if logfile is not None:

        for handler in logging.root.handlers[:]:
            logging.root.removeHandler(handler)

        logging.basicConfig(filename=logfile,
                            filemode='w',
                            level=numeric_level,
                            format=fmt)
    else:
        logging.basicConfig(level=numeric_level)
        coloredlogs.install(level=numeric_level,
                            fmt=fmt,
                            level_styles=level_styles,
                            field_styles=field_styles)

    self._loglevel = numeric_level

    self._logger = logging.getLogger(__name__)

    if len(self.tmp_log) > 0:
        for l in self.tmp_log:
            self._logger.info(l)
    if len(self.tmp_warn) > 0:
        for l in self.tmp_warn:
            self._logger.warning(l)
    if len(self.tmp_err) > 0:
        for l in self.tmp_err:
            self._logger.error(l)
