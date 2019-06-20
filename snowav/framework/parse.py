
import logging
import numpy as np
import os
import copy
import snowav.utils.get_topo_stats as ts
from snowav.utils.snowav_masks import snowav_masks
from snowav.framework.outputs import outputs
import coloredlogs
import netCDF4 as nc
from collections import OrderedDict
from snowav.database.database import connect, make_session
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

    out = snowav_masks(self.dempath, plotorder = self.plotorder,
                       plotlabels = self.plotlabels, logger = self.tmp_log)

    self.dem = out['dem']
    self.masks = out['masks']
    self.nrows = out['nrows']
    self.ncols = out['ncols']
    self.plotorder = out['plotorder']
    self.labels = out['labels']
    self.tmp_log = out['logger']

    # Establish database connection.
    # This will make sqlite database if it doesn't exist
    basins, cnx, self.tmp_log = connect(sqlite = self.sqlite, mysql = self.mysql,
                               plotorder = self.plotorder, user = self.db_user,
                               password = self.db_password, host = self.db_host,
                               port = self.db_port, logger = self.tmp_log)
    self.connector = cnx

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
        raise Exception('Expecting to find {} to get '.format(sfile) +
                        'topo stats, but it is not a valid file')

    # make the bins
    edges = np.arange(self.elev_bins[0],
                      self.elev_bins[1]+self.elev_bins[2],
                      self.elev_bins[2])

    # use for definition
    self.edges = np.arange(self.elev_bins[0]-self.elev_bins[2],
                           self.elev_bins[1],
                           self.elev_bins[2])

    self.xlims = (0,len(edges))

    # Conversion factors and labels
    # Note! If new units are introduced, may need a second look at figure
    # labels, database fields, and writing out variable csv...
    if self.units == 'TAF':
        self.conversion_factor = ((self.pixel**2)*0.000000810713194*0.001)
        self.depth_factor = 0.03937
        self.dem = self.dem * 3.28
        self.ixd = np.digitize(self.dem,edges)
        self.depthlbl = 'in'
        self.vollbl = self.units
        self.elevlbl = 'ft'

    if self.units == 'SI':
        self.conversion_factor = ((self.pixel**2)*0.000000810713194*1233.48/1e9)
        self.depth_factor = 1
        self.ixd = np.digitize(self.dem,edges)
        self.depthlbl = 'mm'
        self.vollbl = '$km^3$'
        self.elevlbl = 'm'

    # get standard run outputs
    out, dirs, rdict = outputs(run_dirs = self.run_dirs,
                               start_date = self.start_date,
                               end_date = self.end_date,
                               filetype = self.filetype,
                               wy = self.wy)
    self.outputs = out
    self.run_dirs = dirs
    self.rundirs_dict = rdict

    # If no dates are specified, use first and last
    if (self.start_date is None) and (self.end_date is None):
        self.start_date = self.outputs['dates'][0]
        self.end_date = self.outputs['dates'][-1]

        self.tmp_log.append(' Config options [run] start_date and/or end_date '
                            'not specified, assigning '
                            '{} and {}'.format(self.start_date,self.end_date))
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
        raise Exception('ERROR! Config option [run] start_date or end_date '
                        'outside of range in run_dirs')

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
    self.figs_path = os.path.join(self.save_path, '%s/'%(ext_shr))

    extf = os.path.splitext(os.path.split(self.config_file)[1])
    ext_shr = self.directory
    self.figs_path = os.path.join(self.save_path, '%s/'%(ext_shr))

    # get forecast outputs
    if self.forecast_flag:

        out, dirs, rdict = outputs(run_dirs = self.for_run_dirs,
                                   start_date = None,
                                   end_date = None,
                                   filetype = self.filetype,
                                   wy = self.wy)
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
            times = p.variables['time'][self.update_numbers]
        p.close()

        flight_dates = []
        pre_flight_dates = []

        for time in times:
            wydate = calculate_date_from_wyhr(int(time), self.wy)
            pre_wydate = calculate_date_from_wyhr(int(time-24), self.wy)
            flight_dates = np.append(flight_dates,wydate)
            pre_flight_dates = np.append(pre_flight_dates,wydate)

        self.flight_outputs = {'swe_z':[],'depth':[], 'dates':[],
                               'time':[], 'density':[]}
        self.pre_flight_outputs = {'swe_z':[],'depth':[], 'dates':[],
                               'time':[], 'density':[]}

        self.flight_rundirs_dict = {}
        flist = copy.deepcopy(self.run_dirs_flight)

        for rd in flist:
            path = rd

            if (any(os.path.isfile(os.path.join(path, i)) for i in os.listdir(path))
               and not (os.path.isfile(path))):

                d = path.split('runs/run')[-1]
                folder_date = datetime(int(d[:4]),int(d[4:6]),int(d[6:8]))

                if ((folder_date.date() in [x.date() for x in flight_dates]) and
                    (folder_date.date() <= self.end_date.date())):

                    output = iSnobalReader(path,
                                           self.filetype,
                                           snowbands = [0,1,2],
                                           wy = self.wy)

                    self.flight_outputs['dates'] = np.append(self.flight_outputs['dates'],output.dates)
                    self.flight_outputs['time'] = np.append(self.flight_outputs['time'],output.time)

                    for t in output.time:
                        self.flight_rundirs_dict[int(t)] = rd

                    for n in range(0,len(output.snow_data[0])):
                        self.flight_outputs['swe_z'].append(output.snow_data[2][n,:,:])
                        self.flight_outputs['depth'].append(output.snow_data[0][n,:,:])
                        self.flight_outputs['density'].append(output.snow_data[1][n,:,:])

                # Get the previous snow.nc as necessary
                if ((folder_date.date() in [x.date()-timedelta(hours = 24) for x in pre_flight_dates]) and
                    (folder_date.date() <= self.end_date.date())):

                    output = iSnobalReader(path,
                                           self.filetype,
                                           snowbands = [0,1,2],
                                           wy = self.wy)
                    self.pre_flight_outputs['dates'] = np.append(self.pre_flight_outputs['dates'],output.dates)
                    self.pre_flight_outputs['time'] = np.append(self.pre_flight_outputs['time'],output.time)

                    for n in range(0,len(output.snow_data[0])):
                        self.pre_flight_outputs['swe_z'].append(output.snow_data[2][n,:,:])
                        self.pre_flight_outputs['depth'].append(output.snow_data[0][n,:,:])
                        self.pre_flight_outputs['density'].append(output.snow_data[1][n,:,:])

                else:
                    self.run_dirs_flight.remove(rd)

            else:
                self.run_dirs_flight.remove(rd)

    # If there are no flights in the period, set to false for the flight diff
    # figure and report
    if not self.run_dirs_flight:
        self.flt_flag = False

    #############################
    # ^ FLIGHT section done ^ #
    #############################

    self.report_date = self.end_date + timedelta(hours=1)
    parts = self.report_name.split('.')
    self.report_name = (parts[0] + self.report_date.date().strftime("%Y%m%d") +
                       '.' + parts[1])

    if not os.path.exists(self.figs_path):
        os.makedirs(self.figs_path)

    config_copy = '{}{}{}'.format(self.figs_path, ext_shr, extf[1])
    generate_config(self.ucfg, config_copy)

    ####################################################
    #             log file                             #
    ####################################################
    if external_logger == None:
        createLog(self)
    else:
        self._logger = external_logger


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

    fmt = '%(levelname)s:%(name)s:%(message)s'

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
