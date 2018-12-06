
import numpy as np
from smrf import ipw
from shutil import copyfile
import os
import copy
import pandas as pd
import datetime
import snowav.utils.wyhr_to_datetime as wy
import snowav.utils.get_topo_stats as ts
from snowav.utils.utilities import get_snowav_path
from snowav.utils.OutputReader import iSnobalReader
from inicheck.tools import get_user_config, check_config
from inicheck.output import generate_config, print_config_report
from inicheck.config import MasterConfig
import logging
import coloredlogs
import netCDF4 as nc
from dateutil.relativedelta import relativedelta
from sqlalchemy import create_engine
from snowav.database.tables import Base, RunMetadata, Watershed, Basin, Results,VariableUnits, Watersheds, Basins
from sqlalchemy.ext.declarative import declarative_base
from sqlalchemy.orm import sessionmaker
from collections import OrderedDict
from snowav import database
from datetime import timedelta


def read_config(self, external_logger=None, awsm=None):
    '''
    Read snowav config file and assign fields.

    Args
        external_logger: awsm logger
        awsm: awsm class, if this is passed, run_dir will be assigned from the
            directory being created in awsm

    '''
    print('Reading SNOWAV config file...')
    ucfg = get_user_config(self.config_file, modules = 'snowav')

    # find path to snowav code directory
    self.snowav_path = get_snowav_path()

    # create blank log and error log because logger is not initialized yet
    self.tmp_log = []
    self.tmp_err = []
    self.tmp_warn = []

    # Check the user config file for errors and report issues if any
    self.tmp_log.append("Reading config file and loading iSnobal outputs...")
    warnings, errors = check_config(ucfg)

    ####################################################
    #             snowav system                        #
    ####################################################
    self.loglevel = ucfg.cfg['snowav system']['log_level'].upper()
    self.log_to_file = ucfg.cfg['snowav system']['log_to_file']
    self.basin = ucfg.cfg['snowav system']['basin']
    self.save_path = ucfg.cfg['snowav system']['save_path']

    if self.save_path is None:
        self.save_path = self.snowav_path + '/snowav/data/'
        print('No save_path specified, using {}'.format(self.save_path))

    self.wy = ucfg.cfg['snowav system']['wy']
    self.units = ucfg.cfg['snowav system']['units']
    self.filetype = ucfg.cfg['snowav system']['filetype']
    self.elev_bins = ucfg.cfg['snowav system']['elev_bins']
    if ucfg.cfg['snowav system']['name_append'] != None:
        self.name_append = ucfg.cfg['snowav system']['name_append']
    else:
        self.name_append = '_gen_' + \
                           datetime.datetime.now().strftime("%Y-%-m-%-d")

    ####################################################
    #           outputs                                #
    ####################################################
    self.dplcs = ucfg.cfg['outputs']['decimals']
    self.start_date = ucfg.cfg['outputs']['start_date']
    self.end_date = ucfg.cfg['outputs']['end_date']

    if (self.start_date is not None and self.end_date is not None):
        self.start_date = self.start_date.to_pydatetime()
        self.end_date = self.end_date.to_pydatetime()

        if self.start_date >= self.end_date:
            print('WARNING! start_date > end_date, needs to be fixed in '
                  'config file, exiting...')
            return

    # Check for forced flight comparison images
    self.flt_start_date = ucfg.cfg['outputs']['flt_start_date']
    self.flt_end_date = ucfg.cfg['outputs']['flt_end_date']

    if self.flt_start_date is not None:
        self.flt_flag = True

        if not isinstance(self.flt_start_date, datetime.date):
            self.flt_start_date = self.flt_start_date.to_pydatetime()
            self.flt_end_date = self.flt_end_date.to_pydatetime()

    else:
        self.flt_flag = False

    # Once we load in self.outputs, will make the indices for these dates the
    # closest to a specified hour

    self.summary = ucfg.cfg['outputs']['summary']
    if type(self.summary) != list:
        self.summary = [self.summary]

    ####################################################
    #           runs                                   #
    ####################################################
    # If True, run_dirs are all sub directories
    self.all_subdirs = ucfg.cfg['runs']['all_subdirs']

    if self.all_subdirs is True:
        self.run_dirs = ([ucfg.cfg['runs']['run_dirs'] + s for s in
                        os.listdir(ucfg.cfg['runs']['run_dirs'])
                        if (os.path.isdir(ucfg.cfg['runs']['run_dirs'] + s)) ])
    else:
        self.run_dirs = ucfg.cfg['runs']['run_dirs']
        if type(self.run_dirs) != list:
            self.run_dirs = [self.run_dirs]

    self.run_dirs.sort()

    ####################################################
    #           validate                               #
    ####################################################
    if (ucfg.cfg['validate']['stations'] != None and
        ucfg.cfg['validate']['labels'] != None and
        ucfg.cfg['validate']['client'] != None):

        self.val_stns = ucfg.cfg['validate']['stations']
        self.val_lbls = ucfg.cfg['validate']['labels']
        self.val_client = ucfg.cfg['validate']['client']
        self.valid_flag = True

    else:
        self.valid_flag = False

    self.pre_val_stns = ucfg.cfg['validate']['pre_stations']
    self.pre_val_lbls = ucfg.cfg['validate']['pre_labels']
    self.precip_figure = ucfg.cfg['validate']['precip_figure']

    # This is being used to combine 2017 HRRR data
    self.offset = int(ucfg.cfg['validate']['offset'])

    ####################################################
    #           basin total                            #
    ####################################################
    if ucfg.cfg['basin total']['summary_swe'] != None:
        self.summary_swe = ucfg.cfg['basin total']['summary_swe']
        self.summary_swi = ucfg.cfg['basin total']['summary_swi']
        self.basin_total_flag = True
    else:
        self.basin_total_flag = False

    if ucfg.cfg['basin total']['netcdf']:
        self.ncvars = ucfg.cfg['Basin Total']['netcdf'].split(',')
        self.nc_flag = True
    else:
        self.nc_flag = False

    self.flight_dates = ucfg.cfg['basin total']['flights']
    if (self.flight_dates is not None) and (type(self.flight_dates) != list):
        self.flight_dates = [self.flight_dates]

    ####################################################
    #           masks                                  #
    ####################################################
    self.dempath = ucfg.cfg['masks']['dempath']

    if ucfg.cfg['masks']['basin_masks'] is not None:
        self.total = ucfg.cfg['masks']['basin_masks'][0]

    ####################################################
    #          plots                                   #
    ####################################################
    self.figsize = (ucfg.cfg['plots']['fig_length'],
                    ucfg.cfg['plots']['fig_height'])
    self.dpi = ucfg.cfg['plots']['dpi']

    self.barcolors = ['xkcd:cobalt',
                      'xkcd:mustard green',
                      'xkcd:lichen', # light grey green
                      'xkcd:pale green',
                      'xkcd:blue green',
                      'xkcd:bluish purple',
                      'xkcd:lightish purple',
                      'xkcd:deep magenta']

    self.annot_x = ucfg.cfg['plots']['annot_x']
    self.annot_y = ucfg.cfg['plots']['annot_y']

    # San SanJoaquin
    # Jose Creek, South Fork, Main, Willow Creek
    # annot_x: 950,1600,375,10
    # annot_y: 1450,500,185,600

    # Tuolumne
    # Tuolumne, Cherry Creek, Eleanor
    # annot_x: 875,175,50
    # annot_y: 375,185,900

    ####################################################
    #          report                                  #
    ####################################################
    self.report_flag = ucfg.cfg['report']['report']
    self.exclude_figs = ucfg.cfg['report']['exclude_figs']

    if type(self.exclude_figs) != list and self.exclude_figs != None:
        self.exclude_figs = [self.exclude_figs]

    if (self.valid_flag is False) and (self.exclude_figs is not None):
        self.exclude_figs.append('VALID')
    if (self.valid_flag is False) and (self.exclude_figs is None):
        self.exclude_figs = ['VALID']

    self.report_name = ucfg.cfg['report']['report_name']
    self.rep_title = ucfg.cfg['report']['report_title']
    self.rep_path = ucfg.cfg['report']['rep_path']
    self.env_path = ucfg.cfg['report']['env_path']
    self.templ_path = ucfg.cfg['report']['templ_path']
    self.tex_file = ucfg.cfg['report']['tex_file']
    self.summary_file = ucfg.cfg['report']['summary_file']
    self.figs_tpl_path = ucfg.cfg['report']['figs_tpl_path']

    # check paths to see if they need default snowav path
    if self.rep_path is None:
        self.rep_path = os.path.join(self.snowav_path,'snowav/data/')
    if self.env_path is None:
        self.env_path = os.path.join(self.snowav_path,
                                     'snowav/report/template/section_text/')
    if self.templ_path is None:
        self.templ_path = os.path.join(self.snowav_path,
                                       'snowav/report/template/')
    if self.summary_file is None:
        self.summary_file = os.path.join(self.snowav_path,
                      'snowav/report/template/section_text/report_summary.txt')
    if self.tex_file is None:
        self.tex_file = os.path.join(self.snowav_path,
                                     'snowav/report/template/snowav_report.tex')
    if self.figs_tpl_path is None:
        self.figs_tpl_path = os.path.join(self.snowav_path,
                                          'snowav/report/figs/')

    ####################################################
    #           hx forecast
    ####################################################
    self.adj_hours = ucfg.cfg['hx forecast']['adj_hours']

    ####################################################
    #           masks
    ####################################################
    for item in ['basin_masks', 'mask_labels']:
        if type(ucfg.cfg['masks'][item]) != list:
            ucfg.cfg['masks'][item] = [ucfg.cfg['masks'][item]]

    masks = ucfg.cfg['masks']['basin_masks']

    ####################################################
    #         results
    ####################################################
    self.db_user = ucfg.cfg['results']['user']
    self.db_password = ucfg.cfg['results']['password']
    self.db_host = ucfg.cfg['results']['host']
    self.db_port = ucfg.cfg['results']['port']
    self.mysql = ucfg.cfg['results']['mysql']
    self.sqlite = ucfg.cfg['results']['sqlite']
    self.write_csv = ucfg.cfg['results']['write_csv']

    if type(self.write_csv) != list and self.write_csv != None:
        self.write_csv = [self.write_csv]

    if self.write_csv != None:
        self.write_csv_flag = True
    else:
        self.write_csv_flag = False

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

    self.run_name = ucfg.cfg['results']['run_name']
    self.write_db = ucfg.cfg['results']['write_db']
    self.plot_runs = ucfg.cfg['results']['plot_runs']
    self.figures_only = ucfg.cfg['results']['figures_only']

    # Establish database connection, create if necessary
    database.database.connect(self)

    # If these fields are specified in the config_file, plot_flag will be set
    # to False, no processing occurs, and no report is generated
    if self.figures_only is True:
        self.plot_flag = True
        self.plot_labels = ucfg.cfg['results']['plot_labels']
        self.plot_wy = ucfg.cfg['results']['plot_wy']
        self.plot_variables = ucfg.cfg['results']['plot_variables']

        if ( (self.plot_labels is None) or
           (self.plot_wy is None) or
           (self.plot_variables is None) ):
            print('Must fill config file options [results] -> plot_labels, '
                 'plot_wy, and plot_variables!')

        print('Config file option [results]->figures_only={}, '
              'creating figures for {}...'.format(self.plot_flag, self.figures_only))
    else:
        self.plot_flag = False

    ########################################################################
    #               done reading config options...                         #
    ########################################################################
    self.plotorder = []
    maskpaths = []

    # ascii/nc
    if (os.path.splitext(self.dempath)[1] != '.nc'):
        for idx, m in enumerate(masks):
            maskpaths.append(m)
            self.plotorder.append(ucfg.cfg['masks']['mask_labels'][idx])

        try:
            self.dem = np.genfromtxt(self.dempath)
            sr = 0
        except:
            self.dem = np.genfromtxt(self.dempath,skip_header = 6)
            sr = 6

        self.nrows = len(self.dem[:,0])
        self.ncols = len(self.dem[0,:])
        blank = np.zeros((self.nrows,self.ncols))

        self.masks = dict()
        for lbl,mask in zip(self.plotorder,maskpaths):
            self.masks[lbl] = {'border': blank,
                               'mask': np.genfromtxt(mask,skip_header=sr),
                               'label': lbl}

    if os.path.splitext(self.dempath)[1] == '.nc':
        self.plotorder = ucfg.cfg['masks']['mask_labels']
        ncf = nc.Dataset(self.dempath, 'r')

        self.dem = ncf.variables['dem'][:]
        self.mask = ncf.variables['mask'][:]

        self.nrows = len(self.dem[:,0])
        self.ncols = len(self.dem[0,:])
        blank = np.zeros((self.nrows,self.ncols))

        self.masks = dict()
        for iters, lbl in enumerate(self.plotorder):
            if iters > 0:
                self.masks[lbl] = {'border': blank,
                                   'mask': ncf[lbl + ' mask'][:],
                                   'label': lbl}
            else:
                self.masks[lbl] = {'border': blank,
                                   'mask': ncf['mask'][:],
                                   'label': lbl}

        ncf.close()

    # Pixel size and elevation bins
    if self.filetype == 'netcdf':
        # split on the last instance of 'output', not very robust
        if 'output' in self.run_dirs[0]:
            fp = os.path.abspath(self.run_dirs[0].rsplit('output/',1)[0] + 'snow.nc')
        else:
            fp = os.path.join(self.run_dirs[0],'snow.nc')

    if self.filetype == 'ipw':
        fp = os.path.join(self.run_dirs[0], os.listdir(self.run_dirs[0])[0])

    topo = ts.get_topo_stats(fp, filetype = self.filetype)

    self.pixel = int(topo['dv'])
    self.edges = np.arange(self.elev_bins[0],
                           self.elev_bins[1]+self.elev_bins[2],
                           self.elev_bins[2])

    # A few remaining basin-specific things
    if self.basin in ['TUOL', 'SJ', 'LAKES']:
        sr = 6
    else:
        sr = 0

    if self.basin == 'LAKES':
        self.imgx = (1200,1375)
        self.imgy = (425,225)

    # Right now this is a placeholder, could edit by basin...
    self.xlims = (0,len(self.edges))

    # Conversion factors and labels
    # Note! If new units are introduced, may need a second look at figure
    # labels, database fields, and writing out variable csv...
    if self.units == 'KAF':
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

        for rd in self.run_dirs:
            if self.filetype == 'ipw':
                path = rd

            if self.filetype == 'netcdf':
                # split on the last instance of 'output', not very robust
                if 'output' in rd:
                    path = rd.rsplit('output',1)[0]
                else:
                    path = rd

            # If the run_dirs isn't empty use it, otherwise remove
            if (any(os.path.isfile(os.path.join(path, i)) for i in os.listdir(path))
               and not (os.path.isfile(path))):
                try:
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
                except:
                    print('Failure loading file in {}, attempting to skip...'.format(path))

            else:
                self.run_dirs.remove(rd)

        # If no dates are specified, use first and last
        if (self.start_date is None) and (self.end_date is None):
            self.start_date = self.outputs['dates'][0]
            self.end_date = self.outputs['dates'][-1]

            print('start_date and/or end_date not specified, '
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
            print('ERROR! [Outputs] -> start_date or end_date outside of '
                  'range in [runs] -> run_dirs')
            print( self.outputs['dates'][0].date(),self.outputs['dates'][-1].date())
            return

        self.report_date = self.end_date + timedelta(hours=1)
        parts = self.report_name.split('.')
        self.report_name = ( parts[0]
                           + self.report_date.date().strftime("%Y%m%d")
                           + '.' + parts[1] )

        # Copy the config file where figs will be saved
        # use name_append if only plotting figures from database and don't
        # have start_date, end_date
        extf = os.path.splitext(os.path.split(self.config_file)[1])
        ext_shr = self.name_append + '_'  + self.start_date.date().strftime("%Y%m%d") + '_' + self.end_date.date().strftime("%Y%m%d")
        self.figs_path = os.path.join(self.save_path, '%s/'%(ext_shr))

    # Otherwise, all we need to do is create the figs_path
    else:
        extf = os.path.splitext(os.path.split(self.config_file)[1])
        # ext_shr = self.name_append + '_'  + self.start_date.date().strftime("%Y%m%d") + '_' + self.end_date.date().strftime("%Y%m%d")
        ext_shr = self.name_append
        self.figs_path = os.path.join(self.save_path, '%s/'%(ext_shr))

    if not os.path.exists(self.figs_path):
        os.makedirs(self.figs_path)


    ####################################################
    #             log file                             #
    ####################################################
    if external_logger == None:
        createLog(self)
    else:
        self._logger = external_logger

    # Only need to store this name if we decide to
    # write more to the copied config file...
    self.config_copy = (self.figs_path + extf[0] + self.name_append
                        + '_%s'%(ext_shr) + extf[1])

    if not os.path.isfile(self.config_copy):
        generate_config(ucfg,self.config_copy)

def createLog(self):
    '''
    Now that the directory structure is done, create log file and print out
    saved logging statements.
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
