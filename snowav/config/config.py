import calendar
import coloredlogs
from copy import deepcopy
from datetime import timedelta, datetime
import logging
import numpy as np
import os
import netCDF4 as nc

from inicheck.tools import get_user_config, check_config, cast_all_variables
from inicheck.output import generate_config, print_config_report
from inicheck.config import MasterConfig
import snowav
from snowav.utils.utilities import masks, get_snowav_path
from snowav.utils.get_topo_stats import get_topo_stats
from snowav.framework.outputs import outputs
from snowav.database.database import connect
from snowav.database.models import AwsmInputsOutputs
from snowav.utils.wyhr import handle_year_stradling, calculate_date_from_wyhr


class UserConfig(object):
    """ Config file class.

    Args
    ------
    config_file: file path
    awsm: awsm object
    end_date: string
    """

    def __init__(self, config_file, awsm=None, end_date=None):

        print('Reading {} and loading files...'.format(config_file))

        self.config_file = config_file
        snowav_mcfg = MasterConfig(modules='snowav')
        ucfg = get_user_config(self.config_file, mcfg=snowav_mcfg)
        ucfg.apply_recipes()
        ucfg = cast_all_variables(ucfg, ucfg.mcfg)
        self.snowav_path = get_snowav_path()

        warnings, errors = check_config(ucfg)
        if errors != [] or warnings != []:
            print_config_report(warnings, errors)

        self.tmp_log = []
        self.tmp_err = []
        self.tmp_warn = []
        self.proc_time_start = datetime.now()

        ####################################################
        #            snowav                                #
        ####################################################
        self.loglevel = ucfg.cfg['snowav']['log_level'].upper()
        self.log_to_file = ucfg.cfg['snowav']['log_to_file']
        self.save_path = ucfg.cfg['snowav']['save_path']
        self.units = ucfg.cfg['snowav']['units']
        self.elev_bins = ucfg.cfg['snowav']['elev_bins']
        self.directory = ucfg.cfg['snowav']['directory']
        self.dempath = ucfg.cfg['snowav']['dempath']
        self.run_name = ucfg.cfg['snowav']['run_name']
        self.plotorder = ucfg.cfg['snowav']['masks']
        self.plotlabels = ucfg.cfg['snowav']['plotlabels']
        self.report_only = ucfg.cfg['snowav']['report_only']

        ####################################################
        #           run                                    #
        ####################################################
        self.dplcs = ucfg.cfg['run']['decimals']
        self.start_date = ucfg.cfg['run']['start_date']
        self.end_date = ucfg.cfg['run']['end_date']

        if end_date is not None:
            self.end_date = end_date
            self.tmp_log.append(' Overriding config end_date with '
                                '{} given with snowav call'.format(end_date))

            if self.end_date <= self.start_date:
                raise Exception('end_date {} earlier than start_date {}'.format(
                    self.end_date, self.start_date))

        if self.start_date is not None and self.end_date is not None:
            self.start_date = self.start_date
            self.end_date = self.end_date

            if self.start_date >= self.end_date:
                self.tmp_log.append(' Error: [run] start_date >= end_date')
                raise Exception('[run] start_date >= [run] end_date')

        else:
            self.tmp_log.append(' [run] start_date and/or end_date was not '
                                'defined in config file, will be assigned '
                                'by available dates in directory')

        self.all_subdirs = ucfg.cfg['run']['all_subdirs']

        if (ucfg.cfg['run']['directory'] is None) and (awsm is not None):
            if self.all_subdirs is True:
                self.run_dirs = ([awsm.pathr + s for s in
                                  os.listdir(awsm.pathr)
                                  if (os.path.isdir(awsm.pathr + s))])
            else:
                self.run_dirs = awsm.pathr
                if type(self.run_dirs) != list:
                    self.run_dirs = [self.run_dirs]

        else:
            directory = ucfg.cfg['run']['directory']

            if len(directory) == 1:
                directory = directory[0]

            if self.all_subdirs is True:
                self.run_dirs = ([directory + s for s in os.listdir(directory)
                                  if (os.path.isdir(directory + s))])
            else:
                self.run_dirs = ucfg.cfg['run']['directory']
                if type(self.run_dirs) != list:
                    self.run_dirs = [self.run_dirs]

        self.run_dirs.sort()

        ####################################################
        #         database
        ####################################################
        self.mysql = ucfg.cfg['database']['mysql']
        self.db_user = ucfg.cfg['database']['user']
        self.db_password = ucfg.cfg['database']['password']
        self.db_host = ucfg.cfg['database']['host']
        self.db_port = ucfg.cfg['database']['port']
        self.db_convert = ucfg.cfg['database']['convert_ws']
        self.add_basins = ucfg.cfg['database']['add_basins']
        self.db_overwrite = ucfg.cfg['database']['overwrite']
        self.properties = ucfg.cfg['database']['properties']
        self.sqlite = ucfg.cfg['database']['sqlite']

        base_bands = ['swi_z', 'evap_z', 'swe_z', 'depth', 'density',
                      'coldcont', 'precip_z']

        for band in base_bands:
            if band not in self.properties:
                self.tmp_log.append(' WARNING! Config option [database] '
                                    'properties does not contain '
                                    '{}'.format(band))

        if ((self.mysql is not None) and
                ((self.db_user is None) or
                 (self.db_password is None) or
                 (self.db_host is None) or
                 (self.db_port is None))):
            raise Exception('If using config option [database] mysql, must '
                            'also supply user, password, host, and port')

        if self.sqlite is not None:
            if not os.path.isdir(os.path.dirname(self.sqlite)):
                raise Exception('{} does not contain a valid base '
                                'path'.format(self.sqlite))
            self.sqlite = 'sqlite:///' + self.sqlite

            if self.mysql is not None:
                raise Exception('Config option [database] section contains '
                                'both "mysql" and "sqlite" entries, pick one.')

        ####################################################
        #           validate                               #
        ####################################################
        self.val_stns = ucfg.cfg['validate']['stations']
        self.val_lbls = ucfg.cfg['validate']['labels']
        self.val_client = ucfg.cfg['validate']['client']
        self.pre_val_stns = ucfg.cfg['validate']['pre_stations']
        self.pre_val_lbls = ucfg.cfg['validate']['pre_labels']
        self.wxdb_user = ucfg.cfg['validate']['user']
        self.wxdb_password = ucfg.cfg['validate']['password']
        self.wxdb_host = ucfg.cfg['validate']['host']
        self.wxdb_port = ucfg.cfg['validate']['port']
        self.point_values = ucfg.cfg['validate']['point_values']
        self.point_values_csv = ucfg.cfg['validate']['point_values_csv']
        self.point_values_date = ucfg.cfg['validate']['point_values_date']
        self.point_values_properties = ucfg.cfg['validate']['point_values_properties']
        self.point_values_heading = ucfg.cfg['validate']['point_values_heading']
        self.point_values_settings = ucfg.cfg['validate']['point_values_settings']

        for n in range(0, 10):
            self.point_values_settings[n] = int(self.point_values_settings[n])

        if self.point_values and self.point_values_csv is None:
            self.point_values = False
            self.tmp_log.append(' Config option [validate] point_values_csv '
                                'was not supplied, point_values being set '
                                'to False')

        if self.point_values and self.point_values_date is None:
            self.point_values = False
            self.tmp_log.append(' Config option [validate] point_values_date '
                                'was not supplied, point_values being set '
                                'to False')

        ####################################################
        #           diagnostics                            #
        ####################################################
        self.diagnostics_flag = ucfg.cfg['diagnostics']['diagnostics']
        self.diag_basins = ucfg.cfg['diagnostics']['basins']
        self.diag_limit = ucfg.cfg['diagnostics']['limit']

        if self.diagnostics_flag:
            if self.diag_basins is not None:
                for basin in self.diag_basins:
                    if basin not in self.plotorder:
                        self.tmp_log.append(' Config [diagnostics] basin: "{}"'
                                            ' does not match [snowav] masks: '
                                            '"{}", diagnostics set to '
                                            'False'.format(basin,
                                                           self.plotorder))
                        self.diagnostics_flag = False

        if 'snow_line' not in self.properties and self.diagnostics_flag:
            self.diagnostics_flag = False
            self.tmp_log.append(' Required properties in [database] properties'
                                ' for [diagnostics] does not exist, setting '
                                'diagnostics: False')

        self.inputs_flag = ucfg.cfg['diagnostics']['inputs_table']
        self.inputs_variables = ucfg.cfg['diagnostics']['inputs_variables']
        self.inputs_percentiles = ucfg.cfg['diagnostics']['inputs_percentiles']
        self.inputs_methods = ucfg.cfg['diagnostics']['inputs_methods']
        self.inputs_basins = ucfg.cfg['diagnostics']['inputs_basins']

        if self.inputs_basins is not None and self.plotorder is not None:
            for basin in self.inputs_basins:
                if basin not in self.plotorder:
                    self.tmp_log.append(' Config option [diagnostics] '
                                        'inputs_basins: {} does not match what '
                                        'was supplied in [snowav] masks: {}, '
                                        'inputs set to '
                                        'False'.format(basin, self.plotorder))
                    self.inputs_flag = False

        if self.inputs_flag:
            s = [x + ', ' for x in self.inputs_variables]
            self.tmp_log.append(' Using variables {} for inputs '
                                'summary'.format(''.join(s)))

            s = [x + ', ' for x in self.inputs_methods]
            self.tmp_log.append(' Using methods {} for inputs '
                                'summary'.format(''.join(s)))

        ####################################################
        #          plots                                   #
        ####################################################
        self.dpi = ucfg.cfg['plots']['dpi']
        self.depth_clip = ucfg.cfg['plots']['depth_clip']
        self.percent_min = ucfg.cfg['plots']['percent_min']
        self.percent_max = ucfg.cfg['plots']['percent_max']
        self.subs_fig = ucfg.cfg['plots']['subs_fig']
        self.density_flag = ucfg.cfg['plots']['density']
        self.swi_flag = ucfg.cfg['plots']['swi']
        self.current_image_flag = ucfg.cfg['plots']['current_image']
        self.image_change_flag = ucfg.cfg['plots']['image_change']
        self.flt_image_change_clims = ucfg.cfg['plots']['flt_image_change_clims']
        self.cold_content_flag = ucfg.cfg['plots']['cold_content']
        self.swe_volume_flag = ucfg.cfg['plots']['swe_volume']
        self.basin_total_flag = ucfg.cfg['plots']['basin_total']
        self.stn_validate_flag = ucfg.cfg['plots']['stn_validate']
        self.nash_sut_flag = ucfg.cfg['plots']['disp_nash_sut']
        self.stns_file = ucfg.cfg['plots']['stns_file']
        self.precip_validate_flag = ucfg.cfg['plots']['precip_validate']
        self.inputs_fig_flag = ucfg.cfg['plots']['inputs']
        self.plots_inputs_variables = ucfg.cfg['plots']['inputs_variables']
        self.compare_runs_flag = ucfg.cfg['plots']['compare_runs']
        self.compare_run_names = ucfg.cfg['plots']['compare_run_names']
        self.compare_run_labels = ucfg.cfg['plots']['compare_run_labels']
        self.compare_run_wys = ucfg.cfg['plots']['compare_run_wys']
        self.precip_depth_flag = ucfg.cfg['plots']['precip_depth']
        self.basin_detail_flag = ucfg.cfg['plots']['basin_detail']
        self.update_file = ucfg.cfg['plots']['update_file']
        self.print_args_dict = ucfg.cfg['plots']['print_args_dict']
        self.figsize = (ucfg.cfg['plots']['fig_length'],
                        ucfg.cfg['plots']['fig_height'])
        self.write_properties = ucfg.cfg['plots']['write_properties']
        self.point_values_flag = ucfg.cfg['plots']['point_values']

        if self.flt_image_change_clims[0] < 0:
            self.flt_image_change_clims[0] = 0
        if self.flt_image_change_clims[1] > 100:
            self.flt_image_change_clims[1] = 100

        if (self.write_properties is not None and
                type(self.write_properties) != list):
            self.write_properties = [self.write_properties]

        numbers = ucfg.cfg['plots']['update_numbers']

        if numbers is not None:
            if type(numbers) != list:
                numbers = [numbers]
            self.update_numbers = [x - 1 for x in numbers]
        else:
            self.update_numbers = None

        if (self.compare_runs_flag and ((self.compare_run_names is None) or
            (self.compare_run_labels is None) or self.compare_run_wys is None)):
            self.tmp_log.append(' Config option [plots] compare_runs set to '
                                'True, but one of compare_run_names, '
                                'compare_run_labels, or compare_run_wys is '
                                'empty, setting compare_runs to False')
            self.compare_runs_flag = False

        if (self.compare_runs_flag and
                (len(self.compare_run_names) != len(self.compare_run_labels))):
            self.tmp_log.append(' Config option [plots] compare_runs set to '
                                'True, must supply equal length '
                                'compare_run_names and compare_run_labels, '
                                'resetting compare_runs to False')
            self.compare_runs_flag = False

        if self.update_file is not None:
            self.flt_flag = True
        else:
            self.flt_flag = False

        if (self.precip_validate_flag and ((self.val_client is None) or
            (self.pre_val_stns is None) or (self.pre_val_lbls is None))):
            self.tmp_log.append(' Config option [plots] precip_validate is '
                                'being set to False')

            self.precip_validate_flag = False

        if (self.stn_validate_flag and (self.val_client is None) or
                (self.val_stns is None) or (self.val_lbls is None) or
                (self.wxdb_user is None) or (self.wxdb_password is None)):
            self.tmp_log.append(' Config option [plots] stn_validate is being '
                                'set to False')

            self.stn_validate_flag = False

        if len(self.point_values_settings) != 14:
            self.tmp_log.append(' Expected [validate] point_values_settings '
                                'to have 14 values, point_values set to False')
            self.point_values_flag = False

        for var in self.plots_inputs_variables:
            if var not in self.inputs_variables:
                self.plots_inputs_variables.remove(var)
                self.tmp_log.append(' Config option [plots] inputs_variables '
                                    'value {} not present in [diagnostics] '
                                    'inputs_variables, being '
                                    'removed'.format(var))

        ####################################################
        #          report                                  #
        ####################################################
        self.report_flag = ucfg.cfg['report']['report']
        self.print_latex = ucfg.cfg['report']['print_latex']
        self.report_name = ucfg.cfg['report']['file']
        self.rep_title = ucfg.cfg['report']['title']
        self.rep_path = ucfg.cfg['report']['save_path']
        self.env_path = ucfg.cfg['report']['env_path']
        self.templ_path = ucfg.cfg['report']['templ_path']
        self.tex_file = ucfg.cfg['report']['tex_file']
        self.summary_file = ucfg.cfg['report']['summary']
        self.figs_tpl_path = ucfg.cfg['report']['figs_tpl_path']
        self.flight_figs = ucfg.cfg['report']['flight_figs']
        self.tables = ucfg.cfg['report']['tables']
        self.report_diagnostics = ucfg.cfg['report']['diagnostics']
        self.report_diagnostics_day = ucfg.cfg['report']['diagnostics_day']
        self.rep_dplcs = ucfg.cfg['report']['decimals']

        if (self.report_diagnostics and
                (not self.inputs_fig_flag or not self.diagnostics_flag)):
            self.tmp_log.append(" [report] diagnostics: True, but must also "
                                "have [plots] inputs: True and [diagnostics] "
                                "diagnostics: True, setting to False")
            self.report_diagnostics = False

        if self.report_diagnostics and self.report_diagnostics_day[0] != 'any':

            if (calendar.day_name[datetime.now().weekday()] not in
                    self.report_diagnostics_day):
                self.report_diagnostics = False
                self.tmp_log.append(" Per [report] diagnostics_day: {}, "
                                    "setting diagnostics: "
                                    "False".format(self.report_diagnostics_day))

        self.rep_swi_flag = ucfg.cfg['report']['swi']
        if not self.swi_flag:
            self.rep_swi_flag = False

        self.rep_image_change_flag = ucfg.cfg['report']['image_change']
        if not self.image_change_flag:
            self.rep_image_change_flag = False

        self.rep_cold_content_flag = ucfg.cfg['report']['cold_content']
        if not self.cold_content_flag:
            self.rep_cold_content_flag = False

        self.rep_swe_volume_flag = ucfg.cfg['report']['swe_volume']
        if not self.swe_volume_flag:
            self.rep_swe_volume_flag = False

        self.rep_basin_total_flag = ucfg.cfg['report']['basin_total']
        if not self.basin_total_flag:
            self.rep_basin_total_flag = False

        self.rep_stn_validate_flag = ucfg.cfg['report']['stn_validate']
        if not self.stn_validate_flag:
            self.rep_stn_validate_flag = False

        self.rep_compare_runs_flag = ucfg.cfg['report']['compare_runs']
        if not self.compare_runs_flag:
            self.rep_compare_runs_flag = False

        self.rep_precip_depth_flag = ucfg.cfg['report']['precip_depth']
        if not self.precip_depth_flag:
            self.rep_precip_depth_flag = False

        # check paths to see if they need default snowav path
        if self.env_path is None:
            self.env_path = os.path.abspath(
                os.path.join(snowav.__path__[0],
                             "report/template/section_text"))
        if self.templ_path is None:
            self.templ_path = os.path.abspath(
                os.path.join(snowav.__path__[0],
                             "report/template"))
        if self.summary_file is None:
            self.summary_file = os.path.abspath(
                os.path.join(snowav.__path__[0],
                             "report/template/section_text/report_summary.txt"))
        if self.tex_file is None:
            self.tex_file = os.path.abspath(
                os.path.join(snowav.__path__[0],
                             "report/template/snowav_report.text"))
        if self.figs_tpl_path is None:
            self.figs_tpl_path = os.path.abspath(
                os.path.join(snowav.__path__[0],
                             "report/figs"))

        ####################################################
        #           forecast                               #
        ####################################################
        self.forecast_flag = ucfg.cfg['forecast']['report']

        if self.forecast_flag:
            self.for_start_date = ucfg.cfg['forecast']['start_date']
            self.for_end_date = ucfg.cfg['forecast']['end_date']
            self.for_run_name = ucfg.cfg['forecast']['run_name']

            if self.for_start_date >= self.for_end_date:
                self.tmp_log.append(' Error: config option [forecast] '
                                    'start_date > end_date')
                raise Exception('Config option [forecast] start_date >'
                                'end_date')

            self.for_run_dir = ([ucfg.cfg['forecast']['run_dir'] + s for s in
                                 os.listdir(ucfg.cfg['forecast']['run_dir'])
                                 if (os.path.isdir(ucfg.cfg['forecast']['run_dir'] + s))])

            self.for_run_dir.sort()

        ####################################################
        #           query                                  #
        ####################################################
        self.query_flag = ucfg.cfg['query']['query']
        self.q_basins = ucfg.cfg['query']['basins']
        self.q_value = ucfg.cfg['query']['value']
        self.q_run_name = ucfg.cfg['query']['run_name']
        self.q_print_all_runs = ucfg.cfg['query']['print_all_runs']
        self.q_start_date = ucfg.cfg['query']['start_date']
        self.q_end_date = ucfg.cfg['query']['end_date']
        self.q_total = ucfg.cfg['query']['total']
        self.q_output = ucfg.cfg['query']['output']
        self.q_csv_base_path = ucfg.cfg['query']['csv_base_path']
        self.q_database = ucfg.cfg['query']['database']

        ####################################################
        #           inflow                                 #
        ####################################################
        self.inflow_flag = ucfg.cfg['inflow']['inflow']
        self.inflow_data = ucfg.cfg['inflow']['inflow_data']
        self.summary_csv = ucfg.cfg['inflow']['summary_csv']
        self.inflow_headings = ucfg.cfg['inflow']['inflow_headings']
        self.basin_headings = ucfg.cfg['inflow']['basin_headings']
        self.sheet_name = ucfg.cfg['inflow']['sheet_name']
        self.skiprows = ucfg.cfg['inflow']['skiprows']
        self.overwrite = ucfg.cfg['inflow']['overwrite']
        self.file_base = ucfg.cfg['inflow']['file_base']
        self.date_idx = ucfg.cfg['inflow']['date_idx']
        self.convert = ucfg.cfg['inflow']['convert']

        self.ucfg = ucfg

    def assign_vars(self, kwargs):
        """ Assign attributes.

        Args
        ------
        kwargs {dict}: {item, value}
        """

        if not isinstance(kwargs, dict):
            raise TypeError("kwargs must be dict")

        for item, v in kwargs.items():
            setattr(self, item, v)

    def figure_names(self):
        """ Assign figure names. """

        et = self.end_date.date().strftime("%Y%-m%-d")
        basin = self.plotorder[0].lower().split(" ")[0]
        stn_validate_name = '{}_validation_{}.png'.format(basin, et)

        self.assign_vars({'stn_validate_fig_name': stn_validate_name})

    def parse(self, external_logger=None):
        """ Parse config options. """

        self.snowav_version = snowav.__version__
        self.cclimit = -5 * 1000 * 1000

        self.barcolors = ['xkcd:cobalt',
                          'xkcd:mustard green',
                          'xkcd:lichen',
                          'xkcd:pale green',
                          'xkcd:blue green',
                          'xkcd:bluish purple',
                          'xkcd:lightish purple',
                          'xkcd:deep magenta',
                          'xkcd:burgundy',
                          'red']

        out = masks(self.dempath, self.db_convert, plotorder=self.plotorder,
                    plotlabels=self.plotlabels)

        self.dem = out['dem']
        self.veg_type = out['veg_type']
        self.masks = out['masks']
        self.nrows = out['nrows']
        self.ncols = out['ncols']
        self.plotorder = out['plotorder']
        self.labels = out['labels']

        for log in out['logger']:
            self.tmp_log.append(log)

        # Establish database connection
        self.basins, cnx, out = connect(sqlite=self.sqlite, sql=self.mysql,
                                        plotorder=self.plotorder, user=self.db_user,
                                        password=self.db_password, host=self.db_host,
                                        port=self.db_port, convert=self.db_convert,
                                        add=self.add_basins)
        self.connector = cnx

        for log in out:
            self.tmp_log.append(log)

        if self.loglevel == 'DEBUG':
            for basin in self.basins:
                self.tmp_log.append(' {}: {}'.format(basin, self.basins[basin]))

        # Check snow.nc file location, get topo stats and water year
        sfile = os.path.join(self.run_dirs[0], 'snow.nc')

        if os.path.isfile(sfile):
            topo = get_topo_stats(sfile)
            self.snow_x = topo['x']
            self.snow_y = topo['y']
            self.pixel = int(topo['dv'])

            ncf = nc.Dataset(sfile)
            t = nc.num2date(ncf.variables['time'][0], ncf.variables['time'].units)
            ncf.close()
            self.wy = handle_year_stradling(t) + 1

        else:
            print('\nGiven config options, expecting to find:\n {}\nto load topo '
                  'stats but is not a valid file\nCheck config [run] options, see '
                  'CoreConfig.ini for details\n'.format(sfile))
            raise Exception('{} not a valid file'.format(sfile))

        # make the bins
        edges = np.arange(self.elev_bins[0],
                          self.elev_bins[1] + self.elev_bins[2],
                          self.elev_bins[2])

        # use for definition
        self.edges = np.arange(self.elev_bins[0] - self.elev_bins[2],
                               self.elev_bins[1],
                               self.elev_bins[2])

        v = self.properties
        if self.inputs_flag:
            v += self.inputs_variables

        # get variables that will be used in processing
        self.variables = AwsmInputsOutputs()
        self.variables.make_variables(v,
                                      self.edges,
                                      self.masks.keys())

        if self.units == 'TAF':
            self.conversion_factor = ((self.pixel ** 2) * 0.000000810713194 * 0.001)
            self.depth_factor = 0.03937
            self.dem = self.dem * 3.28
            self.depthlbl = 'in'
            self.vollbl = self.units
            self.elevlbl = 'ft'

            if max(self.edges) < 5000:
                self.tmp_log.append(" WARNING! Config options [snowav] units: TAF "
                                    "and elev_bins: {} may not match! Consider changing elev_bins "
                                    "values".format(self.elev_bins))

        if self.units == "SI":
            self.conversion_factor = ((self.pixel ** 2) * 0.000000810713194) * 1233.48 / 1e9
            self.depth_factor = 0.01
            self.depthlbl = 'cm'
            self.vollbl = 'M$M^3$'
            self.elevlbl = 'm'

            if max(self.edges) > 5000:
                self.tmp_log.append(" WARNING! Config options [snowav] units: SI "
                                    "and elev_bins: {} may not match! Consider changing elev_bins "
                                    "values".format(self.elev_bins))

        self.ixd = np.digitize(self.dem, edges)
        self.xlims = (0, len(edges))

        if self.loglevel == 'DEBUG' and self.log_to_file is not True:
            print('Reading files in {}...'.format(self.run_dirs[0].split('runs')[0]))

        results = outputs(self.run_dirs, self.wy, self.properties,
                          self.start_date, self.end_date, None, self.loglevel)

        out = results['outputs']
        all_dirs = results['dirs']
        dirs = results['run_dirs']
        rdict = results['rdict']
        log = results['log']

        # If there was an error parsing files catch and log it
        if out == [] and all_dirs == [] and 'not a valid file' in log[-1]:

            self.tmp_log.append(log[-1])
            if self.start_date is not None and self.end_date is not None:
                ext_shr = (self.directory +
                           '_' +
                           self.start_date.date().strftime("%Y%m%d") +
                           '_' +
                           self.end_date.date().strftime("%Y%m%d"))
                self.figs_path = os.path.join(self.save_path, '{}/'.format(ext_shr))

                if external_logger == None:
                    createLog(self)
                else:
                    self._logger = external_logger

            raise Exception(log[-1])

        for l in log:
            self.tmp_log.append(l)

        if out['dates'] == []:
            raise Exception('Supplied [run] directory, start_date, and end_date '
                            'give no valid snow files')

        self.outputs = out
        self.run_dirs = dirs
        self.all_dirs = all_dirs
        self.rundirs_dict = rdict
        self.all_dirs_flt = deepcopy(all_dirs)

        if self.start_date is not None and self.end_date is None:
            self.end_date = self.outputs['dates'][-1]
            self.tmp_log.append(' Config options [run] end_date '
                                'not specified, assigning '
                                '{} and {}'.format(self.start_date, self.end_date))

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
        if self.start_date == datetime(self.wy - 1, 10, 1, 23, 0, 0):
            self.report_start = self.start_date

        else:
            self.report_start = self.start_date + timedelta(hours=1)

        # Copy the config file where figs will be saved
        # use directory if only plotting figures from database and don't
        # have start_date, end_date
        extf = os.path.splitext(os.path.split(self.config_file)[1])
        ext_shr = (self.directory +
                   '_' +
                   self.start_date.date().strftime("%Y%m%d") +
                   '_' +
                   self.end_date.date().strftime("%Y%m%d"))
        self.figs_path = os.path.join(self.save_path, '{}/'.format(ext_shr))

        # get forecast outputs
        if self.forecast_flag:
            results = outputs(self.for_run_dirs, self.wy, self.properties,
                              None, None, None, self.loglevel)

            self.for_outputs = results['outputs']
            self.for_run_dirs = results['run_dirs']
            self.for_rundirs_dict = results['rdict']
            self.for_ixs = 0
            self.for_ixe = len(self.for_outputs['swe_z']) - 1

        # Get outputs for flights
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
                pre_wydate = calculate_date_from_wyhr(int(time - 24), self.wy)
                flight_dates = np.append(flight_dates, wydate)
                pre_flight_dates = np.append(pre_flight_dates, pre_wydate)

            if self.loglevel == 'DEBUG' and self.log_to_file is not True:
                print('Reading files in {} for flight updates'
                      '...'.format(self.run_dirs[0].split('runs')[0]))

            results = outputs(self.all_dirs_flt, self.wy, self.properties,
                              None, None, flight_dates, self.loglevel)

            self.flight_outputs = results['outputs']
            self.run_dirs_flt = results['run_dirs']
            self.flt_rundirs_dict = results['rdict']
            self.flight_diff_dates = results['outputs']['dates']
            self.pre_flight_outputs = results['outputs']

            results = outputs(self.all_dirs_flt, self.wy, self.properties,
                              None, None, pre_flight_dates, self.loglevel)

            self.pre_flight_outputs = results['outputs']

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

    field_styles = {'hostname': {'color': 'magenta'},
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
