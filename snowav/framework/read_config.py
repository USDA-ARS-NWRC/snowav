
from snowav.utils.utilities import get_snowav_path
from inicheck.tools import get_user_config, check_config, cast_all_variables
from inicheck.output import generate_config, print_config_report
from inicheck.config import MasterConfig
import logging
import coloredlogs
import os
import pandas as pd
from datetime import datetime
import copy
import calendar

def read_config(self, external_logger = None, awsm = None, end_date = None):
    '''
    Read snowav config file.

    '''

    print('Reading {} and loading files...'.format(self.config_file))

    snowav_mcfg = MasterConfig(modules = 'snowav')
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
        self.start_date = self.start_date.to_pydatetime()
        self.end_date = self.end_date.to_pydatetime()

        if self.start_date >= self.end_date:
            self.tmp_log.append(' Error: [run] start_date >= [run] end_date')
            raise Exception('Error: [run] start_date >= [run] end_date')

    else:
        self.tmp_log.append(' [run] start_date and/or end_date was not '
                            'defined in config file, will be assigned '
                            'by available dates in directory')

    self.all_subdirs = ucfg.cfg['run']['all_subdirs']

    if (ucfg.cfg['run']['directory'] is None) and (awsm is not None):
        if self.all_subdirs is True:
            self.run_dirs = ([awsm.pathr + s for s in
                            os.listdir(awsm.pathr)
                            if (os.path.isdir(awsm.pathr + s)) ])
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

    base_bands = ['swi_z','evap_z','swe_z','depth','density',
        'coldcont', 'precip_z']

    f = False
    for band in base_bands:
        if band not in self.properties:
            self.tmp_log.append(' WARNING! Config option [database] properties '
                'does not contain {}'.format(band))
            f = True

    if f:
        self.tmp_log.append(' WARNING! Suggest config option [database] '
            'properties contain at least {} for most functionality to '
            'run'.format(base_bands))

    if ((self.mysql is not None) and
        ((self.db_user is None) or
         (self.db_password is None) or
         (self.db_host is None) or
         (self.db_port is None)) ):
        raise Exception('If using config option [database] mysql, must also '
                        'supply user, password, host, and port')

    sqlite = ucfg.cfg['database']['sqlite']

    if (sqlite is not None) and os.path.isdir(os.path.dirname(sqlite)):
        self.sqlite = 'sqlite:///' + sqlite
    elif sqlite is None:
        self.sqlite = None

    if (sqlite is not None) and not os.path.isdir(os.path.dirname(sqlite)):
        raise Exception('Config option [database] sqlite: {} '.format(sqlite) +
                        'contains an invalid base path for the sqlite database')

    if self.mysql is not None and sqlite is not None:
        raise Exception('Config option [database] section contains both mysql '
                        'and sqlite entries, pick one...')

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

    if len(self.point_values_settings) != 14:
        self.tmp_log.append(' Expected [validate] point_values_settings to '
                            'have 14 values, point_values being set to False')
        self.point_values = False

    for n in range(0,10):
        self.point_values_settings[n] = int(self.point_values_settings[n])

    if type(self.point_values_properties) != list:
        self.point_values_properties = [self.point_values_properties]

    if type(self.point_values_heading) != list:
        self.point_values_heading = [self.point_values_heading]

    if type(self.point_values_csv) != list:
        self.point_values_csv = [self.point_values_csv]

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

    if not (len(self.point_values_csv) == len(self.point_values_properties) ==
            len(self.point_values_heading)):
        self.tmp_log.append(' Must supply the same number of [validate] '
                            'point_values_csv, point_values_properties, '
                            'and point_values_heading values, point_values '
                            'being set to False')


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
                    self.tmp_log.append(' Config option [diagnostics] basin: {} '
                                        'does not match what was supplied in '
                                        '[snowav] masks: {}, diagnostics set '
                                        'to False'.format(basin, self.plotorder))
                    self.diagnostics_flag = False

    self.inputs_flag = ucfg.cfg['diagnostics']['inputs_table']
    self.inputs_variables = ucfg.cfg['diagnostics']['inputs_variables']
    self.inputs_percentiles = ucfg.cfg['diagnostics']['inputs_percentiles']
    self.inputs_methods = ucfg.cfg['diagnostics']['inputs_methods']
    self.inputs_basins = ucfg.cfg['diagnostics']['inputs_basins']

    # if self.inputs_flag:
    if self.inputs_basins is not None and self.plotorder is not None:
        for basin in self.inputs_basins:
            if basin not in self.plotorder:
                self.tmp_log.append(' Config option [diagnostics] '
                                    'inputs_basins: {} does not match what '
                                    'was supplied in [snowav] masks: {}, '
                                    'inputs set to '
                                    'False'.format(basin, plotorder))
                self.inputs_flag = False

    if self.inputs_flag:
        s = [x + ', ' for x in self.inputs_variables]
        self.tmp_log.append(' Using variables {} for inputs summary'.format(''.join(s)))

        s = [x + ', ' for x in self.inputs_methods]
        self.tmp_log.append(' Using methods {} for inputs summary'.format(''.join(s)))

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
    if self.write_properties is not None and type(self.write_properties) != list:
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
        self.tmp_log.append(' Config option [plots] compare_runs set to True, '
                            'but one of compare_run_names, compare_run_labels, '
                            'or compare_run_wys is empty, setting '
                            'compare_runs to False')
        self.compare_runs_flag = False

    if (self.compare_runs_flag and
        (len(self.compare_run_names) != len(self.compare_run_labels))):
        self.tmp_log.append(' Config option [plots] compare_runs set to True, '
                            'must supply equal length compare_run_names and  '
                            'compare_run_labels, resetting compare_runs to False')
        self.compare_runs_flag = False

    if self.update_file is not None:
        self.flt_flag = True

    else:
        self.flt_flag = False

    if (self.precip_validate_flag and ((self.val_client is None) or
       (self.pre_val_stns is None) or (self.pre_val_lbls is None))):
        self.tmp_log.append(' Config option [plots] precip_validate is being '
                            'set to False')

        self.precip_validate_flag = False

    if (self.stn_validate_flag and (self.val_client is None) or
       (self.val_stns is None) or (self.val_lbls is None) or
       (self.wxdb_user is None) or (self.wxdb_password is None) ):
        self.tmp_log.append(' Config option [plots] stn_validate is being '
                            'set to False')

        self.stn_validate_flag = False

    for var in self.plots_inputs_variables:
        if var not in self.inputs_variables:
            self.plots_inputs_variables.remove(var)
            self.tmp_log.append(' Config option [plots] inputs_variables '
                                'value {} not present in [diagnostics] '
                                'inputs_variables, being removed'.format(var))


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

    if self.report_diagnostics and (not self.inputs_fig_flag or not self.diagnostics_flag):
        self.tmp_log.append(" [report] diagnostics: True, but must also have "
            "[plots] inputs: True and [diagnostics] diagnostics: True, "
            "setting to False")
        self.report_diagnostics = False

    if self.report_diagnostics and self.report_diagnostics_day[0] != 'any':

        if calendar.day_name[datetime.now().weekday()] not in self.report_diagnostics_day:
            self.report_diagnostics = False
            self.tmp_log.append(" Per [report] diagnostics_day: {}, "
                "setting diagnostics: False".format(self.report_diagnostics_day))

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
        self.env_path = os.path.abspath('./snowav/report/template/section_text')
    if self.templ_path is None:
        self.templ_path = os.path.abspath('./snowav/report/template')
    if self.summary_file is None:
        self.summary_file = os.path.abspath('./snowav/report/template/section_text/report_summary.txt')
    if self.tex_file is None:
        self.tex_file = os.path.abspath('./snowav/report/template/snowav_report.tex')
    if self.figs_tpl_path is None:
        self.figs_tpl_path = os.path.abspath('./snowav/report/figs')

    ####################################################
    #           forecast                               #
    ####################################################

    self.forecast_flag = ucfg.cfg['forecast']['report']

    if self.forecast_flag:
        self.for_start_date = ucfg.cfg['forecast']['start_date'].to_pydatetime()
        self.for_end_date = ucfg.cfg['forecast']['end_date'].to_pydatetime()
        self.for_run_name = ucfg.cfg['forecast']['run_name']

        if self.for_start_date >= self.for_end_date:
            self.tmp_log.append(' Error: config option [forecast] start_date > '
                                'end_date')
            raise Exception('Config option [forecast] start_date is greater '
                            'than end_date')

        self.for_run_dir = ([ucfg.cfg['forecast']['run_dir'] + s for s in
                        os.listdir(ucfg.cfg['forecast']['run_dir'])
                        if (os.path.isdir(ucfg.cfg['forecast']['run_dir'] + s)) ])

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
