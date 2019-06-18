
import logging
import numpy as np
from spatialnc import ipw
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
import coloredlogs
import netCDF4 as nc
from dateutil.relativedelta import relativedelta
from sqlalchemy.orm import sessionmaker
from collections import OrderedDict
from snowav import database
from datetime import timedelta
from shutil import copyfile
from sys import exit


def read_config(self, external_logger=None, awsm=None):
    '''
    Read snowav config file and assign fields.

    Args
        external_logger: awsm logger
        awsm: awsm class, if this is passed, run_dir will be assigned from the
            directory being created in awsm

            awsm_mcfg = MasterConfig(modules = 'awsm')
            smrf_mcfg = MasterConfig(modules = 'smrf')

            # Read in the original users config
            self.ucfg = get_user_config(configFile, mcfg=combined_mcfg)
            self.configFile = configFile

    '''

    snowav_mcfg = MasterConfig(modules = 'snowav')
    ucfg = get_user_config(self.config_file, mcfg=snowav_mcfg)
    self.snowav_path = get_snowav_path()

    warnings, errors = check_config(ucfg)
    if errors != [] or warnings != []:
        print_config_report(warnings, errors)

    # create blank log and error log because logger is not initialized yet
    self.tmp_log = []
    self.tmp_err = []
    self.tmp_warn = []

    ####################################################
    #             snowav                               #
    ####################################################
    self.loglevel = ucfg.cfg['snowav']['log_level'].upper()
    self.log_to_file = ucfg.cfg['snowav']['log_to_file']
    self.save_path = ucfg.cfg['snowav']['save_path']
    self.units = ucfg.cfg['snowav']['units']
    self.filetype = ucfg.cfg['snowav']['filetype']
    self.elev_bins = ucfg.cfg['snowav']['elev_bins']
    self.directory = ucfg.cfg['snowav']['directory']
    self.dempath = ucfg.cfg['snowav']['dempath']
    self.run_name = ucfg.cfg['snowav']['run_name']
    self.plotorder = ucfg.cfg['snowav']['masks']

    if type(self.plotorder) != list:
        self.plotorder = [self.plotorder]

    self.plotlabels = ucfg.cfg['snowav']['plotlabels']

    if type(self.plotlabels) != list:
        self.plotlabels = [self.plotlabels]

    ####################################################
    #           run                                    #
    ####################################################

    self.dplcs = ucfg.cfg['run']['decimals']
    self.start_date = ucfg.cfg['run']['start_date']
    self.end_date = ucfg.cfg['run']['end_date']

    if self.start_date is not None and self.end_date is not None:
        self.start_date = self.start_date.to_pydatetime()
        self.end_date = self.end_date.to_pydatetime()

        if self.start_date >= self.end_date:
            self.tmp_log.append('Error: [run] start_date >= [run] end_date')
            raise Exception('Error: [run] start_date >= [run] end_date')

    else:
        self.tmp_log.append('[run] start_date and/or end_date was not '
                            'defined in config file, will be assigned as '
                            'first and last available dates in run_dirs')

    self.all_subdirs = ucfg.cfg['run']['all_subdirs']

    if self.all_subdirs is True:
        self.run_dirs = ([ucfg.cfg['run']['run_dirs'] + s for s in
                        os.listdir(ucfg.cfg['run']['run_dirs'])
                        if (os.path.isdir(ucfg.cfg['run']['run_dirs'] + s)) ])
    else:
        self.run_dirs = ucfg.cfg['run']['run_dirs']
        if type(self.run_dirs) != list:
            self.run_dirs = [self.run_dirs]

    self.run_dirs.sort()

    # pull from num2date nc snow
    self.wy = 2019

    # if awsm is not None:
        # fmt_cfg = '%Y-%m-%d 23:00'
        # end_date = (datetime.datetime.now() - timedelta(hours=24)).date()
        # end_date = pd.to_datetime(end_date.strftime(fmt_cfg))
        #
        # self.start_date = ucfg.cfg['outputs']['start_date']
        # self.end_date = end_date

    # self.summary = ucfg.cfg['outputs']['summary']
    # if type(self.summary) != list:
    #     self.summary = [self.summary]

    ####################################################
    #           forecast                               #
    ####################################################
    self.forecast_flag = ucfg.cfg['forecast']['report']

    if self.forecast_flag:
        self.for_start_date = ucfg.cfg['forecast']['start_date'].to_pydatetime()
        self.for_end_date = ucfg.cfg['forecast']['end_date'].to_pydatetime()
        self.for_run_name = ucfg.cfg['forecast']['run_name']

        if self.for_start_date >= self.for_end_date:
            self.tmp_log.append('Error: [outputs]->start_date > [outputs]->end_date')
            exit()

        self.for_run_dir = ([ucfg.cfg['forecast']['run_dir'] + s for s in
                        os.listdir(ucfg.cfg['forecast']['run_dir'])
                        if (os.path.isdir(ucfg.cfg['forecast']['run_dir'] + s)) ])

        self.for_run_dir.sort()

    ####################################################
    #         database
    ####################################################

    self.mysql = ucfg.cfg['database']['mysql']
    self.db_user = ucfg.cfg['database']['user']
    self.db_password = ucfg.cfg['database']['password']
    self.db_host = ucfg.cfg['database']['host']
    self.db_port = ucfg.cfg['database']['port']
    if ((self.mysql is not None) and
        ((self.db_user is None) or
         (self.db_password is None) or
         (self.db_host is None) or
         (self.db_port is None)) ):
        raise Exception('If using config option [database] mysql:, must also ' +
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
        raise Exception('Config [databasee] section contains both mysql and ' +
                        'sqlite entries, please pick one...')

    ####################################################
    #           validate                               #
    ####################################################
    self.val_stns = ucfg.cfg['validate']['stations']
    self.val_lbls = ucfg.cfg['validate']['labels']
    self.val_client = ucfg.cfg['validate']['client']
    self.pre_val_stns = ucfg.cfg['validate']['pre_stations']
    self.pre_val_lbls = ucfg.cfg['validate']['pre_labels']

    ####################################################
    #          plots                                   #
    ####################################################
    self.print_args_dict = ucfg.cfg['plots']['print_args_dict']
    self.figsize = (ucfg.cfg['plots']['fig_length'],
                    ucfg.cfg['plots']['fig_height'])
    self.dpi = ucfg.cfg['plots']['dpi']
    self.depth_clip = ucfg.cfg['plots']['depth_clip']
    self.percent_min = ucfg.cfg['plots']['percent_min']
    self.percent_max = ucfg.cfg['plots']['percent_max']
    self.annot_x = ucfg.cfg['plots']['annot_x']
    self.annot_y = ucfg.cfg['plots']['annot_y']
    self.subs_fig = ucfg.cfg['plots']['subs_fig']
    self.flow_file = ucfg.cfg['plots']['flow_file']
    self.density_flag = ucfg.cfg['plots']['density']
    self.inflow_flag = ucfg.cfg['plots']['inflow']
    self.swi_flag = ucfg.cfg['plots']['swi']
    self.current_image_flag = ucfg.cfg['plots']['current_image']
    self.image_change_flag = ucfg.cfg['plots']['image_change']
    self.cold_content_flag = ucfg.cfg['plots']['cold_content']
    self.swe_volume_flag = ucfg.cfg['plots']['swe_volume']
    self.swe_change_flag = ucfg.cfg['plots']['swe_change']
    self.basin_total_flag = ucfg.cfg['plots']['basin_total']
    self.pixel_swe_flag = ucfg.cfg['plots']['pixel_swe']
    self.stn_validate_flag = ucfg.cfg['plots']['stn_validate']
    self.nash_sut_flag = ucfg.cfg['plots']['disp_nash_sut']
    self.stns_file = ucfg.cfg['plots']['stns_file']
    self.precip_validate_flag = ucfg.cfg['plots']['precip_validate']
    self.compare_runs_flag = ucfg.cfg['plots']['compare_runs']
    self.precip_depth_flag = ucfg.cfg['plots']['precip_depth']
    self.basin_detail_flag = ucfg.cfg['plots']['basin_detail']

    self.update_file = ucfg.cfg['plots']['update_file']
    self.update_numbers = ucfg.cfg['plots']['update_numbers']

    self.plot_runs = ucfg.cfg['plots']['plot_runs']
    self.plot_labels = ucfg.cfg['plots']['plot_labels']
    self.plot_variables = ucfg.cfg['plots']['plot_variables']

    if (self.compare_runs_flag is True) and (self.plot_runs is None):
        self.tmp_log.append('No runs listed in [results] -> plot_runs, so [] ' +
        '-> being set to False')
        self.compare_runs_flag = False

    if self.update_file is not None:
        self.flt_flag = True
    else:
        self.flt_flag = False

    if ((self.precip_validate_flag) is True and
       (self.val_client is None) or
       (self.pre_val_stns is None) or
       (self.pre_val_lbls is None)):
        self.tmp_log.append('[plots] -> precip_validate is being set to '
            'to False, either [validate] -> pre_stations, pre_labels or '
            'client is empty')

        self.precip_validate_flag = False

    if ((self.stn_validate_flag) is True and
       (self.val_client is None) or
       (self.val_stns is None) or
       (self.val_lbls is None)):
        self.tmp_log.append('[plots] -> stn_validate is being set to '
            'to False, either [validate] -> stations, labels or '
            'client is empty')

        self.stn_validate_flag = False

    # if self.stn_validate_flag:
    # check that clients, etc., are there before getting all the way to figure

    ####################################################
    #          report                                  #
    ####################################################
    self.report_flag = ucfg.cfg['report']['report']
    self.report_name = ucfg.cfg['report']['report_name']
    self.rep_title = ucfg.cfg['report']['report_title']
    self.rep_path = ucfg.cfg['report']['rep_path']
    self.env_path = ucfg.cfg['report']['env_path']
    self.templ_path = ucfg.cfg['report']['templ_path']
    self.tex_file = ucfg.cfg['report']['tex_file']
    self.summary_file = ucfg.cfg['report']['summary_file']
    self.figs_tpl_path = ucfg.cfg['report']['figs_tpl_path']
    self.swi_flag = ucfg.cfg['report']['swi']
    self.current_image_flag = ucfg.cfg['report']['current_image']
    self.image_change_flag = ucfg.cfg['report']['image_change']
    self.cold_content_flag = ucfg.cfg['report']['cold_content']
    self.swe_volume_flag = ucfg.cfg['report']['swe_volume']
    self.basin_total_flag = ucfg.cfg['report']['basin_total']
    self.pixel_swe_flag = ucfg.cfg['report']['pixel_swe']
    self.stn_validate_flag = ucfg.cfg['report']['stn_validate']
    self.precip_validate_flag = ucfg.cfg['report']['precip_validate']
    self.compare_runs_flag = ucfg.cfg['report']['compare_runs']
    self.precip_depth_flag = ucfg.cfg['report']['precip_depth']

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

    self.ucfg = ucfg
