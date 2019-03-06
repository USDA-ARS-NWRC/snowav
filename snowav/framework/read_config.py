
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
import logging
import coloredlogs
import netCDF4 as nc
from dateutil.relativedelta import relativedelta
# from sqlalchemy import create_engine
# from sqlalchemy.ext.declarative import declarative_base
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

    '''
    # print('Reading SNOWAV config file...')
    if awsm is None:
        ucfg = get_user_config(self.config_file, modules = 'snowav')

    else:
        ucfg = awsm

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
    # self.vars = ucfg.cfg['outputs']['vars']

    if awsm is None:
        self.start_date = ucfg.cfg['outputs']['start_date']
        self.end_date = ucfg.cfg['outputs']['end_date']

    else:
        fmt_cfg = '%Y-%m-%d 23:00'
        end_date = (datetime.datetime.now() - timedelta(hours=24)).date()
        end_date = pd.to_datetime(end_date.strftime(fmt_cfg))

        self.start_date = ucfg.cfg['outputs']['start_date']
        self.end_date = end_date

    if (self.start_date is not None and self.end_date is not None):
        self.start_date = self.start_date.to_pydatetime()
        self.end_date = self.end_date.to_pydatetime()

        if self.start_date >= self.end_date:
            self.tmp_log.append('Error: [outputs]->start_date > [outputs]->end_date')
            exit()

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
    #           forecast                               #
    ####################################################
    self.forecast_flag = ucfg.cfg['forecast']['report']

    if self.forecast_flag is True:
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

    self.pre_val_stns = ucfg.cfg['validate']['pre_stations']
    self.pre_val_lbls = ucfg.cfg['validate']['pre_labels']
    self.val_client = ucfg.cfg['validate']['client']

    # This is being used to combine 2017 HRRR data
    self.offset = int(ucfg.cfg['validate']['offset'])

    ####################################################
    #           basin total                            #
    ####################################################
    self.flight_dates = ucfg.cfg['basin total']['flights']
    if (self.flight_dates is not None) and (type(self.flight_dates) != list):
        self.flight_dates = [self.flight_dates]

    ####################################################
    #           masks                                  #
    ####################################################
    self.dempath = ucfg.cfg['masks']['dempath']
    self.plotorder = ucfg.cfg['masks']['mask_labels']

    if ucfg.cfg['masks']['basin_masks'] is not None:
        self.total = ucfg.cfg['masks']['basin_masks'][0]

    ####################################################
    #          plots                                   #
    ####################################################
    self.figsize = (ucfg.cfg['plots']['fig_length'],
                    ucfg.cfg['plots']['fig_height'])
    self.dpi = ucfg.cfg['plots']['dpi']
    self.annot_x = ucfg.cfg['plots']['annot_x']
    self.annot_y = ucfg.cfg['plots']['annot_y']
    self.subs_fig = ucfg.cfg['plots']['subs_fig']
    self.flow_file = ucfg.cfg['plots']['flow_file']
    self.density_flag = ucfg.cfg['plots']['density']
    self.subbasins_flag = ucfg.cfg['plots']['subbasins']
    self.inflow_flag = ucfg.cfg['plots']['inflow']
    self.accumulated_flag = ucfg.cfg['plots']['accumulated']
    self.current_image_flag = ucfg.cfg['plots']['current_image']
    self.image_change_flag = ucfg.cfg['plots']['image_change']
    self.cold_content_flag = ucfg.cfg['plots']['cold_content']
    self.swe_volume_flag = ucfg.cfg['plots']['swe_volume']
    self.swe_change_flag = ucfg.cfg['plots']['swe_change']
    self.basin_total_flag = ucfg.cfg['plots']['basin_total']
    self.pixel_swe_flag = ucfg.cfg['plots']['pixel_swe']
    self.stn_validate_flag = ucfg.cfg['plots']['stn_validate']
    self.precip_validate_flag = ucfg.cfg['plots']['precip_validate']
    self.compare_runs_flag = ucfg.cfg['plots']['compare_runs']
    self.precip_depth_flag = ucfg.cfg['plots']['precip_depth']
    self.basin_detail_flag = ucfg.cfg['plots']['basin_detail']

    ####################################################
    #          report                                  #
    ####################################################
    self.report_flag = ucfg.cfg['report']['report']
    self.exclude_figs = ucfg.cfg['report']['exclude_figs']

    if type(self.exclude_figs) != list and self.exclude_figs != None:
        self.exclude_figs = [self.exclude_figs]

    elif self.exclude_figs is None:
        self.exclude_figs = []

    else:
        if self.subbasins_flag is False:
            self.exclude_figs.append('SUBBASINS')

        if self.inflow_flag is False:
            self.exclude_figs.append('INFLOW')

        if self.stn_validate_flag is False:
            self.exclude_figs.append('VALID')

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
    self.report_only = ucfg.cfg['results']['report_only']
    self.write_stn_csv_flag = ucfg.cfg['results']['write_stn_csv']

    if self.write_stn_csv_flag is True:
        self.stns_csv = pd.read_csv(ucfg.cfg['results']['stn_csv_file'])
        self.stns_csv.set_index('name')

    self.run_name = ucfg.cfg['results']['run_name']
    self.write_db = ucfg.cfg['results']['write_db']
    self.figures_only = ucfg.cfg['results']['figures_only']
    self.plot_runs = ucfg.cfg['results']['plot_runs']
    self.plot_labels = ucfg.cfg['results']['plot_labels']
    self.plot_variables = ucfg.cfg['results']['plot_variables']

    if (self.compare_runs_flag is True) and (self.plot_runs is None):
        self.tmp_log.append('No runs listed in [results] -> plot_runs, so [] ' +
        '-> being set to False')
        self.compare_runs_flag = False

    if self.figures_only is True:
        self.plot_flag = True

        self.tmp_log.append('Config file option [results]->figures_only={}, '
              'creating figures from database...'.format(self.figures_only))
    else:
        self.plot_flag = False
