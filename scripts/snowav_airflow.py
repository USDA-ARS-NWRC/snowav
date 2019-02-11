#!/usr/bin/env python3
# -*- coding: utf-8 -*-
from inicheck.config import MasterConfig, UserConfig
from inicheck.tools import get_user_config, check_config
from inicheck.output import print_config_report, generate_config
from inicheck.tools import cast_all_variables
from smrf.utils import utils, io
# from awsm.framework.framework import run_awsm
from scripts import snow
from datetime import datetime
from datetime import timedelta
import sys
import os
import argparse
import numpy as np
import copy


def mod_config(config_file):
    """
    Run daily snowav. Currently, this will only edit the end_date in the
    config, so daily snowav runs will happen with the same start_date.
    """
    # define some formats
    fmt_cfg = '%Y-%m-%d 23:00'

    # get config instance
    config = get_user_config(config_file, modules = ['snowav'])
    config.apply_recipes()
    config = cast_all_variables(config, config.mcfg)

    end_date = (datetime.now() - timedelta(hours=24)).date()

    new_config = copy.deepcopy(config)

    new_config.raw_cfg['outputs']['end_date'] = end_date.strftime(fmt_cfg)

    return new_config


def run():
    '''

    '''

    # Parse arguments
    p = argparse.ArgumentParser(description='Run SNOWAV using Airflow scheduler.')

    p.add_argument('-c', '--cfg', required=True,
                   help='Config file that will be modified for the current run date')

    args = p.parse_args()

    fp_cfg = args.cfg

    # set dates and paths
    new_config = mod_config(fp_cfg)

    # apply recipes with new setttings
    # new_config.apply_recipes()
    # new_config = cast_all_variables(new_config, new_config.mcfg)

    # run_awsm(new_config)
    snow.run(config_instance = new_config)


if __name__ == '__main__':
    run()
