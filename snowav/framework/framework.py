from datetime import datetime
import pandas as pd
from sys import exit

from snowav.framework.query import query
from snowav.config.config import UserConfig
from snowav.framework.process import Process
from snowav.framework.figures import figures
from snowav.report.report import report
from snowav.database.database import Database, run_metadata
from snowav.inflow.inflow import excel_to_csv


class Snowav(object):

    def __init__(self, config_file=None, awsm=None, end_date=None):
        """ Read config file, parse config options, process results, put
        results on database, make figures, and make pdf report.

        Args
        -----
        config_file {str}: config file
        awsm {class}: awsm class
        end_date {str}: overwrite of config end_date
        """

        if end_date is not None:
            end_date = pd.to_datetime(end_date)

        # get and parse config options
        cfg = UserConfig(config_file, awsm=awsm, end_date=end_date)
        cfg.parse()
        cfg.figure_names()

        if cfg.report_only:
            report(cfg)
            elapsed = str(datetime.now() - cfg.proc_time_start)
            cfg._logger.info(' Completed snowav processing, '
                             'elapsed time: {}'.format(elapsed))
            exit()

        if cfg.query_flag:
            query(cfg)

        run_metadata(cfg)
        if cfg.db_type == 'sql':
            dbp = cfg.mysql
        else:
            dbp = cfg.sqlite

        db = Database(user=cfg.db_user, password=cfg.db_password,
                      host=cfg.db_host, port=cfg.db_port, db_type=cfg.db_type,
                      database=dbp)
        db.make_connection(logger=cfg._logger)
        db.check_tables(logger=cfg._logger)

        process = Process(cfg, db)

        if cfg.inflow_flag and cfg.inflow_data is not None:
            args = {'path': cfg.inflow_data,
                    'csv_file': cfg.summary_csv,
                    'basin_headings': cfg.basin_headings,
                    'inflow_headings': cfg.inflow_headings,
                    'file_base': cfg.file_base,
                    'sheet_name': cfg.sheet_name,
                    'skiprows': cfg.skiprows,
                    'date_idx': cfg.date_idx,
                    'wy': cfg.wy,
                    'overwrite': cfg.overwrite,
                    'convert': cfg.convert}

            excel_to_csv(args, cfg._logger)

        figures(cfg, process, db)

        if cfg.report_flag:
            report(cfg)

        elapsed = str(datetime.now() - cfg.proc_time_start)

        cfg._logger.info(' Completed snowav processing, '
                         'elapsed time: {}'.format(elapsed))
