
import os
from snowav.framework.query import query
from snowav.framework.read_config import read_config
from snowav.framework.process import process
from snowav.framework.parse import parse
from snowav.framework.figures import figures
from snowav.report.report import report
from snowav.database.database import run_metadata
from snowav.inflow.inflow import excel_to_csv
from datetime import datetime

class snowav(object):

    def __init__(self, config_file = None, external_logger = None, awsm = None,
                 end_date = None):
        '''
        Read config file, parse config options, process results, put results
        on database, make figures, and make pdf report.

        Args
        -----
        config_file : str
            snowav config file
        external_logger : object
            awsm logger
        awsm : class
            awsm class if being run in awsm

        '''

        # Get config options
        if awsm is None and os.path.isfile(config_file):
            self.config_file = config_file
            read_config(self, end_date = end_date)

        elif awsm is not None:
            self.config_file = awsm.configFile
            read_config(self, awsm = awsm, end_date = end_date)

        else:
            raise Exception('No config instance passed, or config file does '
                            'not exist!')

        # query existing database without processing
        if self.query_flag:
            query(self)

        # parse config options
        parse(self)

        # put run metadata on database
        run_metadata(self, self.run_name)

        # process results
        self.pargs['run_id'] = self.run_id
        self.pargs['vid'] = self.vid
        self.precip_flag, out, pre, rain, density = process(self.pargs)

        # gather process() outputs
        # for log in out:
        #     self._logger.info(log)

        self.density = density
        self.rain_total = rain
        self.precip_total = pre

        if not self.precip_flag:
            self.precip_depth_flag = False
            self.precip_validate_flag = False

        if self.inflow_flag and self.inflow_data is not None:
            args = {}
            args['path'] = self.inflow_data
            args['csv_file'] = self.summary_csv
            args['basin_headings'] = self.basin_headings
            args['inflow_headings'] = self.inflow_headings
            args['file_base'] = self.file_base
            args['sheet_name'] = self.sheet_name
            args['skiprows'] = self.skiprows
            args['date_idx'] = self.date_idx
            args['wy'] = self.wy
            args['overwrite'] = self.overwrite
            args['convert'] = self.convert

            excel_to_csv(args, self._logger)

        figures(self)

        # Do additional processing and figures if forecast is supplied. Some
        # field will be overwritten during forecast processing
        if self.forecast_flag:
            self._logger.info(' Starting forecast processing...')

            run_metadata(self, self.for_run_name)

            self.pargs['run_id'] = self.run_id
            self.pargs['vid'] = self.vid
            flags, out, pre, rain, density = process(self.pargs)

            for log in out:
                self._logger.info(log)

            self.density = density
            self.rain_total = rain
            self.precip_total = pre

            # figures for forecast run
            figures(self)

        if self.report_flag:
            report(self)

        elapsed = str(datetime.now() - self.proc_time_start)

        self._logger.info(' Completed snowav processing, elapsed time: {}'.format(elapsed))
