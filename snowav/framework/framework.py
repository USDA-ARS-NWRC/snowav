
import os
from sys import exit
from snowav.framework.read_config import read_config
from snowav.framework.process import process
from snowav.framework.parse import parse
from snowav.framework.figures import figures
from snowav.report.report import report
from snowav.database.database import check_fields,delete,run_metadata,write_csv

class snowav(object):

    def __init__(self, config_file = None, external_logger = None, awsm = None):
        '''
        Initialize with config file, read config file, and run processing,
        results storage, figures, and report with config_file options.

        Args
            config_file: snowav config file
            external_logger: awsm logger
            awsm: awsm class if being run in awsm

        '''

        if awsm is None and os.path.isfile(config_file):
            self.config_file = config_file
            read_config(self)

        elif awsm is not None:
            self.config_file = awsm.filename
            read_config(self, awsm = awsm)

        else:
            raise Exception('No config instance passed, or config file does '
                            'not exist!')

        parse(self)
        run_metadata(self)
        process(self)
        figures(self)

        if self.forecast_flag:
            self._logger.info(' starting forecast processing...')

            # Check for existing fields, delete if necessary
            check_fields(self,
                         self.for_start_date,
                         self.for_end_date,
                         self.plotorder[0],
                         self.for_run_name,
                         'swe_z',
                         forecast=self.for_run_name)
            run_metadata(self, forecast=self.for_run_name)
            process(self, forecast=self.for_run_name)

            # make dict of flags to send to figures? {accumulated: True}
            figures(self)

        if self.report_flag:
            report(self)
