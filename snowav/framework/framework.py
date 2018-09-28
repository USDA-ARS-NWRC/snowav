
import os
import snowav
from snowav import methods
from snowav import database

class SNOWAV(object):

    def __init__(self, config_file = None, external_logger = None, awsm = None):
        '''
        Initialize snowav class with config file, read config file, and
        run processing, results storage, figures, and reports as specified in
        the config file.

        Args:
            config_file: snowav config file for independent snowav processing
            external_logger: awsm logger
            awsm: awsm class if being run with awsm. Currently, this only uses
                the config file and run directory from the awsm class

        '''

        if os.path.isfile(config_file):
            self.config_file = config_file

        else:
            self.config_file = config_file
            print('SNOWAV config file does not exist!')
            return

        # Read config file
        if external_logger is not None:
            methods.read_config.read_config(self,
                                            external_logger = external_logger,
                                            awsm = awsm)
        else:
            methods.read_config.read_config(self)

        # If user just wants plots from database, skip processing and just pull
        # results from the database
        if self.plot_flag is True:
            snowav.plotting.compare_runs.compare_runs(self)

        # Otherwise, start the processing steps
        else:
            # Do any values in this date already already exist?
            flag = database.database.check_fields(self,
                                                  self.start_date,
                                                  self.end_date,
                                                  self.plotorder[0],
                                                  self.run_name,
                                                  'swe_z')

            # Process results and put on the database
            if (flag is True) and (self.db_overwrite_flag is False):
                print('There are existing fields on the database between '
                      '{} and {} with run_name={}, and config file option '
                      '[results] -> overwrite=False, '
                      'skipping processing...'.format(self.start_date.date(),
                                                      self.end_date.date(),
                                                      self.run_name))

            elif (flag is True) and (self.db_overwrite_flag is True):
                print('There are existing fields on the database between '
                      '{} and {} with run_name={}, and config file option '
                      '[results] -> overwrite=True, '
                      'OVERWRITING RESULTS!!!'.format(self.start_date.date(),
                                                      self.end_date.date(),
                                                      self.run_name))

                # Delete existing fields
                for bid in self.plotorder:
                    database.database.delete(self, self.start_date,
                                             self.end_date, bid, self.run_name)

                # Process and put on database
                database.database.run_metadata(self)
                snowav.methods.process.process(self)

            else:
                # Process and put on database
                database.database.run_metadata(self)
                snowav.methods.process.process(self)

            # Write out variables from database to csv if desired
            if self.write_csv_flag is True:
                database.database.write_csv(self)

            # Plots
            snowav.plotting.accumulated.accumulated(self)
            snowav.plotting.current_image.current_image(self)
            snowav.plotting.state_by_elev.state_by_elev(self)
            snowav.plotting.image_change.image_change(self)
            snowav.plotting.swe_change.swe_change(self)
            snowav.plotting.basin_total.basin_total(self)
            snowav.plotting.pixel_swe.pixel_swe(self)
            snowav.plotting.stn_validate.stn_validate(self)

            # The SWI, precip, and rain plot requires process() 
            if self.plot_flag is not True:
                snowav.plotting.precip_depth.precip_depth(self)

            # Make flight difference figure in options in config file
            if self.flt_flag is True:
                snowav.plotting.flt_image_change.flt_image_change(self)

            # Create pdf report
            if self.report_flag is True:
                snowav.report.report.report(self)

        # Close the database session
        self.session.close()
