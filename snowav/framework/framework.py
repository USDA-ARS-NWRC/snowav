
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
        # results from the database. This will also attempt to make the figures
        # if they don't already exist
        if self.plot_flag is True:
            # These currently require loading nc files and read_config()
            # snowav.plotting.accumulated.accumulated(self)
            # snowav.plotting.current_image.current_image(self)
            # snowav.plotting.pixel_swe.pixel_swe(self)
            # snowav.plotting.image_change.image_change(self)
            # snowav.plotting.swe_change.swe_change(self)
            # snowav.plotting.precip_depth.precip_depth(self)

            if (not os.path.isfile('{}swe_elev_{}.png'
                                    .format(self.figs_path,self.name_append))):
                snowav.plotting.state_by_elev.state_by_elev(self)

            # There are two, we just check the first
            if (not os.path.isfile('{}basin_total_{}.png'
                                    .format(self.figs_path,self.name_append))):
                snowav.plotting.basin_total.basin_total(self)

            # There are several, we just check the first
            if (not os.path.isfile('{}density_subs_{}.png'
                                    .format(self.figs_path,self.name_append))):
                snowav.plotting.density.density(self)

            if (not os.path.isfile('{}validation_{}.png'
                                    .format(self.figs_path,self.name_append))):
                snowav.plotting.stn_validate.stn_validate(self)

            # Create pdf report
            if self.report_flag is True:
                snowav.report.report.report(self)

            # This one still needs editing
            # snowav.plotting.compare_runs.compare_runs(self)

        # Otherwise, start the processing steps
        else:
            # Do any values in this date already already exist?
            database.database.check_fields(self,
                                           self.start_date,
                                           self.end_date,
                                           self.plotorder[0],
                                           self.run_name,
                                           'swe_z')

            # Process results and put on the database
            if (self.pflag is True) and (self.write_db is False):
                print('There are existing fields on the database between '
                      '{} and {} with run_name={}, and config file option '
                      '[results] -> overwrite=False\n'
                      'Using iSnobal outputs specified in run_dirs for '
                      'spatial figures, loading processed results '
                      'from the database'.format(self.start_date.date(),
                                                      self.end_date.date(),
                                                      self.run_name))

                # We still send to process() in order to get spatial sums for SWI,
                # precip, etc., but do not insert onto database
                snowav.methods.process.process(self)

            elif (self.pflag is True) and (self.write_db is True):
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
            snowav.plotting.density.density(self)
            snowav.plotting.stn_validate.stn_validate(self)
            snowav.plotting.precip_validate.precip_validate(self)

            # The SWI, precip, and rain plot requires process()
            if ((self.exclude_figs is not None) and
                ('PRECIP_DEPTH' not in self.exclude_figs)):
                snowav.plotting.precip_depth.precip_depth(self)

            # Make flight difference figure in options in config file
            if self.flt_flag is True:
                snowav.plotting.flt_image_change.flt_image_change(self)

            # Create pdf report
            if self.report_flag is True:
                snowav.report.report.report(self)

        # Close the database session
        self.session.close()
