
import os
from sys import exit
from snowav.methods.read_config import read_config
from snowav.methods.process import process
from snowav.report.report import report
from snowav.database.database import check_fields,delete,run_metadata,write_csv
from snowav.plotting.accumulated import accumulated
from snowav.plotting.basin_total import basin_total
from snowav.plotting.cold_content import cold_content
from snowav.plotting.compare_runs import compare_runs
from snowav.plotting.current_image import current_image
from snowav.plotting.density import density
from snowav.plotting.flt_image_change import flt_image_change
from snowav.plotting.image_change import image_change
from snowav.plotting.inflow import inflow
from snowav.plotting.pixel_swe import pixel_swe
from snowav.plotting.point_values import point_values
from snowav.plotting.precip_depth import precip_depth
from snowav.plotting.precip_validate import precip_validate
from snowav.plotting.stn_validate import stn_validate
from snowav.plotting.subbasins import subbasins
from snowav.plotting.swe_change import swe_change
from snowav.plotting.swe_volume import swe_volume

class SNOWAV(object):

    def __init__(self, config_file = None, external_logger = None, awsm = None):
        '''
        Initialize snowav class with config file, read config file, and
        run processing, results storage, figures, and reports as specified in
        the config file.

        Args
            config_file: snowav config file for independent snowav processing
            external_logger: awsm logger
            awsm: awsm class if being run with awsm. Currently, this only uses
                the config file and run directory from the awsm class

        '''

        if awsm is None and os.path.isfile(config_file):
            self.config_file = config_file
            read_config(self)

        elif awsm is not None:
            self.config_file = awsm.filename
            read_config(self, awsm = awsm)

        else:
            print('No config instance passed, or config file does not exist!')
            return

        if self.report_only is True:
            print('Config option [Results] -> report_only=True, ' +
                  'report generated with existing figures.')
            report(self)
            exit()

        # These can be made from pulling from the database
        if self.figures_only is True:

            if self.compare_runs_flag is True:
                compare_runs(self)

            if self.subbasins_flag is True:
                subbasins(self)

            if self.basin_total_flag is True:
                basin_total(self)

            if self.pixel_swe_flag is True:
                pixel_swe(self)

            if self.density_flag is True:
                density(self)

            if self.compare_runs_flag is True:
                compare_runs(self)

            if self.basin_detail_flag is True:
                basin_detail(self)


        # Otherwise, start the processing steps
        else:
            # Do any values in this date already already exist?
            check_fields(self, self.start_date, self.end_date,
                         self.plotorder[0], self.run_name, 'swe_z')

            # Process results and put on the database
            if (self.write_db is False):
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
                # process(self)

            elif (self.pflag is True) and (self.write_db is True):
                print('There are existing fields on the database between '
                      '{} and {} with run_name={}, and config file option '
                      '[results] -> overwrite=True, '
                      'OVERWRITING RESULTS!!!'.format(self.start_date.date(),
                                                      self.end_date.date(),
                                                      self.run_name))

                # Delete existing fields
                for bid in self.plotorder:
                    delete(self, self.start_date, self.end_date,
                           bid, self.run_name)

                # Process and put on database
                run_metadata(self)
                process(self)

            else:
                # Process and put on database
                run_metadata(self)
                process(self)

            # Write out variables from database to csv if desired
            if self.write_csv_flag is True:
                write_csv(self)

            if self.subbasins_flag is True:
                subbasins(self)

            if self.accumulated_flag is True:
                accumulated(self)

            if self.current_image_flag is True:
                current_image(self)

            if self.image_change_flag is True:
                image_change(self)

            if self.cold_content_flag is True:
                cold_content(self)

            if self.swe_volume_flag is True:
                swe_volume(self)

            if self.swe_change_flag is True:
                swe_change(self)

            if self.basin_total_flag is True:
                basin_total(self)

            if self.pixel_swe_flag is True:
                pixel_swe(self)

            if self.density_flag is True:
                density(self)

            if self.stn_validate_flag is True:
                stn_validate(self)

            if self.inflow_flag is True:
                inflow(self)

            if self.compare_runs_flag is True:
                compare_runs(self)

            if self.inflow_flag is True:
                inflow(self)

            if self.precip_depth_flag is True:
                precip_depth(self)

            # Write out current model SWE values at snow course locations
            if self.write_stn_csv_flag is True:
                point_values(self.outputs['swe_z'][-1],
                             self.stns_csv,
                             (self.snow_x, self.snow_y),
                             '{}model_pixel_swe_{}.csv'.format(self.figs_path,
                             self.end_date.date().strftime("%Y%m%d")))

            if self.precip_validate_flag is True:
                precip_validate(self)

            if self.basin_detail_flag is True:
                basin_detail(self)

            # Make flight difference figure in options in config file
            if self.flt_flag is True:
                flt_image_change(self)

            # Create pdf report
            if self.report_flag is True:
                report(self)

        # Close the database session
        self.session.close()
