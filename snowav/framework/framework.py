
import os
import snowav
from snowav import methods

class SNOWAV(object):

    def __init__(self,config_file = None):
        '''
        Initialize snowav class with config file, read config file, and
        run processing, results storage, figures, and reports as specified in
        the config file.

        '''

        if os.path.isfile(config_file):
            self.config_file = config_file

        else:
            self.config_file = config_file
            print('SNOWAV config file does not exist!')
            return

        # Read config file
        methods.read_config.read_config(self)
        # include check in read_config

        # Check if data already exists on database (regardless of config option)
        self.database_flag = False
        if self.database_flag is not True:
            # Finish pre-processing by loading snow.nc, etc...

            # Process
            snowav.methods.process.process(self)

            # Save
            if self.location == 'database':
                # Package results into values dict

                # Insert into database
                print('not done yet')
                # snowav.database.database.insert_results(self,values)

            elif self.location == 'csv':
                print('Make a csv-saver...')

            else:
                print('not sure what to do here yet')

        # If data does already exist on the database
        if self.database_flag is True:
            print('Some nice question about should we overwrite?')

            # Design some queries to get summary data we need
            # Should be able to pull these fields from snow, but write them
            # out explicitly here so that we can call outside of this more easily
            # later on if necessary
            # snowav.database.database.query_basin_value(self.database, start_date,
            #                                            end_date, value)

        snowav.plotting.accumulated.accumulated(self)
        snowav.plotting.current_image.current_image(self)
        snowav.plotting.state_by_elev.state_by_elev(self)
        snowav.plotting.image_change.image_change(self)
        snowav.plotting.swe_change.swe_change(self)
        snowav.plotting.basin_total.basin_total(self)
        snowav.plotting.pixel_swe.pixel_swe(self)
        snowav.plotting.density.density(self)
        snowav.plotting.water_balance.water_balance(self)
        snowav.plotting.stn_validate.stn_validate(self)

        # This may need its own query as well
        if self.flt_flag is True:
            snowav.plotting.flt_image_change.flt_image_change(self)

        # Run reporting module if desired
        if self.report_flag is True:
            snowav.report.report.report(self)
