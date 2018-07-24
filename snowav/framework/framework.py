
import os
import snowav
from snowav import methods
from snowav import database

class SNOWAV(object):

    def __init__(self,config_file = None):
        '''
        Initialize snowav class with config file, read config file, and
        run processing, results storage, figures, and reports as specified in
        the config file.

        To-dos:
        - expand database.database.check_fields functionality for more robust
            field checking

        '''

        if os.path.isfile(config_file):
            self.config_file = config_file

        else:
            self.config_file = config_file
            print('SNOWAV config file does not exist!')
            return

        # Read config file
        methods.read_config.read_config(self)

        # Do any values in this date already already exist?
        flag = database.database.check_fields(self.database, self.start_date,
                                              self.end_date, 'swe')

        # If data does not exist, process and put it on the database
        if (self.db_overwrite_flag is True) or (flag is not True):

            # Process
            snowav.methods.process.process(self)

            # Save data
            if self.location == 'database':

                # Put into format for Results table on database
                # database.package_results.package_results(self)
                print('finished onto database')

                # Insert into database
                # snowav.database.database.insert_results(self,values)

            if self.location == 'csv':
                print('Make a csv-saver...')


        # If data does already exist on the database
        if (flag is True) and (self.db_overwrite_flag is not True):
            print('Database values in the date range '
                  + '%s to %s already exist.\n'%(self.start_date,self.end_date)
                  + 'Config option [results] -> overwrite is set to False, '
                  + 'if you wish to overwrite, set to True.')


        if (flag is True) and (self.db_overwrite_flag is True):
            print('Database values in the date range '
                  + '%s to %s already exist.\n'%(self.start_date,self.end_date)
                  + 'WARNING: config option [results] -> overwrite is set to True, '
                  + 'database values will be overwritten...')

        # Plots
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
