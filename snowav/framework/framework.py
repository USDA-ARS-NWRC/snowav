
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
        flag = database.database.check_fields(self.database,
                                              self.start_date,
                                              self.end_date,
                                              self.plotorder[0],
                                              'swe_z')

        # Process results and put on the database
        if (flag is True) and (self.db_overwrite_flag is False):
            print('There are existing fields on the database for this time '
            'period, and the config file option [results] -> overwrite is set '
            'to False, skipping processing...')
        elif (flag is True) and (self.db_overwrite_flag is True):
            print('There are existing fields on the database for this time '
            'period, and the config file option [results] -> overwrite is set '
            'to True, OVERWRITING RESULTS')

            for bid in self.plotorder:
                database.database.delete(self.database,
                                         self.start_date,
                                         self.end_date,
                                         bid)
            snowav.methods.process.process(self)

        else:
            snowav.methods.process.process(self)

        # Plots
        snowav.plotting.accumulated.accumulated(self)
        snowav.plotting.current_image.current_image(self)
        snowav.plotting.state_by_elev.state_by_elev(self)
        snowav.plotting.image_change.image_change(self)
        snowav.plotting.swe_change.swe_change(self)
        snowav.plotting.basin_total.basin_total(self)
        snowav.plotting.pixel_swe.pixel_swe(self)
        # snowav.plotting.density.density(self)
        # snowav.plotting.water_balance.water_balance(self)
        snowav.plotting.stn_validate.stn_validate(self)

        # This may need its own query as well
        if self.flt_flag is True:
            snowav.plotting.flt_image_change.flt_image_change(self)

        # Run reporting module if desired
        if self.report_flag is True:
            snowav.report.report.report(self)
