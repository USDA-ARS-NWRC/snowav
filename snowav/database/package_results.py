
import datetime

def package_results(self):
    '''
    values = {'basin_id': 1,
              'date_time': datetime.datetime.now(),
              'proc_time': datetime.datetime.now(),
              'version': '1',
              'variable': 'swe',
              'var_units': 'in',
              'value': 2.0,
              'elevation': 'total',
              'elev_units': 'ft'}

    '''

    # Should make a check for table field values...
    values = {'basin_id': 2,
              'date_time': datetime.datetime.now(),
              'proc_time': datetime.datetime.now(),
              'version': '1',
              'variable': 'swi',
              'var_units': 'in',
              'value': 2.0,
              'elevation': 'total',
              'elev_units': 'ft'}

    self.values = values
