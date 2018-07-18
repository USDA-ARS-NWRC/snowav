
import datetime
import snowav

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
    # This is not efficient
    basin_id = {'Boise River Basin':1,'Featherville':2,'Twin Springs':3,'Mores Creek':4}

    # This dict should be created in process()
    dataframes = {'swe':self.swe_volume}
    k = list(dataframes.keys())
    # k[0] is 'swe'
    # var is like 'Boise River Basin'
    # dataframes[k[0]][var] is self.swe_volume['Boise River Basin']

    # By basin
    for var in list(dataframes[k[0]]):

        # By elevation band
        for iters,r in enumerate(dataframes[k[0]][var].values):

            values = {'basin_id': basin_id[var],
                      'date_time': self.dateTo,
                      'proc_time': datetime.datetime.now(),
                      'version': 'snowav'+ snowav.__version__,
                      'variable': 'swe',
                      'var_units': self.vollbl,
                      'value': r,
                      'elevation': str(dataframes[k[0]][var].index[iters]),
                      'elev_units': self.elevlbl}

            snowav.database.database.insert_results(self.database,values)

        # All elevations
        values = {'basin_id': basin_id[var],
                  'date_time': self.dateTo,
                  'proc_time': datetime.datetime.now(),
                  'version': 'snowav'+ snowav.__version__,
                  'variable': 'swe',
                  'var_units': self.vollbl,
                  'value': sum(dataframes[k[0]][var].values),
                  'elevation': 'total',
                  'elev_units': self.elevlbl}

        snowav.database.database.insert_results(self.database,values)


    # self.values = values
