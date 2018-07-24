
import datetime
import snowav

def package_results(self, df, output, dtime):
    '''
    This function sends process() results to the snowav database.

    Args
        df: results DataFrame
        output: output variable ('swe_z')
        dtime: datetime

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
    basin_id = {'Boise River Basin':1,
                'Featherville':2,
                'Twin Springs':3,
                'Mores Creek':4}

    # Make labels
    if ('z' in output) or ('depth' in output):
        lbl = self.depthlbl
    if ('vol' in output) or ('avail' in output):
        lbl = self.vollbl
    if output == 'density':
        lbl = 'kg/m^3'
    if output == 'coldcont':
        lbl = 'MJ'

    # By sub basin
    for var in df:

        # By elevation band
        for iters,r in enumerate(df[var].values):
            values = {'basin_id': basin_id[var],
                      'date_time': dtime,
                      'proc_time': datetime.datetime.now(),
                      'version': 'snowav'+ snowav.__version__,
                      'variable': output,
                      'var_units': lbl,
                      'value': r,
                      'elevation': str(df[var].index[iters]),
                      'elev_units': self.elevlbl}

            snowav.database.database.insert_results(self.database,values)
