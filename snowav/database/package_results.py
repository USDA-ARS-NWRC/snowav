
import datetime
import snowav
from snowav.database.tables import BASINS
from snowav.database.tables import BasinMetadata, Base, Results, RunMetadata, BASINS

def package(self, df, output, dtime):
    '''
    This function sends process() results to the database.

    Args
        df: results DataFrame
        output: output variable ('swe_z')
        dtime: datetime

    '''

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
        for iters,val in enumerate(df[var].values):
            values = {'basin_id': BASINS.basins[var]['basin_id'],
                      'run_id':self.runid,
                      'date_time': dtime,
                      'variable': output,
                      'value': val,
                      'elevation': str(df[var].index[iters])}
            # r = Results(values)

            snowav.database.database.insert(self,'Results',values)
