
import numpy as np
import datetime
import snowav
from snowav.database.tables import Basins
from snowav.database.tables import RunMetadata, Watershed, Basin, Results, VariableUnits, Watersheds, Basins

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
            if np.isnan(val):
                val = None

            values = {'basin_id': Basins.basins[var]['basin_id'],
                      'run_id':self.runid,
                      'date_time': dtime,
                      'variable': output,
                      'variable_id': self.vid[output],
                      'value': val,
                      'elevation': str(df[var].index[iters])}

            snowav.database.database.insert(self,'Results',values)
