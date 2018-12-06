
import numpy as np
import pandas as pd
import datetime
import snowav
from snowav import database
from datetime import datetime
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
            else:
                val = float(val)

            values = {'basin_id': Basins.basins[var]['basin_id'],
                      'run_id':self.runid,
                      'date_time': dtime,
                      'variable': output,
                      'variable_id': self.vid[output],
                      'value': val,
                      'elevation': str(df[var].index[iters])}
                      
            snowav.database.database.insert(self,'Results',values)


def post_process(self, dtime):
    '''
    This function post-processes water year total values for swi and evap.

    '''

    sum_vals = {'swi_vol':'swi_vol_wy','swi_z':'swi_z_wy','evap_z':'evap_z_wy'}

    for val in sum_vals.keys():

        summary = pd.DataFrame(0, index=self.edges, columns=self.masks.keys())

        # Make labels
        if ('z' in val) :
            lbl = self.depthlbl
        if ('vol' in val) :
            lbl = self.vollbl

        for bid in self.plotorder:

            r = database.database.query(self,
                                        datetime(self.wy-1,10,1),
                                        dtime,
                                        self.run_name,
                                        bid,
                                        val)

            for e in self.edges:
                v = r[(r['elevation'] == str(e))]
                summary.loc[e,bid] = np.nansum(v['value'].values)

            for iters, sval in enumerate(summary[bid].values):
                if np.isnan(sval):
                    sval = None
                elif (type(sval) is not str):
                    sval = float(sval)

                values = {'basin_id': Basins.basins[bid]['basin_id'],
                          'run_id':self.runid,
                          'date_time': dtime,
                          'variable': sum_vals[val],
                          'variable_id': self.vid[sum_vals[val]],
                          'value': sval,
                          'elevation': str(summary[bid].index[iters])}

                snowav.database.database.insert(self,'Results',values)

            # add total as the sum
            if np.nansum(summary[bid].values) != np.nan:
                v = float(np.nansum(summary[bid].values))
            else:
                v = np.nansum(summary[bid].values)

            values = {'basin_id': Basins.basins[bid]['basin_id'],
                      'run_id':self.runid,
                      'date_time': dtime,
                      'variable': sum_vals[val],
                      'variable_id': self.vid[sum_vals[val]],
                      'value': v,
                      'elevation': 'total'}

            snowav.database.database.insert(self,'Results',values)
