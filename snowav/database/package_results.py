
import numpy as np
import pandas as pd
from snowav.database.database import insert, query
from datetime import datetime
from snowav.database.tables import Basins
from snowav.database.tables import RunMetadata, Watershed, Basin, Results, VariableUnits, Watersheds, Basins

def package(connector, lbls, basins, df, run_id, vid, output, dtime, run_name):
    '''
    Put process() results on the database.

    Args
    ------
    connector : str
        database connector
    lbls : dict
    basins : dict
    run_id : int
    vid : int
        variable id
    df : DataFrame
    output : str
        output variable (i.e. 'swe_z')
    dtime: datetime

    '''

    if ('z' in output) or ('depth' in output):
        lbl = lbls['depthlbl']
    if ('vol' in output) or ('avail' in output):
        lbl = lbls['vollbl']
    if output == 'density':
        lbl = 'kg/m^3'
    if output == 'coldcont':
        lbl = 'MJ'

    for basin in df:

        for iters, val in enumerate(df[basin].values):
            if np.isnan(val):
                val = None

            else:
                val = float(val)

            values = {'basin_id': basins[basin]['basin_id'],
                      'run_id': run_id,
                      'date_time': dtime,
                      'variable': output,
                      'variable_id': vid[output],
                      'value': val,
                      'elevation': str(df[basin].index[iters])}

            insert(connector,'Results',values)

def post_process(self, dtime):
    '''
    Post-processes water year total values for swi and evap.

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
            r = query(self.connector, datetime(self.wy-1,10,1), dtime,
                      self.run_name, self.basins, bid, val)

            for e in self.edges:
                v = r[(r['elevation'] == str(e))]
                summary.loc[e,bid] = np.nansum(v['value'].values)

            for iters, sval in enumerate(summary[bid].values):
                if np.isnan(sval):
                    sval = None
                elif (type(sval) is not str):
                    sval = float(sval)

                values = {'basin_id': Basins.basins[bid]['basin_id'],
                          'run_id':self.run_id,
                          'date_time': dtime,
                          'variable': sum_vals[val],
                          'variable_id': self.vid[sum_vals[val]],
                          'value': sval,
                          'elevation': str(summary[bid].index[iters])}

                insert(self.connector,'Results',values)

            # add total as the sum
            if np.nansum(summary[bid].values) != np.nan:
                v = float(np.nansum(summary[bid].values))
            else:
                v = np.nansum(summary[bid].values)

            values = {'basin_id': Basins.basins[bid]['basin_id'],
                      'run_id':self.run_id,
                      'date_time': dtime,
                      'variable': sum_vals[val],
                      'variable_id': self.vid[sum_vals[val]],
                      'value': v,
                      'elevation': 'total'}

            insert(self.connector,'Results',values)
