
import numpy as np
import pandas as pd
from snowav.database.database import insert, query
from datetime import datetime
from snowav.database.tables import RunMetadata, Watershed, Basin, Results, VariableUnits

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
    if output == 'snow_line':
        lbl = 'ft'

    for basin in df:

        for iters, val in enumerate(df[basin].values):
            if np.isnan(val):
                val = None

            else:
                val = float(val)

            values = {'basin_id': int(basins[basin]['basin_id']),
                      'run_id': int(run_id),
                      'date_time': dtime,
                      'variable': output,
                      'variable_id': int(vid[output]),
                      'value': val,
                      'elevation': str(df[basin].index[iters])}

            insert(connector,'Results',values)
