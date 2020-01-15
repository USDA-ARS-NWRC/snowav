
from snowav.database.database import collect
import pandas as pd
from datetime import timedelta

def write_properties(args, values, logger):
    '''
    Write daily total snowpack properties to csv.

    Args
    ------
    args : dict
        dictionary with required inputs, see swi() figure for more information.
    values : str
        snowav database value to query and write to csv

    '''

    datestr = args['end_date'] + timedelta(hours=24)
    datestr = datestr.strftime("%Y%m%d")

    for value in values:
        out = collect(args['connector'], args['plotorder'], args['basins'],
                args['wy_start'], args['end_date'], value, args['run_name'],
                'total','daily')

        if 'vol' in value or 'avail' in value:
            unit = args['vollbl']

        if 'z' in value or value == 'depth':
            unit = args['depthlbl']

        if value == 'density':
            unit = 'kg_m3'

        if value == 'coldcont':
            unit = 'MJ'

        path = '{}{}_timeseries_{}_{}.csv'.format(args['figs_path'], value, datestr, unit.lower())

        logger.info(' saving {}'.format(path))

        out.to_csv(path)
