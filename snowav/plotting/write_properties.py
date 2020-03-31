from datetime import datetime
import os

from snowav.database.database import collect


def write_properties(end_date, cnx, plotorder, basins, wy_start, run_name,
                     figs_path, values, vollbl='TAF', logger=None):
    """ Write daily total snowpack properties to csv.

    Args
    ------
    end_date {str}: end_date
    cnx {str}: database connector
    plotorder {list}: basins list
    basins {dict}: basins dict
    wy_start {str}: YYYY1001
    run_name {str}: snowav run_name
    figs_path {str}: path to save files
    values {list}: list of snowav values
    vollbl {str}: volume label
    depthlbl {str}: depth label
    logger {class}: logger
    """

    datestr = end_date.strftime("%Y%m%d")
    now_str = datetime.now().date().strftime("%Y-%m-%d")
    date_col = 'Date generated: {}'.format(now_str)
    unit = vollbl

    for value in values:
        out = collect(cnx, plotorder, basins, wy_start, end_date, value,
                      run_name, 'total', 'daily')

        out.index = out.index.date

        # setting index to date strips the index name
        out.index.name = 'date'

        if value.lower() == 'swe_vol':
            value_line = ('Snow Water Equivalent (SWE) volume in thousands ' +
                          'of acre-feet (TAF)')
        elif value.lower() == 'swi_vol':
            value_line = ('Surface Water Input (SWI) volume in thousands ' +
                          'of acre-feet (TAF)')
        else:
            if logger is not None:
                logger.warning(" Value types other than swe, swi, and "
                               "their derivates are not supported")
            return

        headers = ['USDA Agicultural Research Service Snowpack Summary Data',
                   value_line,
                   'Data provided are daily model results from the iSnobal model',
                   'First column is the date of model result',
                   'Second column is the total basin volume',
                   'Additional columns are the subbasins in the watershed',
                   date_col,
                   'Valid until next reports are generated',
                   'Contact: Scott Havens <scott.havens@usda.gov>'
                   '\n']

        filename = '{}_timeseries_{}_{}.csv'.format(value, datestr, unit.lower())
        path = os.path.join(os.path.abspath(figs_path), filename)

        if os.path.isfile(path):
            os.remove(path)

        with open(path, mode='w', encoding='utf-8') as f:
            f.write('\n'.join(headers))

        out.to_csv(path, encoding='utf-8', mode='a')

        if logger is not None:
            logger.info(' Saved: {}'.format(path))
