
import os
import pandas as pd
from sqlalchemy import and_
from snowav.database.database import connect, make_session
from snowav.database.tables import RunMetadata, Results, VariableUnits
from sys import exit

def query(self):
    '''
    Query snow server snowav database. Requires correct database login
    information in the config file [database] section.

    See README.md and CoreConfig.ini for details and more information.

    '''

    print('\nConfig option [query] = True, running database query...')

    if ((self.q_basins is None) or (self.q_run_name is None) or
        (self.q_start_date is None) or (self.q_end_date is None)):
        raise Exception('Required config [query] fields missing, '
                        'see CoreConfig.ini')

    if type(self.q_basins) != list:
        basins = [self.q_basins]
    else:
        basins = self.q_basins

    bdict, cnx, out = connect(sqlite = self.sqlite, sql = self.mysql,
                              plotorder = basins, user = self.db_user,
                              password = self.db_password, host = self.db_host,
                              port = self.db_port, convert = self.db_convert)

    value = self.q_value
    run_name = self.q_run_name
    start_date = self.q_start_date
    end_date = self.q_end_date
    total = self.q_total
    csv_base_path = self.q_csv_base_path
    output = self.q_output
    db = self.q_database
    print_runs = self.q_print_all_runs

    # Make df from database
    db_val = pd.DataFrame(columns = basins)
    session = make_session(self.q_database)

    for bid in basins:

        # Get available run_names for the watershed
        wid = int(bdict[basins[0]]['watershed_id'])

        # print out all available runs for the time period if desired
        if print_runs is True:

            q = session.query(RunMetadata).join(Results).filter(and_(RunMetadata.watershed_id == wid),
                                                                    (Results.date_time >= start_date),
                                                                    (Results.date_time <= end_date))
            df = pd.read_sql(q.statement, q.session.connection())
            names = df.run_name.unique()

            # run_name and run_id
            rid = {}
            for name in names:
                rid[name] = df[df['run_name'] == name]['run_id'].values

            qry = session.query(Results).join(RunMetadata).filter(and_((Results.date_time >= start_date),
                                                      (Results.date_time <= end_date),
                                                      (RunMetadata.run_name.in_(names)),
                                                      (RunMetadata.watershed_id == wid),
                                                      (Results.basin_id == int(bdict[bid]['basin_id']))))

            df = pd.read_sql(qry.statement, qry.session.connection())

            qry2 = session.query(RunMetadata).join(Results).filter(and_((Results.date_time >= start_date),
                                                      (Results.date_time <= end_date),
                                                      (RunMetadata.run_name.in_(names)) ))
            df2 = pd.read_sql(qry2.statement, qry2.session.connection())

            dash = '-'*100
            print('\n{}'.format(dash))
            print('On snowav database for {}:'.format(bid))
            print(dash)

            for name in rid:
                try:
                    stime = df[df['run_id'].isin(rid[name])]['date_time'].min().date().strftime("%Y-%-m-%-d")
                except:
                    stime = 'nan'

                try:
                    etime = df[df['run_id'].isin(rid[name])]['date_time'].max().date().strftime("%Y-%-m-%-d")
                except:
                    etime  = 'nan'

                print('{:<25s}{:<25s}{:<25s}'.format(name, stime, etime ))

            print('{}\n'.format(dash))

        bid_name = bid.replace(" ", "_")
        csv_file = '{}{}_{}_{}.csv'.format(csv_base_path,bid_name,run_name,value)
        qry = session.query(Results).join(RunMetadata).filter(and_(
                        (Results.date_time >= start_date),
                        (Results.date_time <= end_date),
                        (RunMetadata.run_name == run_name),
                        (RunMetadata.watershed_id == wid),
                        (Results.variable == value),
                        (Results.basin_id == int(bdict[bid]['basin_id']))))

        df = pd.read_sql(qry.statement, qry.session.connection())

        # Use the first call to get the units
        qry2 = session.query(VariableUnits).filter(VariableUnits.id == int(df['variable_id'].values[0]))
        df2 = pd.read_sql(qry2.statement, qry2.session.connection())
        df['units'] = df2['unit'].values[0]

        if total is True:
            fqry = df[df['elevation'] == 'total']

        else:
            fqry = df

        fqry = fqry.sort_values(by=['date_time'])

        db_val.loc['total',bid] = fqry['value'].values

        cols = ['date_time','value','units','elevation']

        if output == 'print':
            with pd.option_context('display.max_rows', None, 'display.max_columns', None):
                print('\n\nsnowav database values \n'
                      'basin: {}\n'
                      'run_name: {}\n'
                      'value: {}\n\n'
                      '{}'.format(bid, run_name, value, fqry[cols]))

        if output == 'csv':
            if total:
                fqry.to_csv(csv_file)
            else:
                for e in fqry.elevation.unique().tolist():
                    banded = fqry.copy()
                    banded = banded.set_index('date_time')
                    banded = banded[banded['elevation'] == e]
                    banded = pd.DataFrame(banded['value'].rename('{} {} [ft]'.format(value, e)))
                    banded.to_csv('{}_{}.csv'.format(csv_file.split('.csv')[0],e))

    session.close()
    print('\nQuery complete, exiting snowav. To process results, set [query] query: False\n')
    exit()
