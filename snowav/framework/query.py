

import os
import pandas as pd
from sqlalchemy import create_engine
from sqlalchemy.orm import sessionmaker
from sqlalchemy import and_
from snowav.database.tables import (Base, RunMetadata, Watershed, Basin, Results, VariableUnits, Watersheds, Basins)
from sys import exit
from inicheck.tools import get_user_config, check_config
from inicheck.output import generate_config, print_config_report
from inicheck.config import MasterConfig
from snowav.utils.utilities import get_snowav_path
import argparse
from sys import exit

def query(self):
    '''
    '''

    print('Config option [query] = True, running database query...')

    if type(self.q_basins) != list:
        basins = [self.q_basins]
    else:
        basins = self.q_basins

    if basins[0] not in Watersheds.watersheds.keys():
        print('First basin in {}: [query] -> basins should be one of '
              '{}'.format(args.config_file,Watersheds.watersheds.keys()))
        exit()

    value = self.q_value
    run_name = self.q_run_name
    start_date = self.q_start_date
    end_date = self.q_end_date
    total = self.q_total
    csv_base_path = self.q_csv_base_path
    output = self.q_output
    db = self.q_database
    print_runs = self.q_print_all_runs

    # connect to database
    engine = create_engine(db)
    connection = engine.connect()
    DBSession = sessionmaker(bind=engine)
    session = DBSession()

    # Make df from database
    db_val = pd.DataFrame(columns = basins)

    for bid in basins:

        # Get available run_names for the watershed
        wid = Basins.basins[basins[0]]['watershed_id']

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
                                                      (Results.basin_id == Basins.basins[bid]['basin_id'])))

            df = pd.read_sql(qry.statement, qry.session.connection())

            qry2 = session.query(RunMetadata).join(Results).filter(and_((Results.date_time >= start_date),
                                                      (Results.date_time <= end_date),
                                                      (RunMetadata.run_name.in_(names)) ))
            df2 = pd.read_sql(qry2.statement, qry2.session.connection())

            session.close()

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
                        (Results.basin_id == Basins.basins[bid]['basin_id'])))

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
                      # '{} to {}\n'
                      'value: {}\n\n'
                      # '{}'.format(bid, run_name, start_date, end_date, value, fqry[cols]))
                      '{}'.format(bid, run_name, value, fqry[cols]))

        if output == 'csv':
            fqry.to_csv(csv_file)

    print('query complete, exiting snowav. To process results, set [query] query = False')
    exit()
