
import os
import pandas as pd
from sqlalchemy import create_engine
from sqlalchemy.orm import sessionmaker
from sqlalchemy import and_
from snowav.database.tables import (Base, RunMetadata, Watershed, Basin,
Results, VariableUnits, Watersheds, Basins)
from sys import exit
from inicheck.tools import get_user_config, check_config
from inicheck.output import generate_config, print_config_report
from inicheck.config import MasterConfig
from snowav.utils.utilities import get_snowav_path
import argparse

def main():
    '''
    Makes simple calls to snowav database and either prints results to the
    terminal or outputs to csv.

    The CoreConfig is ./snowav/config/QueryCore.ini

    snowav_query -f /home/markrobertson/wkspace/code/SNOWAV/scripts/database_query.ini

    '''

    master_config = os.path.join(get_snowav_path(),'snowav/config/QueryCore.ini')

    parser = argparse.ArgumentParser(description="Examine and auto populate")

    parser.add_argument('-f', '--config_file', dest='config_file', type=str,
                        help='Path to configuration file.')

    args = parser.parse_args()

    # read the config file
    mcfg = MasterConfig(path = master_config)
    ucfg = get_user_config(args.config_file, mcfg=mcfg)
    warnings, errors = check_config(ucfg)

    if len(errors) > 0:
        print_config_report(warnings, errors)
        print("Errors in the config file. "
              "See configuration status report above.")
        exit()

    # read in options
    wy = ucfg.cfg['query']['wy']
    basins = ucfg.cfg['query']['basins']
    if type(basins) != list:
        basins = [basins]

    if basins[0] not in Watersheds.watersheds.keys():
        print('First basin in {}: [query] -> basins should be one of '
              '{}'.format(args.config_file,Watersheds.watersheds.keys()))
        exit()

    value = ucfg.cfg['query']['value']
    run_name = ucfg.cfg['query']['run_name']
    start_date = ucfg.cfg['query']['start_date'].to_pydatetime()
    end_date = ucfg.cfg['query']['end_date'].to_pydatetime()
    total = ucfg.cfg['query']['total']
    csv_base_path = ucfg.cfg['query']['csv_base_path']
    output = ucfg.cfg['query']['output']
    db = ucfg.cfg['query']['database']
    print_runs = ucfg.cfg['query']['print_all_runs']

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

if __name__ == '__main__':
    main()
