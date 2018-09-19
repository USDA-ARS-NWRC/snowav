
import sys
import snowav
import argparse
import datetime
import pandas as pd
import os
from sqlalchemy.ext.declarative import declarative_base
from sqlalchemy import create_engine
from sqlalchemy.orm import sessionmaker
from snowav.database.tables import Base, Results, BASINS
from sqlalchemy import and_
from snowav.utils.utilities import get_snowav_path

def run():

    parser = argparse.ArgumentParser(description="Examine and auto populate")

    parser.add_argument('-f', '--config_file', dest='config_file', type=str,
                        help='Path to snowav configuration file.')

    parser.add_argument('-b', '--bid', dest='bid', type=str,
                        help='Basin id, options in snowav/database/tables.py, '
                        'for a database query of existing fields. If -wy and '
                        '-b are passed, will check database.')

    parser.add_argument('-wy', '--wy', dest='wy', type=int,
                        help='Water year, for a database query of existing '
                        'fields. If -wy and -b are passed, will check database.')

    parser.add_argument('-db', '--db', dest='db', type=str,
                        help='Path to snowav database, for a database query of '
                        'existing fields. If -wy and -b are passed, will '
                        'check database.')

    args = parser.parse_args()

    # If config file is passed, do standard snowav processing, figs, and report
    if args.config_file:
        snowav.framework.framework.SNOWAV(config_file = args.config_file)

    # Run a database query for existing fields if desired
    if args.bid and args.wy:

        # Making this a default for now
        value = 'swe_vol'

        # If no database is passed, try the default data location
        if not args.db:
            sbase = get_snowav_path()
            sbase = 'sqlite:///' + sbase
            args.db = sbase + '/snowav/data/model_results.db'

        start_date = datetime.datetime(args.wy - 1,10,1)
        end_date = datetime.datetime(args.wy,9,30)

        try:
            # Connect to database and query
            engine = create_engine(args.db)
            connection = engine.connect()
            DBSession = sessionmaker(bind=engine)
            session = DBSession()
            qry = session.query(Results).filter(and_((Results.date_time >= start_date),
                                                      (Results.date_time <= end_date),
                                                      (Results.basin_id == BASINS.basins[args.bid]['basin_id'])))

            df = pd.read_sql(qry.statement, qry.session.connection())

            session.close()

            print('On database: {},\nAvailable runs for {}: {},\nDate range: {}, {}'.format(
            args.db,args.bid,df.run_name.unique(),df['date_time'].min(),df['date_time'].max()))

        except:
            print('Failed connecting to database {} for query.'.format(args.db))

if __name__ == '__main__':
    run()
