
import sys
import snowav
import argparse
import datetime
import pandas as pd
import os
from sqlalchemy.ext.declarative import declarative_base
from sqlalchemy import create_engine
from sqlalchemy.orm import sessionmaker
from snowav.database.tables import Base, Results, Basins, RunMetadata, Watersheds
from sqlalchemy import and_
from snowav.utils.utilities import get_snowav_path
import urllib.parse

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

    parser.add_argument('-create', '--create', dest='create', type=str,
                        help='Flag for initial database creation for docker.')

    args = parser.parse_args()

    #########################################################################
    #                       SNOWAV run                                      #
    #########################################################################

    # If config file is passed, do standard snowav processing, figs, and report
    if args.config_file:
        snowav.framework.framework.SNOWAV(config_file = args.config_file)

    #########################################################################
    #                       Database creationm                              #
    #########################################################################

    if args.create:
        print('create')

        fp = os.path.abspath(args.create + '/model_results/')

        if not os.path.exists(fp):
            os.makedirs(fp)

        self.database = 'sqlite:///{}'.format(database)

        try:
            print('Creating and using {}'.format(self.database))

            # Make database connection for duration of snowav processing
            engine = create_engine(self.database)
            Base.metadata.create_all(engine)
            DBSession = sessionmaker(bind=engine)
            self.session = DBSession()

            # Initialize basin metadata for all basin defitions
            for basin in BASINS.basins:
                val = BasinMetadata(basin_id = BASINS.basins[basin]['basin_id'],
                                     basin_name = BASINS.basins[basin]['basin_name'],
                                     state = BASINS.basins[basin]['state'])
                self.session.add(val)
                self.session.commit()

        except:

            print('Database {} creation failed...'.format(self.database))

    #########################################################################
    #                       Database query                                  #
    #########################################################################

    # Database query for existing fields
    if args.bid and args.wy:

        # Making this a default for now
        value = 'swe_vol'

        # If no database is passed, try the default data location
        if not args.db:
            sbase = get_snowav_path()
            sbase = 'sqlite:///' + sbase
            args.db = sbase + '/snowav/data/model_results.db'

        start_date = datetime.datetime(args.wy-1,10,1)
        end_date = datetime.datetime(args.wy,9,30)

        try:
            engine = create_engine(args.db)
            connection = engine.connect()
            DBSession = sessionmaker(bind=engine)
            session = DBSession()

            # Get available run_names for the watershed
            wid =  Basins.basins[args.bid]['watershed_id']
            q = session.query(RunMetadata).filter(RunMetadata.watershed_id == wid)
            df = pd.read_sql(q.statement, q.session.connection())
            names = df.run_name.unique()

            qry = session.query(Results).join(RunMetadata).filter(and_((Results.date_time >= start_date),
                                                      (Results.date_time <= end_date),
                                                      (RunMetadata.run_name.in_(names)),
                                                      (Results.basin_id == Basins.basins[args.bid]['basin_id'])))

            df = pd.read_sql(qry.statement, qry.session.connection())

            session.close()

            print('On database: {},\nAvailable runs for {}: {},\nDate range: {}, {}'.format(
            args.db,args.bid,names,df['date_time'].min(),df['date_time'].max()))

        except:
            print('Failed connecting to database {} for query.'.format(args.db))

if __name__ == '__main__':
    run()
