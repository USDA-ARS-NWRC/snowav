

import sys
import snowav
import argparse
import datetime
import pandas as pd
import os
from sqlalchemy.ext.declarative import declarative_base
from sqlalchemy import create_engine
from sqlalchemy.orm import sessionmaker
from sqlalchemy import and_
from snowav.utils.utilities import get_snowav_path
import urllib.parse
import mysql.connector
from snowav.database.tables import Base, RunMetadata, Watershed, Basin, Results, VariableUnits, Watersheds, Basins
from snowav import database
from sys import exit

def run(config_instance = None):

    if config_instance is None:
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

        parser.add_argument('-user', '--user', dest = 'user', type=str,
                            help='Flag for initial database creation for docker.')

        parser.add_argument('-pwd', '--pwd', dest = 'pwd', type=str,
                            help='Flag for initial database creation for docker.')

        parser.add_argument('-host', '--host', dest = 'host', type=str,
                            help='Flag for initial database creation for docker.')

        parser.add_argument('-port', '--port', dest = 'port', type=str,
                            help='Flag for initial database creation for docker.')

        args = parser.parse_args()

        #########################################################################
        #                       SNOWAV run                                      #
        #########################################################################

        # If config file is passed, do standard snowav processing, figs, and report
        if args.config_file:
            snowav.framework.framework.SNOWAV(config_file = args.config_file)

    else:
        args = None
        snowav.framework.framework.SNOWAV(awsm = config_instance)
        exit()
    #########################################################################
    #                       Database creation                               #
    #########################################################################
    # snowav -create model_r -user mark -pwd whatdystm?1 -host 127.0.0.1

    if args.create:
        # First, check if database exists, if not, make it
        cnx = mysql.connector.connect(user=args.user,
                                      password=args.pwd,
                                      host=args.host,
                                      port=args.port,
                                      auth_plugin='mysql_native_password')
        cursor = cnx.cursor()

        # Check if database exists, create if necessary
        query = ("SHOW DATABASES")
        cursor.execute(query)
        dbs = cursor.fetchall()
        dbs = [i[0].decode("utf-8") for i in dbs]
        print(dbs)

        db_engine = 'mysql+mysqlconnector://{}:{}@{}:{}/{}'.format(args.user,
                                                                args.pwd,
                                                                args.host,
                                                                args.port,
                                                                args.create)

        # If the database doesn't exist, create it, otherwise connect
        if args.create not in dbs:
            print('Specified mysql database {} does not exist, it is being '
                  'created...'.format(args.create))
            query = ("CREATE DATABASE {};".format(args.create))
            cursor.execute(query)
            cursor.close()
            cnx.close()

            database.database.create_tables(url = db_engine)

        else:
            try:
                engine = create_engine(db_engine)
                Base.metadata.create_all(engine)
                DBSession = sessionmaker(bind=engine)
                session = DBSession()
                print('Using {} for results...'.format(db_engine))
                session.close()

            except:
                print('Failed trying to make database connection '
                      'to {}'.format(args.create))

    #########################################################################
    #                       Database query                                  #
    #########################################################################
    # snowav -b 'basin' -wy wy -db mysql+mysqlconnector://mark:whatdystm?1@172.17.0.2/snowav

    # Database query for existing fields
    if args.bid and args.wy:

        # Making this a default for now
        value = 'swe_vol'

        # If no database is passed, try the default data location
        if not args.db:
            sbase = get_snowav_path()
            sbase = 'sqlite:///' + sbase
            args.db = sbase + '/snowav/data/snowav.db'

        start_date = datetime.datetime(args.wy-1,10,1)
        end_date = datetime.datetime(args.wy,9,30)

        try:
            engine = create_engine(args.db)
            connection = engine.connect()
            DBSession = sessionmaker(bind=engine)
            session = DBSession()

            # Get available run_names for the watershed
            wid =  Basins.basins[args.bid]['watershed_id']
            q = session.query(RunMetadata).join(Results).filter(and_(RunMetadata.watershed_id == wid),
                                                                    (Results.date_time >= start_date),
                                                                    (Results.date_time <= end_date),                                    )
            df = pd.read_sql(q.statement, q.session.connection())
            names = df.run_name.unique()

            # run_name and run_id
            rid = {}
            for name in names:
                rid[name] = df[df['run_name'] == name]['run_id'].values

            qry = session.query(Results).join(RunMetadata).filter(and_((Results.date_time >= start_date),
                                                      (Results.date_time <= end_date),
                                                      (RunMetadata.run_name.in_(names)),
                                                      (Results.basin_id == Basins.basins[args.bid]['basin_id'])))

            df = pd.read_sql(qry.statement, qry.session.connection())

            qry2 = session.query(RunMetadata).join(Results).filter(and_((Results.date_time >= start_date),
                                                      (Results.date_time <= end_date),
                                                      (RunMetadata.run_name.in_(names)) ))
            df2 = pd.read_sql(qry2.statement, qry2.session.connection())

            session.close()

            dash = '-'*100
            print('\n{}'.format(dash))
            print('On database: {}'.format(args.db))
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
            # print(name)
            # print(df2[df2['run_name'] == name]['data_location'].values)

        except:
            print('Failed database query in snow.py, format ' +
                  'snowav -b "Tuolumne" -wy 2018 -db mysql+mysqlconnector://user:password@172.17.0.2/snowav')

if __name__ == '__main__':
    run()
