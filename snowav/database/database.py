
import os
import datetime
from sqlalchemy import Column, ForeignKey, Integer, String
from sqlalchemy.ext.declarative import declarative_base
from sqlalchemy import create_engine
from sqlalchemy.orm import sessionmaker
from sqlalchemy.orm import backref
from snowav.database.tables import Base, RunMetadata, Watershed, Basin, Results, VariableUnits, Watersheds, Basins
from sqlalchemy import and_
import pandas as pd
import smrf
import awsm
import snowav
from snowav import database
import numpy as np
import copy
import urllib.parse
import mysql.connector


def insert(self, table, values):
    '''
    Inserts standard run results to the database. Uses database session created
    in read_config()

    Args
        table: database table (RunMetadata, Results, or VariableUnits)
        values: dictionary of values

    '''

    dbtable = getattr(snowav.database.tables, table)
    my_dbtable = dbtable()

    for k,v in values.items():
        setattr(my_dbtable, k, v)

    self.session.add(my_dbtable)
    self.session.commit()


def query(self, start_date, end_date, run_name, bid=None, value=None, rid=None):
    '''
    Query and retrieve results from the database, using database session created
    in read_config().

    Args
        start_date: (datetime)
        end_date: (datetime)
        run_name: identifier for run, specified in config file
        bid: basin id in string format ('Boise River Basin')
        value: value to query (i.e. 'swi_z')
        rid: run_id

    Returns
        df: dataframe of query results

    '''

    # Subset query with value and bid if desired
    if (value != None) and (bid != None):
        qry = self.session.query(Results).join(RunMetadata).filter(and_((Results.date_time >= start_date),
                                              (Results.date_time <= end_date),
                                              (RunMetadata.run_name == run_name),
                                              (Results.variable == value),
                                              (Results.basin_id == Basins.basins[bid]['basin_id'])))

    # Report queries here
    elif (value == None) and (bid != None) and rid != None:
        bids = [Basins.basins[n]['basin_id'] for n in bid]
        qry = self.session.query(Results).join(RunMetadata).filter(and_((Results.date_time >= start_date),
                                              (Results.date_time <= end_date),
                                              (Results.run_id == rid),
                                              (RunMetadata.run_name == run_name),
                                              (Results.basin_id.in_(bids))))

    else:
        qry = self.session.query(Results).join(RunMetadata).filter(and_((Results.date_time >= start_date),
                                              (Results.date_time <= end_date),
                                              (RunMetadata.run_name == run_name)))

    # x = self.session.query(Results).get(1)
    # print(x.runmetadata.run_name)
    df = pd.read_sql(qry.statement, qry.session.connection())

    return df

def delete(self, start_date, end_date, bid, run_name):
    '''
    Delete results from the database using database session created
    in read_config(). Currently this deletes all values with basin_id == bid
    and run_name = run_name within the date range, in other words, it removes
    all records of a previous run that share run_name with any date overlap.

    Args
        start_date: (datetime)
        end_date: (datetime)
        bid: basin id in string format ('Boise River Basin')
        run_name: identifier for run, specified in config file

    Table deletion order matters, by relationship.

    '''
    try:
        # Since we know run_name, but not run_id, get the run_id
        qry = self.session.query(Results).join(RunMetadata).\
                        filter(and_((Results.date_time >= start_date),
                        (Results.date_time <= end_date),
                        (RunMetadata.run_name == run_name),
                        (Results.basin_id == Basins.basins[bid]['basin_id'])))

        df = pd.read_sql(qry.statement, qry.session.connection())

        for r in df.run_id.unique():
            self.session.query(Results).\
                                filter(Results.run_id == int(r)).delete()
            self.session.query(VariableUnits).\
                                filter(VariableUnits.run_id == int(r)).delete()
            self.session.query(RunMetadata).\
                                filter(RunMetadata.run_id == int(r)).delete()
            self.session.flush()

        self.session.commit()

    except:
        print('Failed to delete database records in database.delete()...')

def create_tables(self, url=None):
    '''
    This function creates the Watersheds and Basins tables in the database
    from classes defined in /snowav/database/tables.py

    '''

    # Make database connection for duration of snowav processing
    if url is not None:
        engine = create_engine(url)
    else:
        engine = create_engine(self.database)

    Base.metadata.create_all(engine)
    DBSession = sessionmaker(bind=engine)
    self.session = DBSession()

    # Initialize watersheds
    for ws in Watersheds.watersheds:
        wsid = Watersheds.watersheds[ws]['watershed_id']
        wsn = Watersheds.watersheds[ws]['watershed_name']
        wval = Watershed(watershed_id = wsid, watershed_name = wsn)
        self.session.add(wval)

        # Initialize basins within the watershed
        for bid in Basins.basins:
            if Basins.basins[bid]['watershed_id'] == wsid:
                bval = Basin(watershed_id = wsid,
                             basin_id = Basins.basins[bid]['basin_id'],
                             basin_name = Basins.basins[bid]['basin_name'])
                self.session.add(bval)

    self.session.commit()
    print('Completed initializing basins and watersheds')

def run_metadata(self):
    '''
    Create run metadata entry for each snowav run, using database session
    created in read_config().

    '''

    qry = self.session.query(RunMetadata)
    df = pd.read_sql(qry.statement, qry.session.connection())

    # Increase each runid by 1
    if not df.run_id.values.size:
        self.runid = 1
    else:
        self.runid = int(df['run_id'].max() + 1)

    values = {
              'run_id':self.runid,
              'run_name':self.run_name,
              'watershed_id':Watersheds.watersheds[self.plotorder[0]]['watershed_id'],
              'pixel':self.pixel,
              'description':'',
              'smrf_version':'smrf'+ smrf.__version__,
              'awsm_version':'aswm'+ awsm.__version__,
              'snowav_version':'snowav'+ snowav.__version__,
              'data_type':'',
              'data_location':','.join(self.run_dirs),
              'file_type':'',
              'config_file':self.config_file,
              'proc_time':datetime.datetime.now()
                }

    snowav.database.database.insert(self, 'RunMetadata', values)

    # self.vars is defined in read_config(), and is also used in process()
    self.vid = {}
    for v in self.vars.keys():
        if (('vol' in v) or ('avail' in v)):
            u = self.units
        if (('z' in v) or (v == 'depth')):
            u = self.depthlbl
        if v == 'density':
            u = 'kg m^-3'
        if v == 'coldcont':
            u = 'MJ'

        variables = {
                     'run_id':self.runid,
                     'variable':v,
                     'unit':u,
                     'name':self.vars[v]
                    }

        snowav.database.database.insert(self,'VariableUnits',variables)

        # After inserting into VariableUnits, get the existing id
        x = self.session.query(VariableUnits).all()
        self.vid[v] = copy.deepcopy(x[-1].id)


def check_fields(self, start_date, end_date, bid, run_name, value):
    '''
    This functions queries the database and returns True if any value exists
    in the given date range, basin_id = id and run_name = run_name,
    and False if not.

    Args
        start_date: (datetime)
        end_date: (datetime)
        bid: basin id in string format ('Boise River Basin')
        run_name: identifier for run, specified in config file
        value: value to query (i.e. 'swi_z')

    Returns:
        flag: boolean (True if results exist)

    '''

    qry = self.session.query(Results).join(RunMetadata).filter(and_((Results.date_time >= start_date),
                                          (Results.date_time <= end_date),
                                          (RunMetadata.run_name == run_name),
                                          (Results.basin_id == Basins.basins[bid]['basin_id']),
                                          (Results.variable == value))).first()

    if qry is not None:
        flag = True
    else:
        flag = False

    return flag

def connect(self):
    '''
    This establishes a connection with a database for results, and will make
    a default sqlite database in the repo /snowav/data/model_results.db if
    no database is specified.

    '''

    # If no database is specified, make database snowav/data/model_results.db
    # Also, if specified database doesn't exist, it will be created

    # No database specified, create default sqlite in snowav/data
    if self.sqlite is None and self.mysql is None:

        # If this changes, update scripts/snow.py, which uses the default
        # location but does not always call read_config().
        dbpath = os.path.abspath((self.snowav_path + '/snowav/data/'))
        database = os.path.abspath(dbpath + '/model_results.db')
        self.database = 'sqlite:///{}'.format(database)

        # Create metadata tables and keep open self.session
        if not os.path.isfile(database):
            print('No results database specified, creating default sqlite database '
                  '{}'.format(self.database))

            create_tables(self)

        else:
            engine = create_engine(self.database)
            DBSession = sessionmaker(bind=engine)
            self.session = DBSession()
            print('Using {} for results...'.format(self.database))

    # sqlite specified, but doesn't exist
    if self.sqlite is not None:
        fp = urllib.parse.urlparse(self.sqlite)

        if (not os.path.isfile(fp.path)):

            dbpath = os.path.abspath((self.snowav_path + '/snowav/data/'))
            database = os.path.abspath(dbpath + '/model_results.db')
            self.database = 'sqlite:///{}'.format(database)

            print('Database specified in config file does not exist, creating '
                  'default {}'.format(self.database))

            # Create metadata tables and keep open self.session
            create_tables(self)
            engine = create_engine(self.database)

        else:
            engine = create_engine(self.sqlite)

        DBSession = sessionmaker(bind=engine)
        self.session = DBSession()
        print('Using {} for results...'.format(self.sqlite))

    if self.mysql is not None:

        # First, check if database exists, if not, make it
        cnx = mysql.connector.connect(user=self.db_user,
                                      password=self.db_password,
                                      host=self.db_host)
        cursor = cnx.cursor()

        # Check if database exists, create if necessary
        query = ("SHOW DATABASES")
        cursor.execute(query)
        dbs = cursor.fetchall()
        dbs = [i[0].decode("utf-8") for i in dbs]

        db_engine = 'mysql+mysqlconnector://{}:{}@{}/{}'.format(self.db_user,
                                                                self.db_password,
                                                                self.db_host,
                                                                self.mysql)

        # If the database doesn't exist, create it, otherwise connect
        if self.mysql not in dbs:
            print('Specified mysql database {} does not exist, it is being '
                  'created...'.format(self.mysql))
            query = ("CREATE DATABASE {};".format(self.mysql))
            cursor.execute(query)
            cursor.close()
            cnx.close()

            create_tables(self, url = db_engine)

        else:
            try:
                engine = create_engine(db_engine)
                Base.metadata.create_all(engine)
                DBSession = sessionmaker(bind=engine)
                self.session = DBSession()
                print('Using {} for results...'.format(db_engine))

            except:
                print('Failed trying to make database connection '
                      'to {}'.format(self.mysql))

def write_csv(self):
    '''
    This writes out variables, specified in the config file, from the database
    to the figs_path directory

    '''

    # Initialize
    out = pd.DataFrame(columns = self.plotorder)

    # By variable
    for var in self.write_csv:

        # For all subbasins
        for bid in self.plotorder:
            r = database.database.query(self, datetime.datetime(self.wy-1,10,1),
                                        self.end_date, self.run_name, bid, var)

            # For now, just write out totals
            v = r[(r['elevation'] == 'total')]

            for iter,d in enumerate(v['date_time'].values):
                out.loc[d,bid] = v['value'].values[iter]

            # SWI is cumulative
            if var == 'swi_vol':
                out = out.cumsum()

        if self.vollbl == 'KAF':
            out.to_csv(os.path.join(self.figs_path,var + '_' + self.vollbl + '.csv'))

        if self.vollbl == 'SI':
            out.to_csv(os.path.join(self.figs_path,var + '_' + 'm3' + '.csv'))
