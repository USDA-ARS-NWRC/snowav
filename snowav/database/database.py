
import os
import datetime
from sqlalchemy import Column, ForeignKey, Integer, String
from sqlalchemy.ext.declarative import declarative_base
from sqlalchemy import create_engine
from sqlalchemy.orm import sessionmaker
from sqlalchemy.orm import backref
from snowav.database.tables import Basin_Metadata, Base, Results, Run_Metadata, BASINS
from sqlalchemy import and_
import pandas as pd
import smrf
import awsm
import snowav
import numpy as np

def insert(loc,table,values):
    '''
    Inserts standard run results to the database.

    Args
        loc: database location
        table: database table (currently Run_Metadata or Results)
        values: dictionary of values

    '''

    engine = create_engine(loc)
    DeclarativeBase = declarative_base()
    metadata = DeclarativeBase.metadata
    Base.metadata.bind = engine
    DBSession = sessionmaker(bind=engine)
    session = DBSession()

    dbtable = getattr(snowav.database.tables, table)
    my_dbtable = dbtable()

    for k,v in values.items():
        setattr(my_dbtable, k, v)

    session.add(my_dbtable)
    session.commit()
    session.close()


def query(loc, start_date, end_date, run_name, bid = None, value = None):
    '''
    Query and retrieve results from the database.

    Args
        loc: database location
        start_date: (datetime)
        end_date: (datetime)
        run_name: identifier for run, specified in config file
        bid: basin id in string format ('Boise River Basin')
        value: value to query (i.e. 'swi_z')

    Returns:
        dataframe of query results

    '''

    engine = create_engine(loc)
    connection = engine.connect()
    DBSession = sessionmaker(bind=engine)
    session = DBSession()

    # Subset query with value and bid if desired
    if (value != None) and (bid != None):
        qry = session.query(Results).filter(and_((Results.date_time >= start_date),
                                              (Results.date_time <= end_date),
                                              (Results.run_name == run_name),
                                              (Results.variable == value),
                                              (Results.basin_id == BASINS.basins[bid]['basin_id'])))

    # Otherwise you're gonna get a bunch of stuff...
    else:
        qry = session.query(Results).filter(and_((Results.date_time >= start_date),
                                              (Results.date_time <= end_date),
                                              (Results.run_name == run_name)))

    df = pd.read_sql(qry.statement, qry.session.connection())
    session.close()

    return df

def delete(loc, start_date, end_date, bid, run_name):
    '''
    Delete results from the database. Currently this deletes all values
    with basin_id == bid and run_name = run_name within the date range.

    Args
        loc: database location
        start_date: (datetime)
        end_date: (datetime)
        bid: basin id in string format ('Boise River Basin')
        run_name: identifier for run, specified in config file

    '''

    engine = create_engine(loc)
    connection = engine.connect()
    DBSession = sessionmaker(bind=engine)
    session = DBSession()

    try:
        qry = session.query(Results).filter(and_((Results.date_time >= start_date),
                                          (Results.date_time <= end_date),
                                          (Results.run_name == run_name),
                                          (Results.basin_id == BASINS.basins[bid]['basin_id'])))

        # Only ever figured out how to do this one record at a time
        for f in qry:
            session.delete(f)

        session.commit()
        session.close()

    except:
        print('Failed during attempted deletion in database.delete()')
        session.commit()
        session.close()

def run_metadata(self):
    '''
    Create run metadata entry for each snowav run.

    '''

    # Get latest run id
    engine = create_engine(self.database)
    connection = engine.connect()
    DBSession = sessionmaker(bind=engine)
    session = DBSession()

    qry = session.query(Run_Metadata)
    df = pd.read_sql(qry.statement, qry.session.connection())
    session.close()

    # Increase each runid by 1
    if not df.run_id.values.size:
        self.runid = 1
    else:
        self.runid = int(df['run_id'].max() + 1)

    values = {
              'run_id':self.runid,
              'basin_id':BASINS.basins[self.plotorder[0]]['basin_id'],
              'run_name':self.run_name,
              'basin_name':self.plotorder[0],
              'description':'',
              'smrf_version':'smrf'+ smrf.__version__,
              'awsm_version':'aswm'+ awsm.__version__,
              'snowav_version':'snowav'+ snowav.__version__,
              'data_type':'',
              'data_location':','.join(self.run_dirs),
              'file_type':'',
              'config_file':self.config_file,
              'proc_time':datetime.datetime.now(),
              'var_units': (self.depthlbl + ', ' + self.vollbl + ', ' + 'MJ'
                            + ', ' + 'kg/m^3'),
              'elev_units':self.elevlbl
                }

    snowav.database.database.insert(self.database,'Run_Metadata',values)


def check_fields(loc, start_date, end_date, bid, run_name, value):
    '''
    This functions queries the database and returns True if any value exists
    in the given date range, basin_id = id and run_name = run_name,
    and False if not.

    Args
        loc: database location
        start_date: (datetime)
        end_date: (datetime)
        bid: basin id in string format ('Boise River Basin')
        run_name: identifier for run, specified in config file
        value: value to query (i.e. 'swi_z')

    Returns:
        flag: boolean (True if results exist)

    '''

    engine = create_engine(loc)
    connection = engine.connect()
    DBSession = sessionmaker(bind=engine)
    session = DBSession()

    qry = session.query(Results).filter(and_((Results.date_time >= start_date),
                                          (Results.date_time <= end_date),
                                          (Results.run_name == run_name),
                                          (Results.basin_id == BASINS.basins[bid]['basin_id']),
                                          (Results.variable == value))).first()

    if qry is not None:
        flag = True
    else:
        flag = False

    session.close()

    return flag
