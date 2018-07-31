
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

def insert(loc,values):
    '''
    Inserts results to the database.

    Args
        values in the following format:
            values = {'basin_id': 1,
                      'date_time': datetime.datetime.now(),
                      'proc_time': datetime.datetime.now(),
                      'version': '1',
                      'variable': 'swe',
                      'var_units': 'in',
                      'value': 2.0,
                      'elevation': 'total',
                      'elev_units': 'ft'}

    '''

    engine = create_engine('sqlite:///%s'%(loc))
    DeclarativeBase = declarative_base()
    metadata = DeclarativeBase.metadata
    Base.metadata.bind = engine
    DBSession = sessionmaker(bind=engine)
    session = DBSession()

    val = Results(basin_id = values['basin_id'],
                  date_time = values['date_time'],
                  proc_time = values['proc_time'],
                  version = values['version'],
                  variable = values['variable'],
                  var_units = values['var_units'],
                  value = values['value'],
                  elevation = values['elevation'],
                  elev_units = values['elev_units'] )

    session.add(val)
    session.commit()
    session.close()

def query(loc, start_date, end_date, bid = None, value = None):
    '''
    Query and retrieve results from the database.

    Args
        loc: database location
        start_date: (datetime)
        end_date: (datetime)
        bid: basin id in string format ('Boise River Basin')
        value: value to query (i.e. 'swi_z')

    Returns:
        dataframe of query results

    '''

    engine = create_engine('sqlite:///%s'%(loc))
    connection = engine.connect()
    DBSession = sessionmaker(bind=engine)
    session = DBSession()

    # Subset query with value and bid if desired
    if (value != None) and (bid != None):
        qry = session.query(Results).filter(and_((Results.date_time >= start_date),
                                              (Results.date_time <= end_date),
                                              (Results.variable == value),
                                              (Results.basin_id == BASINS.basins[bid]['basin_id'])))
    else:
        qry = session.query(Results).filter(and_((Results.date_time >= start_date),
                                              (Results.date_time <= end_date)))

    df = pd.read_sql(qry.statement, qry.session.connection())
    session.close()

    return df

def delete(loc, start_date, end_date, bid):
    '''
    Delete results from the database. Currently this deletes all values with
    basin_id == bid in the selected date range.

    Args
        loc: database location
        start_date: (datetime)
        end_date: (datetime)
        bid: basin id in string format ('Boise River Basin')

    '''

    engine = create_engine('sqlite:///%s'%(loc))
    connection = engine.connect()
    DBSession = sessionmaker(bind=engine)
    session = DBSession()

    qry = session.query(Results).filter(and_((Results.date_time >= start_date),
                                          (Results.date_time <= end_date),
                                          (Results.basin_id == BASINS.basins[bid]['basin_id'])))

    # print(bid,qry)
    qry.delete()
    session.commit()
    session.close()


def check_fields(loc, start_date, end_date, bid, value):
    '''
    This functions queries the database and returns True if any value exists in
    the given date range, and False if not.

    Args
        loc: database location
        start_date: (datetime)
        end_date: (datetime)
        bid: basin id in string format ('Boise River Basin')
        value: value to query (i.e. 'swi_z')

    Returns:
        flag: boolean (True if results exist)    

    '''

    engine = create_engine('sqlite:///%s'%(loc))
    connection = engine.connect()
    DBSession = sessionmaker(bind=engine)
    session = DBSession()

    qry = session.query(Results).filter(and_((Results.date_time >= start_date),
                                          (Results.date_time <= end_date),
                                          (Results.basin_id == BASINS.basins[bid]['basin_id']),
                                          (Results.variable == value))).first()

    if qry is not None:
        flag = True
    else:
        flag = False

    session.close()

    return flag
