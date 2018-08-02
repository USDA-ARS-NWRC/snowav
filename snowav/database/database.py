
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
    Inserts results to the database.

    Args
        loc: database location
        table: database table (currently Run_Metadata or Results)
        values: dictionary of values

    '''

    engine = create_engine('sqlite:///%s'%(loc))
    DeclarativeBase = declarative_base()
    metadata = DeclarativeBase.metadata
    Base.metadata.bind = engine
    DBSession = sessionmaker(bind=engine)
    session = DBSession()

    if table != 'Run_Metadata':
        val = Results(basin_id = values['basin_id'],
                      run_id = values['run_id'],
                      date_time = values['date_time'],
                      variable = values['variable'],
                      value = values['value'],
                      elevation = values['elevation'])

    else:
        val = Run_Metadata( run_id = values['run_id'],
                            basin_id = values['basin_id'],
                            basin_name = values['basin_name'],
                            description = values['description'],
                            smrf_version = values['smrf_version'],
                            awsm_version = values['awsm_version'],
                            snowav_version = values['snowav_version'],
                            data_type = values['data_type'],
                            data_location = values['data_location'],
                            file_type = values['file_type'],
                            config_file = values['config_file'],
                            proc_time = values['proc_time'],
                            var_units = values['var_units'],
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
    try:
        qry = session.query(Results).filter(and_((Results.date_time >= start_date),
                                          (Results.date_time <= end_date),
                                          (Results.basin_id == BASINS.basins[bid]['basin_id'])))

        # Only ever figured out how to do this one record at a time
        for act in qry:
            session.delete(act)

        # (Results.run_metadata.run_name == run_name)
        session.commit()
        session.close()
    except:
        print('delete exception')
        session.commit()
        session.close()

def run_metadata(self):
    '''
    Create run metadata entry for each snowav run.

    '''

    # Get latest run id
    engine = create_engine('sqlite:///%s'%(self.database))
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
              # 'run_name':self.run_name,
              'basin_name':self.plotorder[0],
              'description':'',
              'smrf_version':'smrf'+ smrf.__version__,
              'awsm_version':'aswm'+ awsm.__version__,
              'snowav_version':'snowav'+ snowav.__version__,
              'data_type':'',
              'data_location':'',
              'file_type':'',
              'config_file':self.config_file,
              'proc_time':datetime.datetime.now(),
              'var_units': (self.depthlbl + ', ' + self.vollbl + ', ' + 'MJ'
                            + ', ' + 'kg/m^3'),
              'elev_units':self.elevlbl
                }

    snowav.database.database.insert(self.database,'Run_Metadata',values)


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
                                          # (Results.run_metadata.run_name == self.run_name),
                                          (Results.basin_id == BASINS.basins[bid]['basin_id']),
                                          (Results.variable == value))).first()

    if qry is not None:
        flag = True
    else:
        flag = False

    session.close()

    return flag
