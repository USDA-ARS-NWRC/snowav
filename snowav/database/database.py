
import os
import datetime
from sqlalchemy import Column, ForeignKey, Integer, String
from sqlalchemy.ext.declarative import declarative_base
from sqlalchemy import create_engine
from sqlalchemy.orm import sessionmaker
from sqlalchemy.orm import backref
from snowav.database.tables import Basin_Metadata, Base, Results, Run_Metadata
from sqlalchemy import and_

# Coming from scripts/snow.py
# values = {}
# snow.db_variables

def insert_results(snow,values):
    '''
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

    loc = snow.database

    engine = create_engine('sqlite:///%s'%(loc))
    DeclarativeBase = declarative_base()
    metadata = DeclarativeBase.metadata
    Base.metadata.bind = engine

    # Bind the engine to the metadata of the Base class so that the
    # declaratives can be accessed through a DBSession instance
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

def pull_results(loc, start_date, end_date, value):
    '''

    '''

    # loc = snow.database

    engine = create_engine('sqlite:///%s'%(loc))

    DBSession = sessionmaker(bind=engine)
    session = DBSession()

    qry = session.query(Results).filter(and_((Results.date_time > start_date),
                                          (Results.date_time < end_date),
                                          (Results.variable == value),
                                          (Results.basin_id == 1)))

    print(len(qry[:]),qry[0].value)
