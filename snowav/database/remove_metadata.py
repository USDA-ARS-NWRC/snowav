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
try:
    import smrf
except:
    print('Could not import smrf')
try:
    import awsm
except:
    print('Could not import awsm')
import snowav
from snowav import database
import numpy as np
import copy
import urllib.parse
import mysql.connector


'''
mysql -h 172.17.0.2 -P 3306 -u mark -p

snowav -b 'basin' -wy wy -db mysql+mysqlconnector://<user>:<pwd>@172.17.0.2/snowav

'''

run_name = 'sj_wy2019_forecast'
start_date = datetime.datetime(2018,10,1)
end_date = datetime.datetime(2019,9,30)
# basins = ['Middle Fork','West Kings','Dinkey Creek','Middle South Fork','South Fork','Mill Creek','North Fork']
basins = ['San Joaquin River Basin','South Fork','Main','Auberry','Redinger']

user = 'user'
pwd = 'password'
host = '172.17.0.2'
port = '3306'
db = 'snowav'

db_engine = 'mysql+mysqlconnector://{}:{}@{}:{}/{}'.format(user,pwd,host,port,db)
engine = create_engine(db_engine)
Base.metadata.create_all(engine)
DBSession = sessionmaker(bind=engine)
session = DBSession()

for bid in basins:
    wid =  Basins.basins[bid]['watershed_id']

    qry = session.query(Results).join(RunMetadata).filter(and_(
                    (Results.date_time >= start_date),
                    (Results.date_time <= end_date),
                    (RunMetadata.run_name == run_name),
                    (Results.basin_id == Basins.basins[bid]['basin_id'])))

    df = pd.read_sql(qry.statement, qry.session.connection())

    # Delete Results by date range, run_id, and basin
    for r in df.run_id.unique():
        session.query(Results).filter(and_(
            (Results.date_time >= start_date),
            (Results.date_time <= end_date),
            (Results.run_id == int(r)),
            (Results.basin_id == Basins.basins[bid]['basin_id']))).delete()

        session.flush()

    # Query Metadata to see what's remaining
    qry = session.query(RunMetadata).join(Results).filter(and_(RunMetadata.watershed_id == wid),
                                                        (RunMetadata.run_name == run_name),
                                                        (Results.date_time >= start_date),
                                                        (Results.date_time <= end_date))


    df = pd.read_sql(qry.statement, qry.session.connection())

    # run_name and run_id
    rid = {}
    rid[run_name] = df[df['run_name'] == run_name]['run_id'].values

    # Toggling

    # for n in rid[run_name]:
    #     session.query(RunMetadata).filter(RunMetadata.run_id == int(n)).delete()
    #     session.flush()
    #     session.commit()


    # for n in rid[run_name]:
    #         session.query(VariableUnits).filter(VariableUnits.run_id == int(n)).delete()
    #         session.flush()
    #         session.commit()

session.close()
