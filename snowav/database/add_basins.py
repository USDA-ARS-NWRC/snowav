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


'''
mysql -h 172.17.0.2 -P 3306 -u mark -p

'''
# ws = 'San Joaquin River Basin'
# watershedid = 3
basins = {
                'Auberry':{'watershed_id':3,
                                'basin_id':27,
                                'basin_name':'Auberry'},

                'Redinger':{'watershed_id':3,
                                'basin_id':28,
                                'basin_name':'Redinger'}
                                }

user = 'mark'
pwd = 'whatdystm?1'
host = '172.17.0.2'
port = '3306'
db = 'snowav'

db_engine = 'mysql+mysqlconnector://{}:{}@{}:{}/{}'.format(user,pwd,host,port,db)
engine = create_engine(db_engine)
Base.metadata.create_all(engine)
DBSession = sessionmaker(bind=engine)
session = DBSession()

wsid = Watersheds.watersheds[ws]['watershed_id']
# wsn = Watersheds.watersheds[ws]['watershed_name']
# wval = Watershed(watershed_id = wsid, watershed_name = wsn)

# Initialize basins within the watershed
for bid in basins:
    print('Adding {} to basins...'.format(bid))
    # if Basins.basins[bid]['watershed_id'] == wsid:
    bval = Basin(watershed_id = wsid,
                 basin_id = basins[bid]['basin_id'],
                 basin_name = bid)
    session.add(bval)


session.commit()
