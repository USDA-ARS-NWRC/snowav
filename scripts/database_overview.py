
import sys
import os
import datetime
from sqlalchemy import Column, ForeignKey, Integer, String
from sqlalchemy.ext.declarative import declarative_base
from sqlalchemy import create_engine
from sqlalchemy.orm import sessionmaker
from sqlalchemy.orm import backref
from snowav.database.tables import BasinMetadata, Base, Results, RunMetadata, BASINS
from sqlalchemy import and_
from snowav import database
from snowav.database.tables import BASINS
import pandas as pd
import smrf
import awsm
import snowav
import numpy as np
from matplotlib import pyplot as plt
import seaborn as sns
import copy
import matplotlib.dates as mdates
import argparse

'''
List available runs on database by basin and water year.

$ python scripts/database_overview.py -b Basin River Basin -wy 2017

'''

if __name__ == '__main__':

    parser = argparse.ArgumentParser(description="Examine and auto populate")

    parser.add_argument('-b', '--bid', dest='bid', type=str,
                        help='Basin id, from snowav/database/tables.py')

    parser.add_argument('-wy', '--wy', dest='wy', type=int,
                        help='Water year')

    args = parser.parse_args()

# Making these defaults for now
location = 'sqlite:////home/markrobertson/wkspace/projects/snowavdb/snowavdb.db'
value = 'swe_vol'

start_date = datetime.datetime(args.wy - 1,10,1)
end_date = datetime.datetime(args.wy,9,30)

# Connect to database and query
engine = create_engine(location)
connection = engine.connect()
DBSession = sessionmaker(bind=engine)
session = DBSession()
qry = session.query(Results).filter(and_((Results.date_time >= start_date),
                                          (Results.date_time <= end_date),
                                          (Results.basin_id == BASINS.basins[args.bid]['basin_id'])))

df = pd.read_sql(qry.statement, qry.session.connection())

session.close()

print('On database: {},\nAvailable runs for {}: {},\nDate range: {}, {}'.format(
location,args.bid,df.run_name.unique(),df['date_time'].min(),df['date_time'].max()))
