
import sys
import os
import datetime
from sqlalchemy import Column, ForeignKey, Integer, String
from sqlalchemy.ext.declarative import declarative_base
from sqlalchemy import create_engine
from sqlalchemy.orm import sessionmaker
from sqlalchemy.orm import backref
from snowav.database.tables import Basin_Metadata, Base, Results, Run_Metadata, BASINS
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

'''
List available runs on database by basin and water year.

python scripts/database_overview.py 'Basin' wy

'''

# Making these defaults for now
location = 'sqlite:////home/markrobertson/wkspace/projects/snowavdb/snowavdb.db'
value = 'swe_vol'

start_date = datetime.datetime(wy - 1,10,1)
end_date = datetime.datetime(wy,9,30)

# Connect to database and query
engine = create_engine(location)
connection = engine.connect()
DBSession = sessionmaker(bind=engine)
session = DBSession()
qry = session.query(Results).filter(and_((Results.date_time >= start_date),
                                          (Results.date_time <= end_date),
                                          (Results.basin_id == BASINS.basins[bid]['basin_id'])))

df = pd.read_sql(qry.statement, qry.session.connection())
session.close()

print('On database: {},\nAvailable runs for {}: {},\nDate range: {}, {}'.format(
location,bid,df.run_name.unique(),df['date_time'].min(),df['date_time'].max()))

if __name__ == '__main__':
    bid = sys.argv[1]

    if len(sys.argv) > 2:
        wy = sys.argv[2]
    else:
        wy = 2018
