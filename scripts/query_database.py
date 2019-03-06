

import numpy as np
import snowav
import datetime
import pandas as pd
from sqlalchemy.ext.declarative import declarative_base
from sqlalchemy import create_engine
from sqlalchemy.orm import sessionmaker
from sqlalchemy import and_
import urllib.parse
import mysql.connector
from snowav.database.tables import (Base, RunMetadata, Watershed, Basin,
Results, VariableUnits, Watersheds, Basins)
from snowav import database
from sys import exit

# snowav -b 'basin' -wy wy -db mysql+mysqlconnector://mark:whatdystm?1@172.17.0.2/snowav
wy = 2019
basins = ['San Joaquin River Basin']
value = 'swe_vol'
run_name = 'sj_wy2019_forecast'
start_date = datetime.datetime(wy-1,10,1)
end_date = datetime.datetime(wy,9,30)
total = False
# print, csv
output = 'prit'

print(vars(Basins))
##############################################################################
db = 'mysql+mysqlconnector://mark:whatdystm?1@172.17.0.2/snowav'
values = ['swe_vol','swe_avail','swe_unavail','swe_z',
          'swi_vol','swi_z','precip_vol','precip_z',
          'depth','density','rain_z','evap_z','coldcont']

if value not in values:
    print('value option must be one of: {}'.format(values))
    exit()

if output not in ['print','csv']:
    print('output option must be "print" or "csv"')
    exit()

engine = create_engine(db)
connection = engine.connect()
DBSession = sessionmaker(bind=engine)
session = DBSession()

# Make df from database
db_val = pd.DataFrame(columns = basins)

for bid in basins:
    qry = session.query(Results).join(RunMetadata).filter(and_(
                    (Results.date_time >= start_date),
                    (Results.date_time <= end_date),
                    (RunMetadata.run_name == run_name),
                    (Results.variable == value),
                    (Results.basin_id == Basins.basins[bid]['basin_id'])))

    df = pd.read_sql(qry.statement, qry.session.connection())

    if total is True:
        fqry = df[df['elevation'] == 'total']
    else:
        fqry = df

    db_val.loc['total',bid] = fqry['value'].values

# print(db_val)
cols = ['date_time','value','elevation']

if output == 'print':
    with pd.option_context('display.max_rows', None, 'display.max_columns', None):
        print('\n\nsnowav database values \n'
              'run_name: {}\n'
              '{} to {}\n'
              'value: {}\n\n'
              '{}'.format(run_name, start_date, end_date, value, fqry[cols]))

if output == 'csv':
    db_val.to_csv
