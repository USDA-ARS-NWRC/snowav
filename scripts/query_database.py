

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
from snowav.database.tables import Base, RunMetadata, Watershed, Basin, Results, VariableUnits, Watersheds, Basins
from snowav import database

# snowav -b 'basin' -wy wy -db mysql+mysqlconnector://mark:whatdystm?1@172.17.0.2/snowav
db = 'mysql+mysqlconnector://mark:whatdystm?1@172.17.0.2/snowav'
wy = 2019
basins = ['San Joaquin River Basin']
value = 'evap_z'
run_name = 'sj_wy2019_ops'

start_date = datetime.datetime(wy-1,10,1)
end_date = datetime.datetime(wy,9,30)
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
    fqry = df[df['elevation'] == 'total']
    db_val.loc['total',bid] = np.nansum(fqry['value'].values)

print(db_val)
