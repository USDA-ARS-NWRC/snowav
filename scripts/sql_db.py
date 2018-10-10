
import numpy as np
from shutil import copyfile
import os
import copy
import pandas as pd
import datetime
from sqlalchemy import create_engine
from snowav.database.tables import Base, RunMetadata, Results, VariableUnits, Watersheds, Basins, Watershed, Basin
from sqlalchemy.ext.declarative import declarative_base
from sqlalchemy.orm import sessionmaker
import urllib.parse
import mysql.connector

user = 'mark'
password = 'whatdystm?1'
host = '127.0.0.1'
db = 'mysql_results'

cnx = mysql.connector.connect(user=user, password=password, host=host)
cursor = cnx.cursor()

# Check if database exists, create if necessary
query = ("SHOW DATABASES")
cursor.execute(query)
dbs = cursor.fetchall()
dbs = [i[0].decode("utf-8") for i in dbs]

if db not in dbs:
    query = ("CREATE DATABASE {};".format(db))
    cursor.execute(query)
    cursor.close()
    cnx.close()

    db_engine = 'mysql+mysqlconnector://{}:{}@{}/{}'.format(user,
                                                            password,
                                                            host,
                                                            db)
    engine = create_engine(db_engine)
    Base.metadata.create_all(engine)
    DBSession = sessionmaker(bind=engine)
    session = DBSession()

    # Initialize watersheds
    for ws in Watersheds.watersheds:
        wsid = Watersheds.watersheds[ws]['watershed_id']
        wsn = Watersheds.watersheds[ws]['watershed_name']
        wval = Watershed(watershed_id = wsid, watershed_name = wsn)
        session.add(wval)

        # Initialize basins within the watershed
        for bid in Basins.basins:
            if Basins.basins[bid]['watershed_id'] == wsid:
                bval = Basin(watershed_id = wsid,
                             basin_id = Basins.basins[bid]['basin_id'],
                             basin_name = Basins.basins[bid]['basin_name'])
                session.add(bval)

    session.commit()
    session.close()
