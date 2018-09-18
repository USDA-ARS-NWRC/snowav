
import os
from sqlalchemy.ext.declarative import declarative_base
from sqlalchemy import create_engine
from sqlalchemy.orm import sessionmaker
from snowav.database.tables import BasinMetadata, Results, RunMetadata, BASINS

'''
This scipt creates the basin and subbasin fields in Basin_Metadata. Those fields
are defined in snowav/database/tables.

To add subbasins:
python scripts/create_metadata.py

'''

loc = '/home/markrobertson/wkspace/projects/snowavdb/snowavdb.db'

def init_basin_metadata(loc, values):

    engine = create_engine('sqlite:///%s'%(loc))
    DBSession = sessionmaker(bind=engine)
    session = DBSession()

    val = Basin_Metadata(basin_id = values['basin_id'],
                         basin_name = values['basin_name'],
                         state = values['state']
                        )

    session.add(val)
    session.commit()
    session.close()

# Initialize basin metadata for all basin defitions
for basin in BASINS.basins:
    # Cleaner way would be to query first
    try:
        init_basin_metadata(loc, BASINS.basins[basin])
    except:
        print('Addition of basin metadata failed on %s'%(basin)
              + ', this may be because it already exists.' )
