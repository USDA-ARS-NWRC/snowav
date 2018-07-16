
import os
from sqlalchemy.ext.declarative import declarative_base
from sqlalchemy import create_engine
from sqlalchemy.orm import sessionmaker
from snowav.database.tables import Basin_Metadata, Results, Run_Metadata

# loc = '/home/markrobertson/mrworkspace/projects/snowavdb/snowav_testy.db'

# Should eventually make a master list for these
values = {'basin_id':4,
          'basin_name':'Mores Creek',
          'state':'Idaho'}

def init_basin_metadata(loc,values):
    '''

    '''

    engine = create_engine('sqlite:///%s'%(loc))

    # Bind the engine to the metadata of the Base class so that the
    # declaratives can be accessed through a DBSession instance
    DBSession = sessionmaker(bind=engine)
    session = DBSession()

    val = Basin_Metadata(basin_id = values['basin_id'],
                         basin_name = values['basin_name'],
                         state = values['state']
                        )

    session.add(val)
    session.commit()
    session.close()

init_basin_metadata(loc,values)
