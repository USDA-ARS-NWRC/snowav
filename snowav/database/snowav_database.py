import numpy as np
import os
import sys
import datetime
from sqlalchemy import Column, ForeignKey, Integer, String
from sqlalchemy.ext.declarative import declarative_base
from sqlalchemy.orm import relationship
from sqlalchemy import create_engine
from sqlalchemy import schema, types
from sqlalchemy.orm import sessionmaker
from sqlalchemy.orm import backref

def update():
    # After db creation...
    DeclarativeBase = declarative_base()
    metadata = DeclarativeBase.metadata

    # Bind the engine to the metadata of the Base class so that the
    # declaratives can be accessed through a DBSession instance
    Base.metadata.bind = engine
    DBSession = sessionmaker(bind=engine)
    session = DBSession()

    basin_names = ['Boise River Basin','Featherville']
    areas = np.array((100,30))

    for bn,a in zip(basin_names,areas):
        brb = Basin_Metadata(basin_name = bn, state = 'Idaho', area = a)
        session.add(brb)

    session.commit()
    session.close()
