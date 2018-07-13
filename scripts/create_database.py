
import sys
import os
from sqlalchemy import Column, ForeignKey, Integer, String
from sqlalchemy.ext.declarative import declarative_base
from sqlalchemy.orm import relationship
from sqlalchemy import create_engine
from sqlalchemy import schema, types
from sqlalchemy.orm import sessionmaker
from sqlalchemy.orm import backref

# db_location = '/home/markrobertson/mrworkspace/projects/snowav_database.db'

def run(db_location):

    if not os.path.isfile(db_location):

        Base = declarative_base()

        class Basin_Metadata(Base):
            __tablename__ = 'Basin_Metadata'

            basin_id = Column(Integer, primary_key=True, autoincrement=True)
            basin_name = Column(String(250), nullable=False, unique=True)
            state = Column(String(250), nullable=True)
            area = Column(types.Float(), nullable=False)

        class Results(Base):
            __tablename__ = 'Results'

            id = Column(Integer, primary_key=True)
            basin_id = Column(Integer, ForeignKey('Basin_Metadata.basin_id'))
            date_time = Column(types.DateTime(),nullable=False)
            proc_time = Column(types.DateTime(),nullable=True)
            version = Column(String(250), nullable=True)
            variable = Column(String(250), nullable=False)
            var_units = Column(String(250), nullable=False)
            value = Column(types.Float(), nullable=False)
            elevation = Column(String(250), nullable=False)
            elev_units = Column(String(250), nullable=False)

            # This puts Basin_Metadata.results and Results.basin_metadata
            basin_metadata = relationship('Basin_Metadata',
                                        backref=backref('results',lazy='dynamic'))

        # Create an engine that stores data in the local directory's
        # sqlalchemy_example.db file.
        engine = create_engine('sqlite:///%s'%(db_location))

        # Create all tables in the engine.
        Base.metadata.create_all(engine)

    else:

        print(
              ('Specified database at %s already exists, '%(db_location)
              + 'if you want to create an empty database in its place '
              + 'delete first and then re-run this script')
              )
        return

if __name__ == '__main__':
    db_location = sys.argv[1]
    run(db_location)
