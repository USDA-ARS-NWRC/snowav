
import sys
import os
import snowav
from sqlalchemy import Column, ForeignKey, Integer, String
from sqlalchemy.ext.declarative import declarative_base
from sqlalchemy.orm import relationship
from sqlalchemy import create_engine
from sqlalchemy import schema, types
from sqlalchemy.orm import sessionmaker
from sqlalchemy.orm import backref
from snowav.database.tables import Base, Basin_Metadata, Results

'''
Currently this isn't working with relative 'db_location'...

'''


def run(db_location):

    if not os.path.isfile(db_location):

        # Base = declarative_base()
        # Basin_Metadata(Base)
        # Results(Base)

        # Create an engine that stores data in the local directory's
        # sqlalchemy_example.db file.
        engine = create_engine('sqlite:////home/markrobertson/mrworkspace/projects/snowavdb/snowav_db_test.db')

        # Create all tables in the engine.
        Base.metadata.create_all(engine)

        # After db creation...
        # DeclarativeBase = declarative_base()
        # metadata = DeclarativeBase.metadata

        # Bind the engine to the metadata of the Base class so that the
        # declaratives can be accessed through a DBSession instance
        # Base.metadata.bind = engine
        # DBSession = sessionmaker(bind=engine)
        # session = DBSession()
        # session.commit()
        # session.close()


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
