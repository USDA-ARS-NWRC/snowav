
import os
import sys
from sqlalchemy import create_engine
from snowav.database.tables import Base, Basin_Metadata, Results, Run_Metadata

'''
This script creates a sqlite database with tables defined in
snowav/database/tables.py

After creation of the database basin metadata needs to be created, to do so
add all necessary definitions in snowav/database/tables.py and run
scripts/create_metadata.py

python scripts/create_database.py <path to database file>

To add subbasin metadata after database already exists:
python scripts/create_metadata.py

'''

def run(db_location):

    if not os.path.isfile(db_location):

        # Create an engine that stores data in the local directory's
        # sqlalchemy_example.db file.
        dbstr = 'sqlite:///%s'%(db_location)
        engine = create_engine(dbstr)

        # Create all tables in the engine.
        Base.metadata.create_all(engine)

    else:
        print('Database already exists. If you really want to replace it, '
        + 'delete manually. I am definitely not doing that for you.')

if __name__ == '__main__':
    db_location = sys.argv[1]
    run(db_location)
