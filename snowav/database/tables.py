
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

Base = declarative_base()

class Basin_Metadata(Base):
    __tablename__ = 'Basin_Metadata'

    basin_id = Column(Integer, primary_key=True, autoincrement=True)
    basin_name = Column(String(250), nullable=False, unique=True)
    state = Column(String(250), nullable=True)
    area = Column(types.Float(), nullable=True)

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

class Run_Metadata(Base):
    __tablename__ = 'Run_Metadata'

    run_id = Column(Integer, primary_key=True)
    basin_id = Column(Integer, ForeignKey('Results.id'))
    basin_name = Column(String(250), nullable=True)
    description = Column(String(250), nullable=True)
    smrf_version = Column(String(250), nullable=True)
    awsm_version = Column(String(250), nullable=True)
    snowav_version = Column(String(250), nullable=True)
    data_type = Column(String(250), nullable=True)
    data_location = Column(String(250), nullable=True)
    file_type = Column(String(250), nullable=True)
    config_file = Column(String(250), nullable=True)

    run_metadata = relationship('Results',
                                backref=backref('run_metadata',lazy='dynamic'))
