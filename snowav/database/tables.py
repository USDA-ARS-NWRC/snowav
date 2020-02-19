
import sys
import os
import numpy as np
from sqlalchemy import Column, ForeignKey, Integer, String
from sqlalchemy.ext.declarative import declarative_base
from sqlalchemy.orm import relationship
from sqlalchemy import create_engine
from sqlalchemy import schema, types
from sqlalchemy.orm import sessionmaker
from sqlalchemy.orm import backref

Base = declarative_base()

class RunMetadata(Base):
    __tablename__ = 'RunMetadata'

    run_id = Column(Integer, primary_key=True)
    run_name = Column(String(250), nullable=False, index=True)
    watershed_id = Column(Integer, ForeignKey('Watershed.watershed_id'), nullable=True)
    pixel = Column(Integer, nullable=True)
    description = Column(String(250), nullable=True)
    smrf_version = Column(String(250), nullable=True)
    awsm_version = Column(String(250), nullable=True)
    snowav_version = Column(String(250), nullable=True)
    data_type = Column(String(250), nullable=True)
    data_location = Column(String(1200), nullable=True)
    file_type = Column(String(250), nullable=True)
    config_file = Column(String(250), nullable=True)
    proc_time = Column(types.DateTime(), nullable=True)

    # This puts Watershed fields available through Basin at Basin.watershed
    run_metadata = relationship('Watershed',
                               backref=backref('RunMetadata',lazy='dynamic'))

class Watershed(Base):
    __tablename__ = 'Watershed'

    watershed_id = Column(Integer, primary_key=True)
    watershed_name = Column(String(250), nullable=False, index=True)
    shapefile = Column(types.LargeBinary, nullable=True)


class Basin(Base):
    __tablename__ = 'Basin'

    watershed_id = Column(Integer,
                          ForeignKey('Watershed.watershed_id', onupdate="cascade"),
                          nullable=False)
    basin_id = Column(Integer, primary_key=True, nullable=False)
    basin_name = Column(String(250), nullable=False)
    shapefile = Column(types.LargeBinary, nullable=True)

    watershed = relationship('Watershed',
                             backref=backref('Basin',lazy='dynamic'))

class Results(Base):
    __tablename__ = 'Results'

    id = Column(Integer, primary_key=True)
    basin_id = Column(Integer, ForeignKey('Basin.basin_id', onupdate="cascade"), index=True)
    run_id = Column(Integer, ForeignKey('RunMetadata.run_id', onupdate="cascade"), index=True)
    date_time = Column(types.DateTime(), nullable=False, index=True)
    variable = Column(String(250), nullable=True, index=True)
    variable_id = Column(Integer, ForeignKey('VariableUnits.id', onupdate="cascade"), nullable=True)
    value = Column(types.Float(), nullable=True)
    elevation = Column(String(250), nullable=True)

    variable_units = relationship('VariableUnits',
                                  backref=backref('Results'))

    runmetadata = relationship('RunMetadata',
                               backref=backref('Results',lazy='dynamic'))

    basinid = relationship('Basin',
                               backref=backref('Results',lazy='dynamic'))

class VariableUnits(Base):
    __tablename__ = 'VariableUnits'

    id = Column(Integer, primary_key=True, autoincrement=True)
    run_id = Column(Integer, ForeignKey('RunMetadata.run_id', onupdate="cascade"),index=True)
    variable = Column(String(250), nullable=True)
    unit = Column(String(250), nullable=True)
    name = Column(String(250), nullable=True)

    run_metadata = relationship('RunMetadata',
                               backref=backref('VariableUnits',lazy='dynamic'))

class Inputs(Base):
    __tablename__ = 'Inputs'

    id = Column(Integer, primary_key=True)
    run_id = Column(Integer, nullable=False)
    run_name = Column(String(250), nullable=False, index=True)
    basin_id = Column(Integer, nullable=False, index=True)
    date_time = Column(types.DateTime(), nullable=False, index=True)
    variable = Column(String(250), nullable=False)
    function = Column(String(250), nullable=False)
    value = Column(types.Float(), nullable=True)
    unit = Column(String(250), nullable=True)
