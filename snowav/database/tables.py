
import sys
import os
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
    # run_id = Column(Integer, ForeignKey('RunMetadata.run_id'))
    # basins = Column(String(250), nullable=True)
    shapefile = Column(types.LargeBinary, nullable=True)


class Basin(Base):
    __tablename__ = 'Basin'

    watershed_id = Column(Integer,
                          ForeignKey('Watershed.watershed_id'),
                          nullable=False)
    basin_id = Column(Integer, primary_key=True, nullable=False)
    basin_name = Column(String(250), nullable=False)
    shapefile = Column(types.LargeBinary, nullable=True)

    watershed = relationship('Watershed',
                             backref=backref('Basin',lazy='dynamic'))

class Results(Base):
    __tablename__ = 'Results'

    id = Column(Integer, primary_key=True)
    basin_id = Column(Integer, ForeignKey('Basin.basin_id'), index=True)
    run_id = Column(Integer, ForeignKey('RunMetadata.run_id'), index=True)
    date_time = Column(types.DateTime(), nullable=False, index=True)
    variable = Column(String(250), nullable=True, index=True)
    variable_id = Column(Integer, ForeignKey('VariableUnits.id'), nullable=True)
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
    run_id = Column(Integer, ForeignKey('RunMetadata.run_id'),index=True)
    variable = Column(String(250), nullable=True)
    unit = Column(String(250), nullable=True)
    name = Column(String(250), nullable=True)

    run_metadata = relationship('RunMetadata',
                               backref=backref('VariableUnits',lazy='dynamic'))

class Watersheds(object):
    '''
    To add a watershed:
        'Basin Name':{'watershed_id':n+1,
                             'watershed_name':'Basin Name short version',
                             'basins': '',
                             'shapefile':'empty'}

    Also fill out Basins with relevant subbasins.

    Basins main key and watershed_id must match a Watershed!

    '''

    watersheds = {
                'Boise River Basin':{'watershed_id':1,
                                     'watershed_name':'Boise',
                                     'basins': '',
                                     'shapefile':'empty'},

                'Extended Tuolumne':{'watershed_id':2,
                                     'watershed_name':'Extended Tuolumne',
                                     'basins': '',
                                     'shapefile':'empty'},

                'San Joaquin River Basin':{'watershed_id':3,
                                           'watershed_name':'San Joaquin',
                                           'basins': '',
                                           'shapefile':'empty'},

                'Lakes Basin':{'watershed_id':4,
                               'watershed_name':'Lakes',
                               'basins': '',
                               'shapefile':'empty'},

                'Merced River Basin':{'watershed_id':5,
                                      'watershed_name':'Merced',
                                      'basins': '',
                                      'shapefile':'empty'},

                'Kaweah River Basin':{'watershed_id':6,
                                      'watershed_name':'Kaweah',
                                      'basins': '',
                                      'shapefile':'empty'},

                'Kings River Basin':{'watershed_id':7,
                                     'watershed_name':'Kings',
                                     'basins': '',
                                     'shapefile':'empty'},

                'Reynolds Creek':{'watershed_id':8,
                                  'watershed_name':'Reynolds Creek',
                                  'basins': '',
                                  'shapefile':'empty'}
                }


class Basins(object):
    '''
    This class defines the basin metadata for all of the subbasins that are
    used in the database.

    To add basins, Basins main key and watershed_id must match the
    correct Watershed!

    '''
    basins = {
                # BRB
                'Boise River Basin':{'watershed_id':1,
                                     'basin_id':1,
                                     'basin_name':'Boise River Basin'},
                'Featherville':{'watershed_id':1,
                                'basin_id':2,
                                'basin_name':'Featherville'},
                'Twin Springs':{'watershed_id':1,
                                'basin_id':3,
                                'basin_name':'Twin Springs'},
                'Mores Creek':{'watershed_id':1,
                               'basin_id':4,
                               'basin_name':'Mores Creek'},

                # Tuolumne
                'Extended Tuolumne':{'watershed_id':2,
                                     'basin_id':5,
                                     'basin_name':'Extended Tuolumne'},
                'Tuolumne':{'watershed_id':2,
                            'basin_id':6,
                            'basin_name':'Tuolumne'},
                'Cherry Creek':{'watershed_id':2,
                                'basin_id':7,
                                'basin_name':'Cherry Creek'},
                'Eleanor':{'watershed_id':2,
                           'basin_id':8,
                           'basin_name':'Eleanor'},

                # San Joaquin
                'San Joaquin River Basin':{'watershed_id':3,
                                           'basin_id':9,
                                           'basin_name':'San Joaquin'},

                'South Fork':{'watershed_id':3,
                              'basin_id':10,
                              'basin_name':'South Fork'},

                'Main':{'watershed_id':3,
                        'basin_id':11,
                        'basin_name':'Main'},

                'Jose Creek':{'watershed_id':3,
                              'basin_id':12,
                              'basin_name':'Jose Creek'},
                              
                'Willow Creek':{'watershed_id':3,
                                'basin_id':13,
                                'basin_name':'Willow Creek'},

                # Reynolds
                'Reynolds Creek':{'watershed_id':8,
                                  'basin_id':14,
                                  'basin_name':'Reynolds Creek'},
                'Tollgate':{'watershed_id':8,
                            'basin_id':15,
                            'basin_name':'Tollgate'},

                # Lakes
                'Lakes Basin':{'watershed_id':4,
                               'basin_id':16,
                               'basin_name':'Lakes'},

                # Merced
                'Merced River Basin':{'watershed_id':5,
                                      'basin_id':17,
                                      'basin_name':'Merced'},

                # Kaweah
                'Kaweah River Basin':{'watershed_id':6,
                                      'basin_id':18,
                                      'basin_name':'Kaweah'},

                # Kings
                'Kings River Basin':{'watershed_id':7,
                                     'basin_id':19,
                                     'basin_name':'Kings'},

                'Middle Fork':{'watershed_id':7,
                              'basin_id':20,
                              'basin_name':'Kings'},

                'West Kings':{'watershed_id':7,
                              'basin_id':21,
                              'basin_name':'Kings'},

                'Middle South Fork':{'watershed_id':7,
                              'basin_id':22,
                              'basin_name':'Kings'},

                'South Fork':{'watershed_id':7,
                              'basin_id':23,
                              'basin_name':'Kings'},

                'Mill Creek':{'watershed_id':7,
                              'basin_id':24,
                              'basin_name':'Kings'},

                'North Fork':{'watershed_id':7,
                              'basin_id':25,
                              'basin_name':'Kings'},

                'Dinkey Creek':{'watershed_id':7,
                              'basin_id':26,
                              'basin_name':'Kings'}


                }
