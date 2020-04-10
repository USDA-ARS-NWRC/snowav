
from sqlalchemy import types, Column, ForeignKey, Integer, String, Float, DateTime
from sqlalchemy.ext.declarative import declarative_base
from sqlalchemy.orm import backref, relationship

Base = declarative_base()


class RunMetadata(Base):
    __tablename__ = 'RunMetadata'

    run_id = Column(Integer, primary_key=True)
    run_name = Column(String(250), nullable=False, index=True)
    watershed_id = Column(Integer,
                          ForeignKey('Watershed.watershed_id'),
                          nullable=True)
    pixel = Column(Integer, nullable=True)
    description = Column(String(250), nullable=True)
    smrf_version = Column(String(250), nullable=True)
    awsm_version = Column(String(250), nullable=True)
    snowav_version = Column(String(250), nullable=True)
    data_type = Column(String(250), nullable=True)
    data_location = Column(String(1200), nullable=True)
    file_type = Column(String(250), nullable=True)
    config_file = Column(String(250), nullable=True)
    proc_time = Column(DateTime, nullable=True)

    # This puts Watershed fields available through Basin at Basin.watershed
    run_metadata = relationship('Watershed',
                                backref=backref('RunMetadata', lazy='dynamic'))


class Watershed(Base):
    __tablename__ = 'Watershed'

    watershed_id = Column(Integer, primary_key=True)
    watershed_name = Column(String(250), unique=True, nullable=False,
                            index=True)
    shapefile = Column(types.LargeBinary, nullable=True)


class Basin(Base):
    __tablename__ = 'Basin'

    watershed_id = Column(Integer,
                          ForeignKey('Watershed.watershed_id',
                                     onupdate="cascade"),
                          nullable=False)
    basin_id = Column(Integer, primary_key=True, nullable=False)
    basin_name = Column(String(250), nullable=False)
    shapefile = Column(types.LargeBinary, nullable=True)

    watershed = relationship('Watershed',
                             backref=backref('Basin', lazy='dynamic'))


class Results(Base):
    __tablename__ = 'Results'

    id = Column(Integer, primary_key=True)
    basin_id = Column(Integer,
                      ForeignKey('Basin.basin_id', onupdate="cascade"),
                      index=True)
    run_id = Column(Integer,
                    ForeignKey('RunMetadata.run_id', onupdate="cascade"),
                    index=True)
    date_time = Column(DateTime, nullable=False, index=True)
    variable = Column(String(250), nullable=True, index=True)
    variable_id = Column(Integer,
                         ForeignKey('VariableUnits.id', onupdate="cascade"),
                         nullable=True)
    value = Column(Float, nullable=True)
    elevation = Column(String(250), nullable=True)

    variable_units = relationship('VariableUnits',
                                  backref=backref('Results'))

    runmetadata = relationship('RunMetadata',
                               backref=backref('Results', lazy='dynamic'))

    basinid = relationship('Basin',
                           backref=backref('Results', lazy='dynamic'))


class VariableUnits(Base):
    __tablename__ = 'VariableUnits'

    id = Column(Integer, primary_key=True, autoincrement=True)
    run_id = Column(Integer,
                    ForeignKey('RunMetadata.run_id', onupdate="cascade"),
                    index=True)
    variable = Column(String(250), nullable=True)
    unit = Column(String(250), nullable=True)
    name = Column(String(250), nullable=True)

    run_metadata = relationship('RunMetadata',
                                backref=backref('VariableUnits', lazy='dynamic'))


class Inputs(Base):
    __tablename__ = 'Inputs'

    id = Column(Integer, primary_key=True)
    run_id = Column(Integer, nullable=False)
    run_name = Column(String(250), nullable=False, index=True)
    basin_id = Column(Integer, nullable=False, index=True)
    date_time = Column(DateTime, nullable=False, index=True)
    variable = Column(String(250), nullable=False)
    function = Column(String(250), nullable=False)
    value = Column(Float, nullable=True)
    unit = Column(String(250), nullable=True)


class Pixels(Base):
    __tablename__ = 'Pixels'

    id = Column(Integer, primary_key=True, unique=True)
    location = Column(String(250), nullable=False)
    model_row = Column(Integer, nullable=False)
    model_col = Column(Integer, nullable=False)
    description = Column(String(250), nullable=True)
    name = Column(String(250), nullable=True)
    utm_x = Column(Float, nullable=False)
    utm_y = Column(Float, nullable=False)
    elevation = Column(Float, nullable=False)


class PixelsData(Base):
    __tablename__ = 'PixelsData'

    pid = Column(Integer, primary_key=True)
    pixel_id = Column(Integer, ForeignKey('Pixels.id'), nullable=False)
    date_time = Column(DateTime, nullable=False, index=True)
    pixels = relationship(Pixels, foreign_keys=[pixel_id])

    # smrf
    air_temp = Column(Float, nullable=True)
    albedo_vis = Column(Float, nullable=True)
    albedo_ir = Column(Float, nullable=True)
    clear_ir_beam = Column(Float, nullable=True)
    clear_ir_diffuse = Column(Float, nullable=True)
    cloud_factor = Column(Float, nullable=True)
    cloud_ir_beam = Column(Float, nullable=True)
    cloud_ir_diffuse = Column(Float, nullable=True)
    cloud_vis_beam = Column(Float, nullable=True)
    cloud_vis_diffuse = Column(Float, nullable=True)
    dew_point = Column(Float, nullable=True)
    flatwind = Column(Float, nullable=True)
    net_solar = Column(Float, nullable=True)
    percent_snow = Column(Float, nullable=True)
    precip = Column(Float, nullable=True)
    precip_temp = Column(Float, nullable=True)
    snow_density = Column(Float, nullable=True)
    storm_days = Column(Float, nullable=True)
    thermal = Column(Float, nullable=True)
    thermal_clear = Column(Float, nullable=True)
    thermal_veg = Column(Float, nullable=True)
    thermal_cloud = Column(Float, nullable=True)
    vapor_pressure = Column(Float, nullable=True)
    veg_ir_diffuse = Column(Float, nullable=True)
    veg_ir_beam = Column(Float, nullable=True)
    veg_vis_beam = Column(Float, nullable=True)
    veg_vis_diffuse = Column(Float, nullable=True)
    wind_direction = Column(Float, nullable=True)
    wind_speed = Column(Float, nullable=True)

    # em.nc
    net_rad = Column(Float, nullable=True)
    sensible_heat = Column(Float, nullable=True)
    latent_heat = Column(Float, nullable=True)
    snow_soil = Column(Float, nullable=True)
    precip_advected = Column(Float, nullable=True)
    sum_EB = Column(Float, nullable=True)
    evaporation = Column(Float, nullable=True)
    snowmelt = Column(Float, nullable=True)
    SWI = Column(Float, nullable=True)
    cold_content = Column(Float, nullable=True)

    # snow.nc
    thickness = Column(Float, nullable=True)
    density = Column(Float, nullable=True)
    specific_mass = Column(Float, nullable=True)
    liquid_water = Column(Float, nullable=True)
    temp_surf = Column(Float, nullable=True)
    temp_lower = Column(Float, nullable=True)
    temp_snowcover = Column(Float, nullable=True)
    thickness_lower = Column(Float, nullable=True)
    water_saturation = Column(Float, nullable=True)
