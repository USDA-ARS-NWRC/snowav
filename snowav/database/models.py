from collections import OrderedDict
from copy import deepcopy
import numpy as np
import pandas as pd


class AwsmInputsOutputs(object):
    """ smrf input and iSnobal output variable definitions for processing,
    including bands, units, database table, and derivatives. """

    def __init__(self):

        self.vars = {'depth':
                         {'band': 0,
                          'nc_name': 'thickness',
                          'derivatives': {'requires': [], 'products': []},
                          'unit_type': 'depth',
                          'description': 'snow depth',
                          'units': 'm',
                          'calculate': 'mean',
                          'table': 'Results',
                          'file': 'snow.nc'},
                     'density':
                         {'band': 1,
                          'nc_name': 'snow_density',
                          'derivatives': {'requires': [], 'products': []},
                          'unit_type': 'density',
                          'description': 'density',
                          'units': 'kg/m^3',
                          'calculate': 'mean',
                          'table': 'Results',
                          'file': 'snow.nc'},
                     'swe_z':
                         {'band': 2,
                          'nc_name': 'specific_mass',
                          'derivatives': {'requires': ['coldcont'],
                                          'products': ['swe_vol', 'swe_avail',
                                                       'swe_unavail']},
                          'unit_type': 'depth',
                          'description': 'snow water equivalent depth',
                          'units': 'mm',
                          'calculate': 'mean',
                          'table': 'Results',
                          'file': 'snow.nc'},
                     'lwc': {'band': 3,
                             'nc_name': 'liquid_water',
                             'derivatives': {'requires': [], 'products': []},
                             'unit_type': 'depth',
                             'description': 'liquid water content',
                             'units': 'mm',
                             'calculate': 'mean',
                             'table': 'Results',
                             'file': 'snow.nc'},
                     'temp_surface':
                         {'band': 4,
                          'nc_name': 'temp_surf',
                          'derivatives': {'requires': [], 'products': []},
                          'unit_type': 'temperature',
                          'description': 'surface layer temperature',
                          'units': 'C',
                          'calculate': 'mean',
                          'table': 'Results',
                          'file': 'snow.nc'},
                     'temp_lower':
                         {'band': 5,
                          'nc_name': 'temp_lower',
                          'derivatives': {'requires': [], 'products': []},
                          'unit_type': 'temperature',
                          'description': 'lower layer temperature',
                          'units': 'C',
                          'calculate': 'mean',
                          'table': 'Results',
                          'file': 'snow.nc'},
                     'temp_bulk':
                         {'band': 6,
                          'nc_name': 'temp_snowcover',
                          'derivatives': {'requires': [], 'products': []},
                          'unit_type': 'temperature',
                          'description': 'bulk snowpack temperature',
                          'units': 'C',
                          'calculate': 'mean',
                          'table': 'Results',
                          'file': 'snow.nc'},
                     'depth_lower_layer':
                         {'band': 7,
                          'nc_name': 'thickness_lower',
                          'derivatives': {'requires': [], 'products': []},
                          'unit_type': 'depth',
                          'description': 'lower layer depth',
                          'units': 'm',
                          'calculate': 'mean',
                          'table': 'Results',
                          'file': 'snow.nc'},
                     'h20_sat':
                         {'band': 8,
                          'nc_name': 'water_saturation',
                          'derivatives': {'requires': [], 'products': []},
                          'unit_type': 'percent',
                          'description': '% liquid water saturation',
                          'units': '%',
                          'calculate': 'mean',
                          'table': 'Results',
                          'file': 'snow.nc'},
                     'R_n':
                         {'band': 0,
                          'nc_name': 'net_rad',
                          'derivatives': {'requires': [], 'products': []},
                          'unit_type': 'flux',
                          'description': 'net all-wave radiation',
                          'units': 'W/m^2',
                          'calculate': 'mean',
                          'table': 'Results',
                          'file': 'em.nc'},
                     'H':
                         {'band': 1,
                          'nc_name': 'sensible_heat',
                          'derivatives': {'requires': [], 'products': []},
                          'unit_type': 'flux',
                          'description': 'sensible heat transfer',
                          'units': 'W/m^2',
                          'calculate': 'mean',
                          'table': 'Results',
                          'file': 'em.nc'},
                     'L_v_E':
                         {'band': 2,
                          'nc_name': 'latent_heat',
                          'derivatives': {'requires': [], 'products': []},
                          'unit_type': 'flux',
                          'description': 'latent heat exchange',
                          'units': 'W/m^2',
                          'calculate': 'mean',
                          'table': 'Results',
                          'file': 'em.nc'},
                     'G':
                         {'band': 3,
                          'nc_name': 'snow_soil',
                          'derivatives': {'requires': [], 'products': []},
                          'unit_type': 'flux',
                          'description': 'snow/soil heat exchange',
                          'units': 'W/m^2',
                          'calculate': 'mean',
                          'table': 'Results',
                          'file': 'em.nc'},
                     'M':
                         {'band': 4,
                          'nc_name': 'precip_advected',
                          'derivatives': {'requires': [], 'products': []},
                          'unit_type': 'flux',
                          'description': 'advected heat from precip',
                          'units': 'W/m^2',
                          'calculate': 'mean',
                          'table': 'Results',
                          'file': 'em.nc'},
                     'delta_Q':
                         {'band': 5,
                          'nc_name': 'sum_EB',
                          'derivatives': {'requires': [], 'products': []},
                          'unit_type': 'flux',
                          'description': 'sum of e.b. terms',
                          'units': 'W/m^2',
                          'calculate': 'mean',
                          'table': 'Results',
                          'file': 'em.nc'},
                     'evap_z':
                         {'band': 6,
                          'nc_name': 'evaporation',
                          'derivatives': {'requires': [], 'products': []},
                          'unit_type': 'depth',
                          'description': 'evaporation',
                          'units': 'mm',
                          'calculate': 'mean',
                          'table': 'Results',
                          'file': 'em.nc'},
                     'melt':
                         {'band': 7,
                          'nc_name': 'snowmelt',
                          'derivatives': {'requires': [], 'products': []},
                          'unit_type': 'depth',
                          'description': 'melt',
                          'units': 'mm',
                          'calculate': 'mean',
                          'table': 'Results',
                          'file': 'em.nc'},
                     'swi_z':
                         {'band': 8,
                          'nc_name': 'SWI',
                          'derivatives': {'requires': [],
                                          'products': ['swi_vol']},
                          'unit_type': 'depth',
                          'description': 'surface water input depth',
                          'units': 'mm',
                          'calculate': 'mean',
                          'table': 'Results',
                          'file': 'em.nc'},
                     'coldcont':
                         {'band': 9,
                          'nc_name': 'cold_content',
                          'derivatives': {'requires': [], 'products': []},
                          'unit_type': 'flux',
                          'description': 'cold content',
                          'units': 'J/m^2',
                          'calculate': 'mean',
                          'table': 'Results',
                          'file': 'em.nc'},
                     'air_temp':
                         {'band': None,
                          'nc_name': 'air_temp',
                          'derivatives': {'requires': [], 'products': []},
                          'unit_type': 'temperature',
                          'description': 'air temperature',
                          'units': 'C',
                          'calculate': 'mean',
                          'table': 'Inputs',
                          'file': 'air_temp.nc'},
                     'cloud_factor':
                         {'band': None,
                          'nc_name': 'cloud_factor',
                          'derivatives': {'requires': [], 'products': []},
                          'unit_type': 'dimensionless',
                          'description': 'cloud factor',
                          'units': 'unitless',
                          'calculate': 'mean',
                          'table': 'Inputs',
                          'file': 'cloud_factor.nc'},
                     'net_solar':
                         {'band': None,
                          'nc_name': 'net_solar',
                          'derivatives': {'requires': [], 'products': []},
                          'unit_type': 'flux',
                          'description': 'net solar',
                          'units': 'W/m^2',
                          'calculate': 'mean',
                          'table': 'Inputs',
                          'file': 'net_solar.nc'},
                     'storm_days':
                         {'band': None,
                          'nc_name': 'storm_days',
                          'derivatives': {'requires': [], 'products': []},
                          'unit_type': 'time',
                          'description': 'days since last precip',
                          'units': 'day',
                          'calculate': 'mean',
                          'table': 'Inputs',
                          'file': 'storm_days.nc'},
                     'precip_z':
                         {'band': None,
                          'nc_name': None,
                          'derivatives': {'requires': ['percent_snow'],
                                          'products': ['precip_vol',
                                                       'rain_z']},
                          'unit_type': 'depth',
                          'description': 'precipitation',
                          'units': 'mm',
                          'calculate': 'mean',
                          'table': 'Results',
                          'file': 'precip.nc'},
                     # currently we have 'precip' that goes to Inputs table
                     # and 'precip_z' on Results
                     'precip':
                         {'band': None,
                          'nc_name': 'precip',
                          'derivatives': {'requires': [],
                                          'products': []},
                          'unit_type': 'depth',
                          'description': 'precipitation',
                          'units': 'mm',
                          'calculate': 'mean',
                          'table': 'Inputs',
                          'file': 'precip.nc'},

                     'percent_snow':
                         {'band': None,
                          'nc_name': 'percent_snow',
                          'derivatives': {'requires': [],
                                          'products': ['precip_vol',
                                                       'precip_z',
                                                       'rain']},
                          'unit_type': 'percent',
                          'description': 'percent of precipitation that is snow',
                          'units': '%',
                          'calculate': 'mean',
                          'table': 'Inputs',
                          'file': 'percent_snow.nc'},
                     'snow_density':
                         {'band': None,
                          'nc_name': 'snow_density',
                          'derivatives': {'requires': [], 'products': []},
                          'unit_type': 'density',
                          'description': 'density of storm snow',
                          'units': 'kg/m^3',
                          'calculate': 'mean',
                          'table': 'Inputs',
                          'file': 'snow_density.nc'},
                     'thermal':
                         {'band': None,
                          'nc_name': 'thermal',
                          'derivatives': {'requires': [], 'products': []},
                          'unit_type': 'flux',
                          'description': 'thermal radiation',
                          'units': 'W/m^2',
                          'calculate': 'mean',
                          'table': 'Inputs',
                          'file': 'thermal.nc'},
                     'vapor_pressure':
                         {'band': None,
                          'nc_name': 'vapor_pressure',
                          'derivatives': {'requires': [], 'products': []},
                          'unit_type': 'pressure',
                          'description': 'vapor pressure',
                          'units': 'Pa',
                          'calculate': 'mean',
                          'table': 'Inputs',
                          'file': 'vapor_pressure.nc'},
                     'wind_speed':
                         {'band': None,
                          'nc_name': 'wind_speed',
                          'derivatives': {'requires': [], 'products': []},
                          'unit_type': 'speed',
                          'description': 'wind_speed',
                          'units': 'm/s',
                          'calculate': 'mean',
                          'table': 'Inputs',
                          'file': 'wind_speed.nc'},

                     # these are the derivative products
                     'swi_vol':
                         {'band': None,
                          'nc_name': None,
                          'derivatives': {'requires': ['swi_z'],
                                          'products': []},
                          'unit_type': 'volumne',
                          'description': 'surface water input volume',
                          'units': 'm^3',
                          'calculate': 'sum',
                          'table': 'Results',
                          'file': None},
                     'swe_avail':
                         {'band': None,
                          'nc_name': None,
                          'derivatives': {'requires': ['swe_z, coldcont'],
                                          'products': []},
                          'unit_type': 'volumne',
                          'description': 'snow water equivalent available for melt',
                          'units': 'm^3',
                          'calculate': 'sum',
                          'table': 'Results',
                          'file': None},
                     'swe_unavail':
                         {'band': None,
                          'nc_name': None,
                          'derivatives': {'requires': ['swe_z, coldcont'],
                                          'products': []},
                          'unit_type': 'volumne',
                          'description': 'snow water equivalent unavailable for melt',
                          'units': 'm^3',
                          'calculate': 'sum',
                          'table': 'Results',
                          'file': None},
                     'swe_vol':
                         {'band': None,
                          'nc_name': None,
                          'derivatives': {'requires': ['swe_z'],
                                          'products': []},
                          'unit_type': 'volumne',
                          'description': 'snow water equivalent volume',
                          'units': 'm^3',
                          'calculate': 'sum',
                          'table': 'Results',
                          'file': None},
                     'precip_vol':
                         {'band': None,
                          'nc_name': None,
                          'derivatives': {'requires': ['precip_z'],
                                          'products': []},
                          'unit_type': 'volume',
                          'description': 'precipitation volume',
                          'units': 'm^3',
                          'calculate': 'sum',
                          'table': 'Results',
                          'file': None},
                     'rain_z':
                         {'band': None,
                          'nc_name': None,
                          'derivatives': {'requires': ['precip_z',
                                                       'percent_snow'],
                                          'products': []},
                          'unit_type': 'depth',
                          'description': 'rain depth',
                          'units': 'mm',
                          'calculate': 'mean',
                          'table': 'Results',
                          'file': None},
                     'snow_line':
                         {'band': None,
                          'nc_name': None,
                          'derivatives': {'requires': ['swe_z'],
                                          'products': []},
                          'unit_type': 'elevation',
                          'description': 'snow line',
                          'units': 'm',
                          'calculate': 'mean',
                          'table': 'Results',
                          'file': None}
                     }

        # snow.nc and em.nc variables
        self.awsm_variables = []

        # smrf variables
        self.smrf_variables = []

        # all named snowav variables
        self.snowav_named_variables = []

        # nc variables
        self.nc_variables = []

        for v in self.vars.keys():
            self.snowav_named_variables.append(v)
            if self.vars[v]['band'] is not None:
                self.awsm_variables.append(v)

            if (self.vars[v]['file'] is not None and
                    self.vars[v]['file'] not in ['em.nc', 'snow.nc']):
                self.smrf_variables.append(v)

            if self.vars[v]['nc_name'] is not None:
                self.nc_variables.append(self.vars[v]['nc_name'])

        # Results table values
        self.snowav_results_variables = []
        for v in self.vars.keys():
            if self.vars[v]['table'] == 'Results':
                self.snowav_results_variables.append(v)
                if self.vars[v]['derivatives']['products'] != []:
                    for d in self.vars[v]['derivatives']['products']:
                        self.snowav_results_variables.append(d)

        # Inputs table values
        self.snowav_inputs_variables = []
        for v in self.vars.keys():
            if self.vars[v]['table'] == 'Inputs':
                self.snowav_inputs_variables.append(v)
                if self.vars[v]['derivatives']['products'] != []:
                    for d in self.vars[v]['derivatives']['products']:
                        self.snowav_inputs_variables.append(d)

        self.cumulative_sum_variables = ['swi_vol', 'precip_vol', 'precip_z',
                                         'swi_z']

        self.process_depth_units = ['coldcont', 'density', 'depth', 'evap_z',
                                    'L_v_E', 'lwc', 'temp_surface',
                                    'temp_lower', 'temp_bulk',
                                    'depth_lower_layer', 'h20_sat', 'R_n', 'H',
                                    'L_v_E', 'G', 'M', 'delta_Q']

    def make_variables(self, properties, edges, columns):
        """ Make a subset of vars with only the requested properties.

        This uses an OrderedDict, and expects properties to be in the
        correct processing order.

        Args
        ------
            properties: list, composed of keys from self.vars
            edges: list, elevation bins
            columns: list, basin names
        """

        self.variables = OrderedDict()
        self.snowav_inputs_variables = []
        self.snowav_results_variables = []

        # make a variable for each desired from the master list, and add
        # an empty dataframe for results
        for p in properties:

            # make lists for each table
            # inputs are processed first, but precip gets placed on Results
            # table
            if self.vars[p]['table'] == 'Inputs':
                self.snowav_inputs_variables.append(p)
                for d in self.vars[p]['derivatives']['products']:
                    self.snowav_inputs_variables.append(d)

            if self.vars[p]['table'] == 'Results':
                self.snowav_results_variables.append(p)
                for d in self.vars[p]['derivatives']['products']:
                    self.snowav_results_variables.append(d)

            # catch invalid properties
            if p not in list(self.vars.keys()):
                raise ValueError('{} not in available snowav '
                                 'properties'.format(p))

            self.variables[p] = deepcopy(self.vars[p])

            if p == 'snow_line':
                index = ['total']
            else:
                index = edges

            self.variables[p]['df'] = pd.DataFrame(np.nan,
                                                   index=index,
                                                   columns=columns)

            # add derivates as well
            if self.vars[p]['derivatives']['products'] != []:
                for d in self.vars[p]['derivatives']['products']:
                    self.variables[d] = deepcopy(self.vars[d])
                    self.variables[d]['df'] = pd.DataFrame(np.nan,
                                                           index=edges,
                                                           columns=columns)
