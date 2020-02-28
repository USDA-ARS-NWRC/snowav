
from collections import OrderedDict
from copy import deepcopy
import numpy as np
import pandas as pd

class AwsmInputsOutputs(object):
    def __init__(self):
        """ smrf input and iSnobal output variable definitions for processing,
        including bands, units, database table, and derivatives.
        """

        self.vars = {'depth':
                            {'band': 0,
                            'derivatives': {'requires': [], 'products': []},
                            'unit_type': 'depth',
                            'description': 'snow depth',
                            'units': 'm',
                            'calculate': 'mean',
                            'table': 'Results',
                            'file': 'snow.nc'},
                     'density':
                            {'band': 1,
                            'derivatives': {'requires': [], 'products': []},
                            'unit_type': 'density',
                            'description': 'density',
                            'units': 'kg/m^3',
                            'calculate': 'mean',
                            'table': 'Results',
                            'file': 'snow.nc'},
                     'swe_z':
                            {'band': 2,
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
                             'derivatives': {'requires': [], 'products': []},
                             'unit_type': 'depth',
                             'description': 'liquid water content',
                             'units':'mm',
                             'calculate': 'mean',
                             'table': 'Results',
                             'file': 'snow.nc'},
                     'temp_surface':
                             {'band': 4,
                             'derivatives': {'requires': [], 'products': []},
                             'unit_type': 'temperature',
                             'description': 'surface layer temperature',
                             'units': 'C',
                             'calculate': 'mean',
                             'table': 'Results',
                             'file': 'snow.nc'},
                     'temp_lower':
                             {'band': 5,
                             'derivatives': {'requires': [], 'products': []},
                             'unit_type': 'temperature',
                             'description': 'lower layer temperature',
                             'units': 'C',
                             'calculate': 'mean',
                             'table': 'Results',
                             'file': 'snow.nc'},
                     'temp_bulk':
                             {'band': 6,
                             'derivatives': {'requires': [], 'products': []},
                             'unit_type': 'temperature',
                             'description': 'bulk snowpack temperature',
                             'units': 'C',
                             'calculate': 'mean',
                             'table': 'Results',
                             'file': 'snow.nc'},
                     'depth_lower_layer':
                             {'band': 7,
                             'derivatives': {'requires': [], 'products': []},
                             'unit_type': 'depth',
                             'description': 'lower layer depth',
                             'units': 'm',
                             'calculate': 'mean',
                             'table': 'Results',
                             'file': 'snow.nc'},
                     'h20_sat':
                             {'band': 8,
                             'derivatives': {'requires': [], 'products': []},
                             'unit_type': 'percent',
                             'description': '% liquid water saturation',
                             'units': 'percent',
                             'calculate': 'mean',
                             'table': 'Results',
                             'file': 'snow.nc'},
                     'R_n':
                             {'band': 0,
                             'derivatives': {'requires': [], 'products': []},
                             'unit_type': 'flux',
                             'description': 'net all-wave radiation',
                             'units': 'W/m^2',
                             'calculate': 'mean',
                             'table': 'Results',
                             'file': 'em.nc'},
                     'H':
                             {'band': 1,
                             'derivatives': {'requires': [], 'products': []},
                             'unit_type': 'flux',
                             'description': 'sensible heat transfer',
                             'units': 'W/m^2',
                             'calculate': 'mean',
                             'table': 'Results',
                             'file': 'em.nc'},
                     'L_v_E':
                             {'band': 2,
                             'derivatives': {'requires': [], 'products': []},
                             'unit_type': 'flux',
                             'description': 'latent heat exchange',
                             'units': 'W/m^2',
                             'calculate': 'mean',
                             'table': 'Results',
                             'file': 'em.nc'},
                     'G':
                             {'band': 3,
                             'derivatives': {'requires': [], 'products': []},
                             'unit_type': 'flux',
                             'description': 'snow/soil heat exchange',
                             'units': 'W/m^2',
                             'calculate': 'mean',
                             'table': 'Results',
                             'file': 'em.nc'},
                     'M':
                             {'band': 4,
                             'derivatives': {'requires': [], 'products': []},
                             'unit_type': 'flux',
                             'description': 'advected heat from precip',
                             'units': 'W/m^2',
                             'calculate': 'mean',
                             'table': 'Results',
                             'file': 'em.nc'},
                     'delta_Q':
                             {'band': 5,
                             'derivatives': {'requires': [], 'products': []},
                             'unit_type': 'flux',
                             'description': 'sum of e.b. terms',
                             'units': 'W/m^2',
                             'calculate': 'mean',
                             'table': 'Results',
                             'file': 'em.nc'},
                     'evap_z':
                             {'band': 6,
                             'derivatives': {'requires': [], 'products': []},
                             'unit_type': 'depth',
                             'description': 'evaporation',
                             'units': 'mm',
                             'calculate': 'mean',
                             'table': 'Results',
                             'file': 'em.nc'},
                     'melt':
                             {'band': 7,
                             'derivatives': {'requires': [], 'products': []},
                             'unit_type': 'depth',
                             'description': 'melt',
                             'units': 'mm',
                             'calculate': 'mean',
                             'table': 'Results',
                             'file': 'em.nc'},
                     'swi_z':
                             {'band': 8,
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
                             'derivatives': {'requires': [], 'products': []},
                             'unit_type': 'flux',
                             'description': 'cold content',
                             'units': 'J/m^2',
                             'calculate': 'mean',
                             'table': 'Results',
                             'file': 'em.nc'},
                     'air_temp':
                             {'band': None,
                             'derivatives': {'requires': [], 'products': []},
                             'unit_type': 'temperature',
                             'description': 'air temperature',
                             'units': 'C',
                             'calculate': 'mean',
                             'table': 'Inputs',
                             'file': 'air_temp.nc'},
                     'cloud_factor':
                             {'band': None,
                             'derivatives': {'requires': [], 'products': []},
                             'unit_type': 'dimensionless',
                             'description': 'cloud factor',
                             'units': 'dimensionless',
                             'calculate': 'mean',
                             'table': 'Inputs',
                             'file': 'cloud_factor.nc'},
                     'net_solar':
                             {'band': None,
                             'derivatives': {'requires': [], 'products': []},
                             'unit_type': 'flux',
                             'description': 'net solar',
                             'units': 'W/m^2',
                             'calculate': 'mean',
                             'table': 'Inputs',
                             'file': 'net_solar.nc'},
                     'storm_days':
                             {'band': None,
                             'derivatives': {'requires': [], 'products': []},
                             'unit_type': 'time',
                             'description': 'days since last precip',
                             'units': 'day',
                             'calculate': 'mean',
                             'table': 'Inputs',
                             'file': 'storm_days.nc'},
                     'precip_z':
                             {'band': None,
                             'derivatives': {'requires': ['percent_snow'],
                                             'products': ['precip_vol',
                                                          'rain_z']},
                             'unit_type': 'depth',
                             'description': 'precipitation',
                             'units': 'mm',
                             'calculate': 'mean',
                             'table': 'Results',
                             'file': 'precip.nc'},
                     'percent_snow':
                             {'band': None,
                             'derivatives': {'requires': [],
                                             'products': ['precip_vol',
                                                          'precip_z',
                                                          'rain_z']},
                             'unit_type': 'percent',
                             'description': 'percent of precipitation that is snow',
                             'units': 'mm',
                             'calculate': 'mean',
                             'table': 'Inputs',
                             'file': 'percent_snow.nc'},
                     'snow_density':
                             {'band': None,
                             'derivatives': {'requires': [], 'products': []},
                             'unit_type': 'density',
                             'description': 'density of storm snow',
                             'units': 'kg/m^3',
                             'calculate': 'mean',
                             'table': 'Inputs',
                             'file': 'snow_density.nc'},
                     'thermal':
                             {'band': None,
                             'derivatives': {'requires': [], 'products': []},
                             'unit_type': 'flux',
                             'description': 'thermal radiation',
                             'units': 'W/m^2',
                             'calculate': 'mean',
                             'table': 'Inputs',
                             'file': 'thermal.nc'},
                     'vapor_pressure':
                             {'band': None,
                             'derivatives': {'requires': [], 'products': []},
                             'unit_type': 'pressure',
                             'description': 'vapor pressure',
                             'units': 'Pa',
                             'calculate': 'mean',
                             'table': 'Inputs',
                             'file': 'vapor_pressure.nc'},
                     'wind_speed':
                             {'band': None,
                             'derivatives': {'requires': [], 'products': []},
                             'unit_type': 'speed',
                             'description': 'wind_speed',
                             'units': 'm/s',
                             'calculate': 'mean',
                             'table': 'Inputs',
                             'file': 'smrf'},

                     # these are the derivative products
                     'swi_vol':
                             {'band': None,
                              'derivatives': {'requires': ['swi_z'],
                                              'products': []},
                              'unit_type': 'volumne',
                              'description': 'surface water input volume',
                              'units': 'm^3',
                              'calculate': 'sum',
                              'table': 'Results'},
                     'swe_avail':
                             {'band': None,
                              'derivatives': {'requires': ['swe_z, coldcont'],
                                              'products': []},
                              'unit_type': 'volumne',
                              'description': 'snow water equivalent available for melt',
                              'units': 'm^3',
                              'calculate': 'sum',
                              'table': 'Results'},
                     'swe_unavail':
                             {'band': None,
                              'derivatives': {'requires': ['swe_z, coldcont'],
                                              'products': []},
                              'unit_type': 'volumne',
                              'description': 'snow water equivalent unavailable for melt',
                              'units': 'm^3',
                              'calculate': 'sum',
                              'table': 'Results'},
                     'swe_vol':
                             {'band': None,
                              'derivatives': {'requires': ['swe_z'],
                                              'products': []},
                              'unit_type': 'volumne',
                              'description': 'snow water equivalent volume',
                              'units': 'm^3',
                              'calculate': 'sum',
                              'table': 'Results'},
                     'precip_vol':
                             {'band': None,
                              'derivatives': {'requires': ['precip_z'],
                                              'products': []},
                              'unit_type': 'volume',
                              'description': 'precipitation volume',
                              'units': 'm^3',
                              'calculate': 'sum',
                              'table': 'Results'},
                     'rain_z':
                             {'band': None,
                              'derivatives': {'requires': ['precip_z',
                                                           'percent_snow'],
                                              'products': []},
                              'unit_type': 'depth',
                              'description': 'rain depth',
                              'units': 'mm',
                              'calculate': 'mean',
                              'table': 'Results'},
                     'snow_line':
                             {'band': None,
                              'derivatives': {'requires': ['swe_z'],
                                              'products': []},
                              'unit_type': 'elevation',
                              'description': 'snow line',
                              'units': 'm',
                              'calculate': 'mean',
                              'table': 'Results'}
                          }

        # all named snowav variables
        self.snowav_named_variables = []

        for v in self.vars.keys():
            self.snowav_named_variables.append(v)

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

        # make a variable for each desired from the master list, and add
        # an empty dataframe for results
        for p in properties:

            # catch invalid properties
            if p not in list(self.vars.keys()):
                raise ValueError('{} not in available snowav '
                                 'properties'.format(p))

            self.variables[p] = deepcopy(self.vars[p])
            self.variables[p]['df'] = pd.DataFrame(np.nan,
                                                   index = edges,
                                                   columns = columns)

            # add derivates as well
            if self.vars[p]['derivatives']['products'] != []:
                for d in self.vars[p]['derivatives']['products']:
                    self.variables[d] = deepcopy(self.vars[d])
                    self.variables[d]['df'] = pd.DataFrame(np.nan,
                                                           index = edges,
                                                           columns = columns)
