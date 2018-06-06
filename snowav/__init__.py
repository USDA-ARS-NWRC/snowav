import os
__core_config__ = os.path.abspath(os.path.dirname(__file__) + '/config/CoreConfig.ini')
__recipes__ = os.path.abspath(os.path.dirname(__file__) + '/config/recipes.ini')

__config_titles__ = {'Basin': 'Overview information for SNOWAV, including'
                              ' file names, save directory, and basin'
                              ' summary information',
                      'Outputs': 'Start and end files for report period',
                      'Runs': 'iSnobal output run directories to include'
                              ' in processing',
                      'Validate': 'Stations and labels for snow pillow'
                              ' validation figure',
                      'Plots': 'Parameters for figure size',
                      'Masks': 'Paths and labels for subbasins and DEM',
                      'Report': 'Paths for report LaTex template and summary'
                              ' files'
                      }
from . import plotting
from . import methods
from . import report
from . import utils

__config_header__ =# utils.utilities.get_config_header()
