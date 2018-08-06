import os
__version__ = '0.3.0'
__core_config__ = os.path.abspath(os.path.dirname(__file__) + '/config/CoreConfig.ini')
__recipes__ = os.path.abspath(os.path.dirname(__file__) + '/config/recipes.ini')

__config_titles__ = {'basin': 'Overview information for SNOWAV, including'
                              ' file names, save directory, and basin'
                              ' summary information',
                      'outputs': 'Start and end files for report period',
                      'runs': 'iSnobal output run directories to include'
                              ' in processing',
                      'validate': 'Stations and labels for snow pillow'
                              ' validation figure',
                      'plots': 'Parameters for figure size',
                      'masks': 'Paths and labels for subbasins and DEM',
                      'report': 'Paths for report LaTex template and summary'
                              ' files'
                      }
from . import framework
from . import database
from . import plotting
from . import methods
from . import report
from . import utils
__config_header__ = utils.utilities.get_config_header()
