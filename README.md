# SNOw and Water model Analysis and Visualization

[![GitHub version](https://badge.fury.io/gh/USDA-ARS-NWRC%2Fsnowav.svg)](https://badge.fury.io/gh/USDA-ARS-NWRC%2Fsnowav)
[![Build Status](https://travis-ci.org/USDA-ARS-NWRC/snowav.svg?branch=devel)](https://travis-ci.org/USDA-ARS-NWRC/snowav)

SNOWAV was developed at the USDA Agricultural Research Service in Boise, Idaho, and processes AWSM model outputs into formats and figures for use by water resource managers. See CoreConfig.ini for details on config options.

![image](https://raw.githubusercontent.com/USDA-ARS-NWRC/awsm/master/docs/_static/ModelSystemOverview_new.png)

## Install
The user can either clone this repository or pip install from PyPi.

## Requirements
Currently snowav requires:
- awsm model results in awsm_daily format, including output files in the paths ```.../runs/runYYYYMMDD/snow.nc``` and ```.../runs/runYYYYMMDD/em.nc```, and input files in the paths ```.../data/dataYYYYMMDD/smrfOutputs/```
- topo.nc files that have been created by the basin_setup package
- correct date tags in all snow.nc files
- input files in the paths ```.../data/dataYYYYMMDD/smrfOutputs/precip.nc``` and ```.../data/dataYYYYMMDD/smrfOutputs/percent_snow.nc```. If these do not exist, any precip figures will be set to False.
- the ```lidar_depths_wyYYYY.nc``` file, specified in the config file in [plots] update_file, contains all relevant flights that have been flown. Flights in this file can be removed from processing using [plots] *update_numbers*, but they can't be added.

## Standard Use
The standard snowav model processing, database placement, and figures run is:
```
$ snowav -f config.ini
```

Below are notes for some of the config file fields, for additional information for all options see ```snowav/config/CoreConfig.ini```.

#### [snowav]
```save_path:``` Path to save results to.

```directory:``` Name of folder that will be created under ```save_path```.

```run_name:``` **This is important!** The ```run_name``` field is how results get stored on the database, so only one ```run_name``` field can exist per basin and date. Results will be overwritten if ```[database] overwrite: True```.

```masks:``` This can be used to subset masks that appear in the ```topo.nc``` file in ```dempath```, but is not necessary. If used, the names must match what is found in ```topo.nc```.

```plotlabels:```  This can be used to make different plot labels for the subbasins listed in ```topo.nc```. This must be used with a ```masks``` list.


#### [run]
```directory:``` Load everything from this directory, must end with ```.../runs```. Can subset a date range by using ```start_date``` and ```end_date``` if desired. Can also give a list of single directories that end with ```.../runs/runYYYYMMDD``` and set ```all_subdirs: False```.

```start_date:``` Subsets ```directory```.

```end_date:``` Subsets ```directory```. Can also include a command line call to overwrite this field, which is used in docker-airflow application

```
$ snowav -f config.ini -end_date "2019-12-30 23:00"
```

#### [database]

```overwrite:``` If True, will overwrite database values with same basin and ```run_name```. Set to True if results need to be updated, otherwise leave False for processing speed up.

```convert_ws:``` Converts watershed names to original snowav database version from wy2019, set to True for all but the Tuolumne basin.

#### [diagnostics]
```diagnostics:``` Creates a simple model state diagnostics figure.

```inputs_table:``` Includes smrf field output summaries on the snowav database.

#### [plots]
Turn figures on and off. Some figures require additional fields in other sections to be set in order to work, including ```point_values```, ```inputs```, and ```stn_validate```.

```update_file:``` If a path to a  ```lidar_depths.nc``` is supplied, flight difference figures will be made.

Currently the figures listed below are created by default in both the ```[plots]``` and ```[reports]``` sections. If any of them are set to *False* in ```[plots]``` they will be set to *False* in ```[reports]```. For additional figures to be added to a report they must be *True* in both ```[plots]``` and ```[reports]```.
  - swi
  - image_change
  - cold_content  
  - swe_volume
  - basin_total
  - precip_depth

#### [report]
Figures and tables that are included in the report can be turned on and off here, but they need to be created in ```[plots]``` if you want them in the report.

#### [inflow]
```inflow: True``` triggers reading in operator-generated inflow and the inflow figure. Currently only configured to work with the Tuolumne, and 'FORM11' excel sheets. San Joaquin requires editing of the .xls sheets.


## Additional Features

To process a single snow.nc file and display simple SWE volume figure, without putting results on a database:

```
$ snowav -t <topo.nc> -d <snow.nc>
```

By default this will display SWE volume and save the figure in the snowav repo. Specify a figure save location with ```-p <path>``` and include ``` -s ``` to display the figure.

For a simple difference between two snow.nc files:

```
$ snowav -t <topo.nc> -A <snow.nc> -B <snow.nc>
```

To query and output existing database records, either to the terminal or csv, see CoreConfig.ini ```[query]``` section and example below. If ```[query] query: True``` no other functionality other than the query will be run. This currently requires database login created from ```[database]``` ```mysql: snowav```.
```
[query]
query:              True
basins:             Kaweah River Basin, North Fork
value:              swe_vol
run_name:           kaweah_wy2019_ops
print_all_runs:     False
start_date:         2019-3-15
end_date:           2019-3-17
total:              True
output:             print
csv_base_path:      /mnt/volumes/wkspace/projects/csv_output/
database:           mysql+mysqlconnector://<user>:<pwd>@172.17.0.2/snowav
```

## Notes and Considerations
- **Be deliberate with the [snowav] *run_name* field!** This is how processing runs are identified on the database, and existing records with the same ```run_name```, basin, and date range will be deleted and replaced during re-runs of the same snow.nc files. The ```run_name``` field should, in most cases and normal use, be connected with a specific ```[run]``` ```directory```. We suggest using fields such as *tuol_wy2019_ops* and re-running snowav with *tuol_wy2019_ops* over the proper directory if modifications are made mid-season.

- If depth updates have been applied, year-to-date precipitation values in report Table 1 are no longer valid, although they will still appear. This table includes precipitation and rain derived from the HRRR inputs and does not account for mass that may be added or removed via snow depth updates.

- To help debugging the pdf report if it fails rendering to latex, set config option ```[report]``` ```print_latex: True```. This will print the full latex file to the screen immediately prior to rendering.

- The SWE pillow validation figure ```[plots]``` ```stn_validate``` requires SWE fields on the existing weather database and ```[validate]``` fields to be filled in.

- Topo.nc, masks, and plot labels <br/>
```[snowav] masks:``` can be left blank, and will default to the long_name fields in the ```topo.nc``` file. To subset the number of basins processed and plotted, use ```masks:``` with a list: <br/> *masks: San Joaquin River Basin, Main, South Fork* <br/>
To replace the plot labels, use the *plotlabels* field in combination with *masks*: <br/> *masks: San Joaquin River Basin, Main, South Fork* <br/> *plotlabels: San Joaquin, Mammoth, South Fork*


## SNOWAV Processing Utility
The snowav processing utility snowav.utils.utilities.calculate can be used to make simple calculations on snow.nc and em.nc files. See scripts.sample_process.py for an example.

## DataFrame from Existing Database Records
This simple sample script shows pulling a DataFrame of existing database records. Use this in combination with a snowav query to find and pull records that you want.

```
from datetime import datetime
from snowav.database.database import connect, collect
from snowav.utils.utilities import masks

# These settings will change depending on the basin, time frame, run, and value
# you are interested in
dempath = '/home/ops/wy2019/kings/topo/topo.nc'
start_date = datetime(2019,3,1)
end_date = datetime(2019,7,30)
value = 'swe_vol'
run_name = 'kings_wy2019_ops'

# These are logins for the snowav database, and shouldn't need to change
sql = 'snowav'
user = ''
password = ''
host = '172.17.0.2'
port  = '3306'

value_options = ['swe_vol','swe_z','swi_vol','swi_z','precip_vol','precip_z',
                 'density','coldcont','depth','evap_z']

if value not in value_options:
    raise Exception("'value' must be one of {}".format(value_options))

# Get the list of basins we want by reading the topo.nc file
out = masks(dempath, False)
basin_list = out['plotorder']

# Establish the snowav database connection and get the 'basins' dictionary we
# need for pulling results
basins, cnx, out = connect(None,sql,basin_list,user,password,host,port,True)

# Get snowav database results in a DataFrame
df = collect(cnx, basin_list, basins, start_date, end_date, value, run_name,
             'total', 'daily')
```

## Figures from Existing Database Records
If results have already been processed and put onto a database, figures can be created outside of a snowav processing run (see also scripts/sample_figure.py). See snowav.framework.figures for templates for additional figure creation. Also, if a standard snowav run is processed with [plots] *print_args_dict: True*, the full input dictionary for each figure will be printed to the screen.

```
from datetime import datetime
import netCDF4 as nc
from snowav.plotting.swe_volume import swe_volume
from snowav.plotting.plotlims import plotlims as plotlims
from snowav.database.database import connect, make_session, collect
from snowav.utils.utilities import masks

""" See CoreConfig.ini for more options and details """

connector = 'sqlite:////<database.db>'
dempath = '<topo.nc>'
plotorder = ['Extended Tuolumne', 'Tuolumne', 'Cherry Creek', 'Eleanor']
lims = plotlims(plotorder)
edges = [3000, 4000, 5000, 6000, 7000, 8000, 9000, 10000, 11000, 12000, 13000]

ncf = nc.Dataset('<snow.nc>')
image = ncf['specific_mass'][:]
ncf.close()

out = masks(dempath)
dem = out['dem']
masks = out['masks']
nrows = out['nrows']
ncols = out['ncols']
plotorder = out['plotorder']
labels = out['labels']
start_date = datetime(2019,6,1,23,0,0)
end_date = datetime(2019,6,2,23,0,0)
run_name = 'tuol_wy2019_devel'

barcolors = ['xkcd:cobalt', 'xkcd:mustard green', 'xkcd:lichen',
             'xkcd:pale green', 'xkcd:blue green']

basins, cnx, out = connect(sqlite=connector, mysql=None, plotorder=plotorder,
                           user=None, password=None, host=None, port=None)

args = {'report_start': '',
        'report_date': '',
        'print': False,
        'run_name': run_name,
        'start_date': start_date,
        'end_date': end_date,
        'directory': 'test',
        'figs_path': '/<path>/',
        'edges': edges,
        'plotorder': plotorder,
        'labels': labels,
        'lims': lims,
        'masks': masks,
        'figsize': (10,5),
        'dpi': 300,
        'depthlbl': 'in',
        'vollbl': 'TAF',
        'elevlbl': 'ft',
        'dplcs': 2,
        'barcolors': barcolors,
        'xlims': (0,len(edges)),
        'depth_clip': 0.01,
        'percent_min': 0.5,
        'percent_max': 99.5,
        'basins': basins,
        'title': 'SAMPLE'}

df = collect(connector, plotorder, basins, start_date, end_date, 'swe_vol',
             run_name, edges, 'end')

args['df'] = df
args['image'] = image[0,:,:]

swe_volume(args, logger = None)

```
