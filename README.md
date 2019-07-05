# SNOw and Water model Analysis and Visualization

[![GitHub version](https://badge.fury.io/gh/USDA-ARS-NWRC%2Fsnowav.svg)](https://badge.fury.io/gh/USDA-ARS-NWRC%2Fsnowav)

SNOw and Water model Analysis and Visualization was developed at the USDA Agricultural Research Service in Boise, Idaho. It processes model
outputs from AWSM into formats and figures for water resource managers. See the CoreConfig.ini in this repo for more information on config options.

![image](https://raw.githubusercontent.com/USDA-ARS-NWRC/awsm/master/docs/_static/ModelSystemOverview_new.png)

## Install

## Requirements
Currently snowav requires:
- awsm model results in awsm_daily format, including output files in the paths ```.../runs/runYYYYMMDD/snow.nc``` and ```.../runs/runYYYYMMDD/em.nc```
- topo.nc files that have been created by the basin_setup package
- correct date information in all snow.nc files

snowav expects:

- input files in the paths ```.../data/dataYYYYMMDD/smrfOutputs/precip.nc``` and ```.../data/dataYYYYMMDD/smrfOutputs/percent_snow.nc```. If these do not exist the precip figures will be turned off.
- the lidar_depths_wy2019.nc file, specified in the config file in [plots] update_file, contains all relevant flights that have been flown. Flights in this file can be removed from processing using [plots] update_numbers, but they can't be added.

## Usage
The standard snowav model results database processing and figures run is:
```
$ snowav -f config.ini
```

To process a single snow.nc file and display simple SWE volume figure, without putting results on a database:

```
$ snowav -t <topo.nc> -d <snow.nc>
```

By default this will display SWE volume and save the figure in the snowav repo. Specify a figure save location with ```-p <path>``` and include ``` -s ``` to display the figure.

For a simple difference between two snow.nc files:

```
$ snowav -t <topo.nc> -A <snow.nc> -B <snow.nc>
```

## Notes
- if depth updates have been applied, year-to-date precipitation values are no longer valid

## Figures from existing database
If results have already been processed and put onto a database, figures can be created outside of a snowav processing run (see also script/sample_figure.py). See snowav.framework.figures for templates for additional figure creation.

```
from datetime import datetime
import netCDF4 as nc
from snowav.plotting.swe_volume import swe_volume
from snowav.plotting.plotlims import plotlims as plotlims
from snowav.database.database import connect, make_session, collect
from snowav.utils.utilities import masks

""" See CoreConfig.ini for more options and details """

connector = 'sqlite:////<path>.db'
dempath = '/<path>/topo.nc'
plotorder = ['Extended Tuolumne', 'Tuolumne', 'Cherry Creek', 'Eleanor']
lims = plotlims(plotorder)
edges = [3000, 4000, 5000, 6000, 7000, 8000, 9000, 10000, 11000, 12000, 13000]

ncf = nc.Dataset('/<path>/snow.nc')
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

## To do:
- non-awsm_daily
- supply dates only rather than 23:00
