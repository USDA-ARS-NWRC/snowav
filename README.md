# SNOw and Water model Analysis and Visualization

[![GitHub version](https://badge.fury.io/gh/USDA-ARS-NWRC%2Fsnowav.svg)](https://badge.fury.io/gh/USDA-ARS-NWRC%2Fsnowav)

SNOw and Water model Analysis and Visualization was developed at the USDA
Agricultural Research Service in Boise, Idaho. It processes model
outputs from AWSM into formats and figures for
water resource managers. See the CoreConfig.ini in this repo for more information on config options.

![image](https://raw.githubusercontent.com/USDA-ARS-NWRC/awsm/master/docs/_static/ModelSystemOverview_new.png)

To do:
- non-awsm_daily
- supply dates only rather than 23:00 

Currently snowav requires:
- awsm model results in awsm_daily format, including  output files in the paths ```.../runs/runYYYYMMDD/snow.nc``` and ```.../runs/runYYYYMMDD/em.nc```
- topo.nc files that have been created by the basin_setup package


snowav expects:

- input files in the paths ```.../data/dataYYYYMMDD/smrfOutputs/precip.nc``` and ```.../data/dataYYYYMMDD/smrfOutputs/percent_snow.nc```. If these do not exist the precip figures will be turned off.
- the lidar_depths_wy2019.nc file, specified in the config file in [plots] update_file, contains all relevant flights that have been flown. Flights in this file can be removed using [plots] update_numbers, but they can't be added.

The standard manual snowav processing run is:
```
snowav -f config.ini
```

Notes and considerations:
- if depth updates have been applied, year-to-date precipitation values are no longer valid


If results have already been processed and put onto a database, figures
can be created outside of a snowav processing run.

```
import netcdf as nc
from snowav.plotting.swe_volume import swe_volume
from snowav.plotting.plotlims import plotlims as plotlims
from snowav.database.database import connect, make_session
from snowav.utils.utilities import masks

""" See CoreConfig.ini for more options and details """

connector = '< path_to_database >'
dempath = '< path_to_topo.nc >'
plotorder = ['Extended Tuolumne', 'Tuolumne', 'Cherry Creek', 'Eleanor']
lims = plotlims(plotorder)
edges = [3000 4000 5000 6000 7000 8000 9000 10000 11000 12000 13000]

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

barcolors = ['xkcd:cobalt', 'xkcd:mustard green', 'xkcd:lichen',
             'xkcd:pale green', 'xkcd:blue green', 'xkcd:bluish purple',
             'xkcd:lightish purple', 'xkcd:deep magenta']

basins, cnx, out = connect(sqlite=connector, mysql=None, plotorder=plotorder,
                           user=None, password=None, host=None, port=None)             

# suggested args with: [data type] (suggested)
args = {'report_start': [datestr],
        'report_date': [datestr],
        'print': [bool] (False),
        'run_name': [str] ('tuol_wy2019_ops'),
        'start_date': [datetime],
        'end_date': [datetime],
        'directory': [str] ('test')
        'figs_path': [str] ('/<path>/'),
        'edges': [np.array], (edges)
        'plotorder': [list] (plotorder),
        'labels': [dict] (labels),
        'lims': [object] (lims),
        'masks': [dict] (masks),
        'figsize': [array] (10,5),
        'dpi': [int] (250),
        'depthlbl': [str] ('in'),
        'vollbl': [str] ('TAF'),
        'elevlbl': [str] ('ft'),
        'dplcs': [int] (2),
        'barcolors': [list], (barcolors)
        'xlims': [array] (0,len(edges)),
        'depth_clip': [float] (0.01),
        'percent_min': [float] (0.5),
        'percent_max': [float] (99.5),
        'basins': [dict] (basins)}

df = collect(connector, plotorder, start_date, end_date, 'swe_vol',
             run_name, edges, 'end')

args['df'] = df
args['image'] = image
args['title'] = title

swe_volume(args, logger = None)

```
