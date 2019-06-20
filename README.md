# SNOw and Water model Analysis and Visualization

[![GitHub version](https://badge.fury.io/gh/USDA-ARS-NWRC%2Fsnowav.svg)](https://badge.fury.io/gh/USDA-ARS-NWRC%2Fsnowav)

SNOw and Water model Analysis and Visualization was developed at the USDA
Agricultural Research Service (ARS) in Boise, ID. SNOWAV processes model
outputs from AWSM into basin and sub-basins formats that are meaningful for
water resource managers and scientists.

![image](https://raw.githubusercontent.com/USDA-ARS-NWRC/awsm/master/docs/_static/ModelSystemOverview_new.png)

The standard manual snowav processing run is:
```
snowav -f config.ini
```

If results have already been processed and put onto a database, figures
can be created outside of a snowav processing run.

```
import netcdf as nc
from snowav.plotting.swe_volume import swe_volume
from snowav.plotting.plotlims import plotlims as plotlims

""" See CoreConfig.ini for more options and details """

ncf = nc.Dataset('/<path>/snow.nc')
image = ncf['specific_mass'][:]
ncf.close()

connector = '< path_to_database >'
dempath = '< path_to_topo.nc'
plotorder = ['Extended Tuolumne', 'Tuolumne', 'Cherry Creek', 'Eleanor']
lims = plotlims(plotorder)
edges = [3000 4000 5000 6000 7000 8000 9000 10000 11000 12000 13000]

out = snowav_masks(dempath)
dem = out['dem']
masks = out['masks']
nrows = out['nrows']
ncols = out['ncols']
plotorder = out['plotorder']
labels = out['labels']

barcolors = ['xkcd:cobalt', 'xkcd:mustard green', 'xkcd:lichen',
             'xkcd:pale green', 'xkcd:blue green', 'xkcd:bluish purple',
             'xkcd:lightish purple', 'xkcd:deep magenta']

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
        'figsize': [array] (10,8),
        'dpi': [int] (250),
        'depthlbl': [str] ('in'),
        'vollbl': [str] ('TAF'),
        'elevlbl': [str] ('ft'),
        'dplcs': [int] (2),
        'barcolors': [list], (barcolors)
        'xlims': [array] (0,len(edges)),
        'depth_clip': [float] (0.01),
        'percent_min': [float] (0.5),
        'percent_max': [float] (99.5)}

df = collect(connector, plotorder, start_date, end_date, 'swe_vol',
             run_name, edges, 'end')

args['df'] = df
args['image'] = image
args['title'] = title

swe_volume(args, logger = None)

```
