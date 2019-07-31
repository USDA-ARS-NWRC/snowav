# SNOw and Water model Analysis and Visualization

[![GitHub version](https://badge.fury.io/gh/USDA-ARS-NWRC%2Fsnowav.svg)](https://badge.fury.io/gh/USDA-ARS-NWRC%2Fsnowav)

SNOWAV was developed at the USDA Agricultural Research Service in Boise, Idaho. It processes AWSM model outputs into formats and figures for water resource managers. See CoreConfig.ini for details on config options.

![image](https://raw.githubusercontent.com/USDA-ARS-NWRC/awsm/master/docs/_static/ModelSystemOverview_new.png)

## Requirements
Currently snowav requires:
- awsm model results in awsm_daily format, including output files in the paths ```.../runs/runYYYYMMDD/snow.nc``` and ```.../runs/runYYYYMMDD/em.nc```
- topo.nc files that have been created by the basin_setup package
- correct date information in all snow.nc files

snowav expects:

- input files in the paths ```.../data/dataYYYYMMDD/smrfOutputs/precip.nc``` and ```.../data/dataYYYYMMDD/smrfOutputs/percent_snow.nc```. If these do not exist the precip figures will be turned off.
- the lidar_depths_wy2019.nc file, specified in the config file in [plots] update_file, contains all relevant flights that have been flown. Flights in this file can be removed from processing using [plots] *update_numbers*, but they can't be added.

## Usage
The standard snowav model processing, database placement, and figures run is:
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

For awsm >= 0.9.26, if awsm config field [awsm master] snowav_config is given a valid snowav config file, a full snowav run with that config will be run at the end of all awsm functionality.

To query and output existing database records, either to the terminal or csv, see CoreConfig.ini [query] section and example below. If [query] *query: True* no other functionality other than the query will be run. This currently requires database login created from [database] section and *mysql: snowav*.
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
- **Be deliberate with the [snowav] *run_name* field!** This is how processing runs are identified on the database, and existing records with the same *run_name*, basin, and date range will be deleted and replaced during re-runs of the same snow.nc files. The *run_name* field should, in most cases and normal use, be connected with a specific [run] *directory*. We suggest using fields such as *tuol_wy2019_ops* and re-running snowav with *tuol_wy2019_ops* over the proper directory if modifications are made mid-season.

- Currently the figures listed below are created by default in both the [plots] and [reports] section. If any of them are set to *False* in [plots] they will be set to *False* in [reports] internally. For additional figures to be added to a report they must be *True* in both [plots] and [reports].
  - swi
  - image_change
  - cold_content  
  - swe_volume
  - basin_total
  - precip_depth


- If depth updates have been applied, year-to-date precipitation values in report Table 1 are no longer valid, although they will still appear. This table includes precipitation and rain derived from the HRRR inputs and does not account for mass that may be added or removed via snow depth updates.

- To help debugging the pdf report if it fails rendering to latex, set config option [report] *print_latex: True*. This will print the full latex file to the screen immediately prior to rendering.

- The SWE pillow validation figure [plots] *stn_validate* requires SWE fields on the existing weather database and [validate] fields to be filled in.

- Topo.nc, masks, and plot labels <br/>
Config field [snowav] *masks* can be left blank, and will default to the long_name fields in the topo.nc file. To subset the number of basins processed and plotted, use *masks* with a list: <br/> *masks: San Joaquin River Basin, Main, South Fork* <br/>
To replace the plot labels, use the *plotlabels* field in combination with *masks*: <br/> *masks: San Joaquin River Basin, Main, South Fork* <br/> *plotlabels: San Joaquin, Mammoth, South Fork*

Setting config option [inflow] *inflow: True* triggers reading in operator-generated inflow and the inflow figure. Currently only configured to work with the Tuolumne, and 'FORM11' excel sheets.

## SNOWAV processing utility
The snowav processing utility snowav.utils.utilities.calculate can be used to make simple calculations on snow.nc and em.nc files. See scripts.sample_process.py for an example.

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
