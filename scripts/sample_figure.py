
from datetime import datetime
import netCDF4 as nc
from snowav.plotting.swe_volume import swe_volume
from snowav.plotting.plotlims import plotlims as plotlims
from snowav.database.database import connect, make_session, collect
from snowav.utils.utilities import masks

"""
See README.md and CoreConfig.ini for more options and details.

"""

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
        'figs_path': '<path>',
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
