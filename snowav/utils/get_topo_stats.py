
import os
import netCDF4 as nc
from spatialnc import ipw

def get_topo_stats(fp, filetype='netcdf'):
    """
    Get stats about topo from the topo file
    Returns:
        ts - dictionary of topo header data
    """

    fp = os.path.abspath(fp)

    ts = {}

    if filetype == 'netcdf':
        ds = nc.Dataset(fp, 'r')
        # ts['units'] = ds.variables['y'].units
        ts['y'] = ds.variables['y'][:]
        ts['x'] = ds.variables['x'][:]
        ts['nx'] = len(ts['x'])
        ts['ny'] = len(ts['y'])
        ts['du'] = ts['y'][1] - ts['y'][0]
        ts['dv'] = ts['x'][1] - ts['x'][0]
        ts['v'] = ts['x'][0]
        ts['u'] = ts['y'][0]
        ds.close()

    if filetype == 'ipw':
        i = ipw.IPW(fp)
        ts['nx'] = i.nsamps
        ts['ny'] = i.nlines
        ts['units'] = i.bands[0].units
        ts['du'] = i.bands[0].dline
        ts['dv'] = i.bands[0].dsamp
        ts['v'] = float(i.bands[0].bsamp)
        ts['u'] = float(i.bands[0].bline)
        ts['csys'] = i.bands[0].coord_sys_ID

    return ts
