
import os
import netCDF4 as nc

def get_topo_stats(path):
    """
    Get stats about topo from the topo file.

    Returns
    ------
    data - dictionary of topo header data
    """

    path = os.path.abspath(path)
    data = nc.Dataset(path, 'r')

    stats = {}
    stats['y'] = data.variables['y'][:]
    stats['x'] = data.variables['x'][:]
    stats['nx'] = len(stats['x'])
    stats['ny'] = len(stats['y'])
    stats['du'] = stats['y'][1] - stats['y'][0]
    stats['dv'] = stats['x'][1] - stats['x'][0]
    stats['v'] = stats['x'][0]
    stats['u'] = stats['y'][0]

    data.close()

    return stats
