
from snowav.utils.utilities import calculate
import netCDF4 as nc

"""
Sample processing script for snowav processing utility. Calculate values from
iSnobal output snow.nc and em.nc files and convert to desired units, including:
    mean SWE depth
    total SWE volume
    mean snow depth
    total SWI volume
    mean SWI depth
    mean evaporation
    total cold content
    mean density

Applies all masks supplied, pass multiple masks as a list.

"""

# Specify the iSnobal output file and topo.nc for the basin mask (optional).
file_path = '/data/blizzard/tuolumne/ops/wy2019/ops/runs/run20190101/snow.nc'
topo_path = '/home/ops/wy2019/tuolumne/topo/topo.nc'
value = 'specific_mass'
pixel = 50

# Get desired output value
ncf = nc.Dataset(file_path)
array = ncf.variables[value][:]
ncf.close()

# Optional - get total basin mask from the topo.nc file.
# Can also pass additional masks, grouped in a list
ncf = nc.Dataset(topo_path)
mask = ncf['mask'][:]
ncf.close()

# SWE volume [TAF]
swe_vol = calculate(array[0], pixel, mask, 'sum', 'volume', 'TAF')

# mean SWE depth [in]
swe_in = calculate(array[0], pixel, mask, 'mean', 'depth', 'TAF')
