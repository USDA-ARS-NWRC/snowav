import netCDF4 as nc
import numpy as np
from matplotlib import pyplot as plt

file = '/data/snowdrift/kaweah/ops/wy2019/ops/runs/run20181218/snow.nc'
out = '/home/markrobertson/wkspace/kaweah_20181218_swe.asc'

snow = nc.Dataset(file,'r')

swe = snow['specific_mass'][0]
np.savetxt(out,swe,delimiter ='\t')

swe = np.loadtxt(out)

f = plt.figure()
ax = plt.gca()

ax.imshow(swe)
plt.show()
