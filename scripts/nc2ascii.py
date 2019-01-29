import netCDF4 as nc
import numpy as np
from matplotlib import pyplot as plt

basin = 'kaweah'

file = '/data/snowdrift/{}/ops/wy2019/ops/runs/run20181220/snow.nc'.format(basin)
out = '/home/markrobertson/wkspace/{}_20181220_swe.asc'.format(basin)

snow = nc.Dataset(file,'r')

swe = snow['specific_mass'][0]
np.savetxt(out,swe,delimiter ='\t')

swe = np.loadtxt(out)

f = plt.figure()
ax = plt.gca()

ax.imshow(swe)
plt.show()
