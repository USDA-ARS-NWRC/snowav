
import numpy as np
import matplotlib
from matplotlib import pyplot as plt
from mpl_toolkits.axes_grid1 import make_axes_locatable
import seaborn as sns
import utm
import netCDF4 as nc
import pandas as pd
from datetime import datetime
import copy
from datetime import date
import os

def point_values(outputs, value, ncpath, x, y, savepath, append):
    '''
    Args
        outputs: dict of nc outputs from all the run_dirs (snow.outputs)
        value: value to plot ('swe_z', 'depth', etc)
        ncpath: path to one of the nc from run_dirs, used for x,y coordinates
        savepath: base figure save path (self.figs_path)
        append: name to append to file (self.name_append)

    '''

    # outputs['swe_z']

    ncf = nc.Dataset(ncpath, 'r')
    ncxvec = ncf.variables['x'][:]
    ncyvec = ncf.variables['y'][:]

    ll = utm.from_latlon(x,y)
    # ll      = utm.from_latlon(37.641922,-119.055443)
    # xind = np.where(abs(ncxvec-ll[0]) == min(abs(ncxvec-ll[0])))[0]
    # yind = np.where(abs(ncyvec-ll[1]) == min(abs(ncyvec-ll[1])))[0]
    xind = 1000
    yind = 500

    pt = np.array([])
    for n in range(0,len(outputs['swe_z'][:])):
        np.append(pt,outputs['swe_z'][n][xind,yind])

    fig = plt.figure(num=15, figsize=(8,6))
    ax = plt.gca()

    ax.plot(pt)

    # plt.show()

    # snow._logger.info('saving figure to {}point_validation_{}.png'.format(savepath,append))
    plt.savefig('{}point_validation_{}.png'.format(savepath,append))
