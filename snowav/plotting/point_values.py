
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

def point_values(output, imagexy, pointxy, savepath, filename):
    '''

    Inputs:
        output: array of output values (self.outputs from read_config.py,
            generated from iSnobalReader)
        imagexy: array of x and y coordinates for the image
        pointxy: array of x and y points to extract
        savepath: path to save csv file
        filename: name for csv file

    '''


    ncpath = rname.split('output')[0]
    ncf = nc.Dataset(os.path.join(ncpath,'snow.nc'), 'r')
    nctvec = ncf.variables['time'][:]
    vswe = ncf.variables['specific_mass']
    ncxvec = ncf.variables['x'][:]
    ncyvec = ncf.variables['y'][:]
    ll = utm.from_latlon(meta_sno.ix[stn,'latitude'],meta_sno.ix[stn,'longitude'])
    # ll      = utm.from_latlon(37.641922,-119.055443)

    xind = np.where(abs(ncxvec-ll[0]) == min(abs(ncxvec-ll[0])))[0]
    yind = np.where(abs(ncyvec-ll[1]) == min(abs(ncyvec-ll[1])))[0]
    swe = pd.Series(vswe[:,yind+m,xind+n].flatten(),index=nctvec)
