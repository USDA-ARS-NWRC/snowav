
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
import warnings

def point_values(output, stn_df, imgxy, filename):
    '''
    Writes to csv the model pixel values at selected locations.

    Inputs:
        output: array of output values (self.outputs['swe_z'][-1] from
            read_config.py, generated from iSnobalReader)
        stn_df: DataFrame, generated from csv file
            name        latitude    longitude
            Bishop Pass 37.1        -118.557
            Bench Lake  36.958      -118.445
        imgxy: array of x and y coordinates for the image, generated from
            get_topo_stats.py
        filename: name for csv file

    '''

    # dataframe and constants
    c = 1.0/25.4
    cols = ['x','y','xind','yind','swe [in]','median swe [in]']
    stns = stn_df['name']
    pixel_swe = pd.DataFrame(index=stns, columns=cols)

    # plt.close(50)
    # f = plt.figure(num = 50)
    # a = plt.gca()
    # a.imshow(output)

    # get closest model pixel
    for i,sta in enumerate(stns):
        lat = stn_df[stn_df['name'] == sta]['latitude'].values
        long = stn_df[stn_df['name'] == sta]['longitude'].values
        ll = utm.from_latlon(lat,long)
        xind = np.where(abs(imgxy[0]-ll[0]) == min(abs(imgxy[0]-ll[0])))[0]
        yind = np.where(abs(imgxy[1]-ll[1]) == min(abs(imgxy[1]-ll[1])))[0]

        if (xind == len(output[0,:])) or (yind == len(output[:,0])):
            warnings.warn('Latitude and/or longitude listed in snow course ' +
                          'csv are likely outside the model domain!')

        # Assign values
        pixel_swe.loc[sta,'x'] = float(stn_df[stn_df['name'] == sta]['latitude'].values)
        pixel_swe.loc[sta,'y'] = float(stn_df[stn_df['name'] == sta]['longitude'].values)
        pixel_swe.loc[sta,'xind'] = int(xind)
        pixel_swe.loc[sta,'yind'] = int(yind)

        # Get the median of the nine closest pixels
        val = np.array([])
        for n in range(-1,1):
            for m in range(-1,1):
                val = np.append(val,output[yind+n,xind+m]*c)
                if n == 0 and m == 0:
                    pixel_swe.loc[sta,'swe [in]'] = np.round(float(output[yind+n,xind+m]*c),decimals = 2)

        pixel_swe.loc[sta,'median swe [in]'] = np.round(np.nanmedian(val),decimals = 2)

        # a.plot(xind,yind,'gx',markersize = 10)

    pixel_swe.to_csv(filename)
    # plt.savefig('{}tmp_snowcourse_{}.png'.format(snow.figs_path,snow.name_append))
