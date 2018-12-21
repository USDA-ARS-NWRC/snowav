
import numpy as np
import matplotlib
from matplotlib import pyplot as plt
from mpl_toolkits.axes_grid1 import make_axes_locatable
import seaborn as sns
import utm
import netCDF4 as nc
import mysql.connector
import pandas as pd
from datetime import datetime
import copy
from datetime import date
import os

def air_temp(rundirs,dem,stns,lbls,client,wy,end_date,figs_path,name_append,start_date=None):
    '''
    '''

    # rundirs = snow.run_dirs
    # stns = snow.val_stns
    # lbls = snow.val_lbls
    # client = snow.val_client

    if start_date is None:
        start_date = datetime(wy - 1,10,1)

    # get metadata from the data base from snotel sites
    qry = ('SELECT tbl_metadata.* FROM tbl_metadata '
           + 'INNER JOIN tbl_stations ON tbl_metadata.primary_id ='
           + 'tbl_stations.station_id WHERE tbl_stations.client = '
           + '"'"{}"'" ;'.format(client))

    cnx = mysql.connector.connect(user='markrobertson',
                                  password='whatdystm?1',
                                  host='10.200.28.137',
                                  database='weather_db')

    meta_sno = pd.read_sql(qry, cnx)
    meta_sno.index = meta_sno['primary_id']
    meas = pd.DataFrame(index = pd.date_range(start_date,
                                              end_date,
                                              freq='H'),
                                              columns = stns)
    mod = pd.DataFrame(index = pd.date_range(start_date,
                                              end_date,
                                              freq='H'),
                                              columns = stns)

    tbl = 'tbl_level1'
    var = 'air_temp'
    st_time = start_date.date().strftime("%Y-%-m-%-d")
    end_time = end_date.date().strftime("%Y-%-m-%-d")

    # Get Snotel station results
    for iters,stn in enumerate(stns):
        cnx = mysql.connector.connect(user='markrobertson',
                                      password='whatdystm?1',
                                      host='10.200.28.137',
                                      port='32768',database='weather_db')

        # var_qry = ('SELECT weather_db.{0}.date_time, weather_db.{0}.{1} FROM ' +
        #           'weather_db.{0} WHERE weather_db.{0}.date_time between '+
        #           ''"{2} "' and '"{3}"' AND weather_db.{0}.station_id IN '+
        #           ''"'{4}'"';').format(tbl,var,st_time,end_time,stn)

        var_qry = ('SELECT weather_db.%s.date_time, weather_db.%s.%s ' % (tbl,tbl,var) +
                    'FROM weather_db.%s ' % tbl +
                    "WHERE weather_db.%s.date_time between '" % tbl + st_time+ "' and '"+end_time+"'"
                    "AND weather_db.%s.station_id IN ('" % tbl + stn + "');")

        try:
            data = pd.read_sql(var_qry, cnx, index_col=bytes(bytearray(b'date_time')))
        except:
            data = pd.read_sql(var_qry, cnx, index_col='date_time')

        data.index.names=['date_time']
        dind = pd.date_range(st_time,end_time,freq='H')
        meas[stn] = data.reindex(dind)

    sns.set_style('darkgrid')
    sns.set_context('notebook')

    plt.close(30)

    if len(stns) <= 6:
        fig, axs = plt.subplots(num=30, figsize=(10,10), nrows=3, ncols=2)
        fig1, axs1 = plt.subplots(num=31, figsize=(10,10), nrows=3, ncols=2)

    if (len(stns) > 6) and (len(stns) <= 9):
        fig, axs = plt.subplots(num=30, figsize=(10,10), nrows=3, ncols=3)
        fig1, axs1 = plt.subplots(num=31, figsize=(10,10), nrows=3, ncols=3)

    if (len(stns) > 9) and (len(stns) <= 12):
        fig, axs = plt.subplots(num=30, figsize=(10,10), nrows=4, ncols=3)
        fig1, axs1 = plt.subplots(num=31, figsize=(10,10), nrows=4, ncols=3)

    axs = axs.flatten()
    axs1 = axs1.flatten()
    z = {}

    for iters,stn in enumerate(stns):

        iswe = 0

        for rname in rundirs:
            sf = rname.replace('runs','data')
            sf = sf.replace('run','data')
            ncpath = sf.split('output')[0] + '/smrfOutputs/air_temp.nc'
            ncf = nc.Dataset(ncpath, 'r')
            nctvec = ncf.variables['time'][0:23]
            val = ncf.variables['air_temp']

            ncxvec = ncf.variables['x'][:]
            ncyvec = ncf.variables['y'][:]
            ll = utm.from_latlon(meta_sno.ix[stn,'latitude'],meta_sno.ix[stn,'longitude'])
            # ll      = utm.from_latlon(37.641922,-119.055443)

            xind = np.where(abs(ncxvec-ll[0]) == min(abs(ncxvec-ll[0])))[0]
            yind = np.where(abs(ncyvec-ll[1]) == min(abs(ncyvec-ll[1])))[0]
            z[stn] = dem[yind,xind]
            swe = pd.Series(val[0:23,yind,xind].flatten(),index=nctvec)

            try:
                mod.loc[iswe:(iswe + len(swe.values)),stn] = swe.values

            except:
                sv = mod[stn].values
                lx = len(sv[iswe::])
                mod.loc[iswe:(iswe + lx),stn] = swe.values[0:lx]

            ncf.close()
            iswe = iswe + len(swe.values)

        # z = snow.dem[yind,xind]
        axs[iters].plot(meas[stn],'k',label='station')
        axs[iters].plot(mod[stn],'b',linewidth = 0.75,label='model')
        axs[iters].set_title(lbls[iters])
        # axs[iters].set_xlim((datetime(wy - 1, 10, 1),end_date))

        mask = pd.notna(meas[stn])
        axs1[iters].plot(meas[stn][mask].values,mod[stn][mask].values,'ko',markersize=5)

    if len(stns) <= 6:
        for n in (1,3,5):
            axs[n].yaxis.tick_right()

        for n in (4,5):
            for tick in axs[n].get_xticklabels():
                tick.set_rotation(30)
        for n in (0,1,2,3):
            axs[n].set_xticklabels('')

    if (len(stns) > 6) and (len(stns) <=9):
        for n in (2,5,8):
            axs[n].yaxis.tick_right()

        for n in (6,7,8):
            for tick in axs[n].get_xticklabels():
                tick.set_rotation(30)
        for n in (0,1,2,3,4,5):
            axs[n].set_xticklabels('')

    if (len(stns) > 9) and (len(stns) <=12):
        for n in (2,5,8,11):
            axs[n].yaxis.tick_right()

        for n in (9,10,11):
            for tick in axs[n].get_xticklabels():
                tick.set_rotation(30)
        for n in (0,1,2,3,4,5,6,7,8):
            axs[n].set_xticklabels('')

    # Plot
    meas = meas.replace('[]',np.nan)

    for iters in range(0,len(stns)):
        axs[iters].set_xlim((start_date,end_date))

    axs[0].legend(loc='upper left')
    axs[0].set_ylabel('SWE [mm]')

    plt.suptitle('Validation at Measured Sites')
    plt.subplots_adjust(top=0.92)

    # snow._logger.info('saving figure to {}air_temp_{}.png'.format(snow.figs_path,snow.name_append))
    plt.savefig('{}air_temp_{}.png'.format(figs_path,name_append))
    # plt.show()

    # mask = pd.notna(meas[stn])
    # plt.savefig('{}air_temp_test{}.png'.format(figs_path,name_append))
