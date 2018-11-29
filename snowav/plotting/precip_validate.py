
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

def precip_validate(snow):

    # if snow.valid_flag == False:
    #     snow._logger.debug('No stations listed in config file for validation figure!')
    #     return

    rundirs = snow.run_dirs
    stns = snow.pre_val_stns
    lbls = snow.pre_val_lbls
    client = snow.val_client

    # get metadata from the data base from snotel sites
    qry = ('SELECT tbl_metadata.* FROM tbl_metadata '
           + 'INNER JOIN tbl_stations ON tbl_metadata.primary_id ='
           + 'tbl_stations.station_id WHERE tbl_stations.client = '
           + '"'"%s"'" ;'%client)

    cnx = mysql.connector.connect(user='markrobertson',
                                  password='whatdystm?1',
                                  host='10.200.28.137',
                                  database='weather_db')

    meta_sno = pd.read_sql(qry, cnx)
    meta_sno.index = meta_sno['primary_id']
    pre_meas = pd.DataFrame(index = pd.date_range(datetime(snow.wy - 1,10,1),
                                                     snow.end_date,
                                                     freq='D'),columns = stns)
    pre_mod = pd.DataFrame(index = pd.date_range(datetime(snow.wy - 1,10,1),
                                                     snow.end_date,
                                                     freq='D'),columns = stns)

    tbl = 'tbl_level1'
    var = 'precip_accum'
    st_time = '%s-10-1 00:00:00'%(str(snow.wy - 1))
    end_time = snow.end_date.date().strftime("%Y-%-m-%-d")

    # Get Snotel station results
    for iters,stn in enumerate(stns):
        cnx = mysql.connector.connect(user='markrobertson',
                                      password='whatdystm?1',
                                      host='10.200.28.137',
                                      port='32768',database='weather_db')

        var_qry = ('SELECT weather_db.%s.date_time, weather_db.%s.%s ' % (tbl,tbl,var) +
                    'FROM weather_db.%s ' % tbl +
                    "WHERE weather_db.%s.date_time between '" % tbl + st_time+ "' and '"+end_time+"'"
                    "AND weather_db.%s.station_id IN ('" % tbl + stn + "');")

        try:
            data = pd.read_sql(var_qry, cnx, index_col=bytes(bytearray(b'date_time')))
        except:
            data = pd.read_sql(var_qry, cnx, index_col='date_time')

        data.index.names=['date_time']
        dind = pd.date_range(st_time,end_time,freq='D')
        pre_meas[stn] = data.reindex(dind)

    sns.set_style('darkgrid')
    sns.set_context('notebook')
    plt.close(9)

    if len(stns) <= 6:
        fig, axs = plt.subplots(num=9, figsize=(10,10), nrows=3, ncols=2)

    if (len(stns) > 6) and (len(stns) <= 9):
        fig, axs = plt.subplots(num=9, figsize=(10,10), nrows=3, ncols=3)

    if (len(stns) > 9) and (len(stns) <= 12):
        fig, axs = plt.subplots(num=9, figsize=(10,10), nrows=4, ncols=3)

    axs = axs.flatten()

    if ( (snow.offset == 0)
        and (snow.outputs['dates'][0].date() != datetime(snow.wy-1,10,1)) ):
        iswe = (snow.outputs['dates'][0].date() -  datetime(snow.wy-1,10,1).date()).days
    else:
        iswe = snow.offset

    for rname in rundirs:

        ncpath = rname.split('output')[0]
        path = ncpath.replace('runs', 'data')
        path = path.replace('run','data')

        ncf = nc.Dataset(os.path.abspath(path + '/smrfOutputs/precip.nc'), 'r')
        nctvec = ncf.variables['time'][:]

        precip = np.zeros((snow.nrows,snow.ncols))
        for s in range(0,len(nctvec)):
            precip = precip + ncf.variables['precip'][s]

        ncxvec = ncf.variables['x'][:]
        ncyvec = ncf.variables['y'][:]

        for stn in stns:
            ll = utm.from_latlon(meta_sno.ix[stn,'latitude'],meta_sno.ix[stn,'longitude'])
            xind = np.where(abs(ncxvec-ll[0]) == min(abs(ncxvec-ll[0])))[0]
            yind = np.where(abs(ncyvec-ll[1]) == min(abs(ncyvec-ll[1])))[0]
            pre = pd.Series(precip[yind,xind])

            try:
                pre_mod.loc[iswe:(iswe + len(pre.values)),stn] = pre.values

            except:
                sv = pre_mod[stn].values
                lx = len(sv[iswe::])
                pre_mod.loc[iswe:(iswe + lx),stn] = pre.values[0:lx]

        ncf.close()
        iswe = iswe + len(pre.values)

        z = snow.dem[yind,xind]

    for iters, stn in enumerate(stns):
        axs[iters].plot(pre_meas[stn],'k',label='measured')
        axs[iters].plot(pre_mod[stn].cumsum(skipna=True),'b',linewidth = 0.75,label='model')
        axs[iters].set_title(lbls[iters])
        axs[iters].set_xlim((datetime(snow.wy - 1, 10, 1),snow.end_date))

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
    pre_meas = pre_meas.replace('[]',np.nan)
    maxm = np.nanmax(np.nanmax(pre_mod.cumsum(skipna=True).values))

    for iters in range(0,len(stns)):
        axs[iters].set_ylim((-0.1,maxm + maxm*0.3))

    axs[0].legend(['measured','HRRR'],loc='upper left')
    axs[0].set_ylabel('precip [mm]')

    plt.suptitle('Precipitation Validation')
    plt.subplots_adjust(top=0.92)

    snow._logger.info('saving figure to %sprecip_validation_%s.png'%(snow.figs_path,snow.name_append))
    plt.savefig('%sprecip_validation_%s.png'%(snow.figs_path,snow.name_append))
