
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
from snowav.utils.wyhr import calculate_date_from_wyhr
from spotpy import objectivefunctions

def stn_validate(snow):

    rundirs = snow.lrdirs
    stns = snow.val_stns
    lbls = snow.val_lbls
    client = snow.val_client
    factor = 25.4

    # get metadata from the data base from snotel sites
    if snow.basin == 'BRB':
        qry = ('SELECT tbl_metadata.* FROM tbl_metadata '
               + 'INNER JOIN tbl_stations ON tbl_metadata.primary_id = '
               + 'tbl_stations.station_id WHERE tbl_stations.client = '
               + ' "'"%s"'" HAVING network_name = "'"SNOTEL"'";'%client)
    else:
        qry = ('SELECT tbl_metadata.* FROM tbl_metadata '
               + 'INNER JOIN tbl_stations ON tbl_metadata.primary_id ='
               + 'tbl_stations.station_id ;')

    cnx = mysql.connector.connect(user='markrobertson',
                                  password='whatdystm?1',
                                  host='10.200.28.137',
                                  database='weather_db')

    meta_sno = pd.read_sql(qry, cnx)
    meta_sno.index = meta_sno['primary_id']
    meta_sno = meta_sno[~meta_sno.index.duplicated(keep='first')]

    swe_meas = pd.DataFrame(index = pd.date_range(datetime(snow.wy - 1,10,1),
                            snow.end_date,freq='D'),
                            columns = stns)

    stns_pixel = []
    for sta in stns:
        stns_pixel = stns_pixel + ['{}_'.format(sta) + str(s) for s in range(0,9)]

    swe_mod = pd.DataFrame(index = pd.date_range(datetime(snow.wy - 1,10,1),
                           snow.end_date,freq='D'),
                           columns = stns_pixel)

    tbl = 'tbl_level1'
    var = 'snow_water_equiv'
    st_time = '{}-10-1 00:00:00'.format(str(snow.wy - 1))
    end_time = snow.end_date.date().strftime("%Y-%-m-%-d")
    ncxvec = snow.snow_x
    ncyvec = snow.snow_y

    # station results
    for iters,stn in enumerate(stns):
        cnx = mysql.connector.connect(user='markrobertson',
                                      password='whatdystm?1',
                                      host='10.200.28.137',
                                      port='32768',database='weather_db')

        var_qry = ('SELECT weather_db.{0}.date_time, weather_db.{0}.{1} FROM '
                   'weather_db.{0} WHERE weather_db.{0}.date_time between "{2}" '
                   'and "{3}" AND weather_db.{0}.station_id IN '
                   ' ("{4}") ;'.format(tbl,var,st_time,end_time,stn))

        data = pd.read_sql(var_qry, cnx, index_col=bytes(bytearray(b'date_time')))
        data.index.names=['date_time']
        dind = pd.date_range(st_time,end_time,freq='D')
        swe_meas[stn] = data.reindex(dind)

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

    px = (1,1,1,0,0,0,-1,-1,-1)
    py = (1,0,-1,1,0,-1,1,0,-1)
  
    for rname in rundirs:

        d = rname.split('runs/run')[-1]
        folder_date = datetime(int(d[:4]),int(d[4:6]),int(d[6:8]))

        # Only load the rundirs that we need
        if (folder_date.date() <= snow.end_date.date()):

            ncf = nc.Dataset(os.path.join(rname,'snow.nc'), 'r')
            # nctvec = ncf.variables['time'][:]

            for iters,stn in enumerate(stns):
                ll = utm.from_latlon(meta_sno.ix[stn,'latitude'],meta_sno.ix[stn,'longitude'])
                xind = np.where(abs(ncxvec-ll[0]) == min(abs(ncxvec-ll[0])))[0]
                yind = np.where(abs(ncyvec-ll[1]) == min(abs(ncyvec-ll[1])))[0]

                for ix, (n,m) in enumerate(zip(px,py)):
                    land_stn = stn + '_{}'.format(str(ix))
                    swe = ncf.variables['specific_mass'][:,yind+m,xind+m][0][0][0]
                    swe_mod.loc[folder_date.date(),land_stn] = swe

                    # z = snow.dem[yind,xind]
                    if n == 0 and m == 0:
                        axs[iters].plot(swe_meas[stn]/factor,'k',label='measured')

                    axs[iters].plot(swe_mod[land_stn]/factor,'b',linewidth = 0.75,label='model')
                    axs[iters].set_title(lbls[iters])
                    axs[iters].set_xlim((datetime(snow.wy - 1, 10, 1),snow.end_date))

            ncf.close()

    # Nash-Sutcliffe on pixel 0,0
    if snow.nash_sut_flag:
        for iters,sta in enumerate(stns):
            pz_sta = sta + '_5'
            
            nsv = objectivefunctions.nashsutcliffe(swe_meas[swe_meas.index > datetime(snow.wy - 1, 10, 15)][stn].values.tolist(),
                                                            swe_mod[swe_mod.index > datetime(snow.wy - 1, 10, 15)][pz_sta].values.tolist())
            nstr = 'Nash-Sutcliffe:\n{}'.format(str(np.round(nsv,2)))

            nsx = 0.05
            nsy = 0.85

            if iters == 0:
               nsy = 0.675

            if nsv > 0:
                axs[iters].text(nsx,nsy,nstr,horizontalalignment='left',
                                transform=axs[iters].transAxes,fontsize = 8)

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

        for n in (1,4,7):
            axs[n].set_yticklabels('')            

    if (len(stns) > 9) and (len(stns) <=12):
        for n in (2,5,8,11):
            axs[n].yaxis.tick_right()

        for n in (9,10,11):
            for tick in axs[n].get_xticklabels():
                tick.set_rotation(30)
        for n in (0,1,2,3,4,5,6,7,8):
            axs[n].set_xticklabels('')

    # Plot
    swe_meas = swe_meas.replace('[]',np.nan)
    maxm = np.nanmax(swe_meas.max().values/factor)
    maxi = np.nanmax(swe_mod.max().values/factor)

    if maxm > maxi:
        maxswe = maxm
    else:
        maxswe = maxi

    for iters in range(0,len(stns)):
        axs[iters].set_ylim((-0.1,maxswe + maxswe*0.05))

    axs[0].legend(['measured','modeled'],loc='upper left')
    axs[0].set_ylabel('SWE [in]')

    plt.suptitle('Validation at Measured Sites')
    plt.subplots_adjust(top=0.92)

    snow._logger.info(' saving {}validation_{}.png'.format(snow.figs_path,snow.name_append))
    plt.savefig('{}validation_{}.png'.format(snow.figs_path,snow.name_append))
