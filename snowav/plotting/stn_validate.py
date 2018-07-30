
import numpy as np
from matplotlib import pyplot as plt
from mpl_toolkits.axes_grid1 import make_axes_locatable
import seaborn as sns
import utm
import netCDF4 as nc
import mysql.connector
import pandas as pd
from datetime import datetime
import copy

def stn_validate(snow):

    if snow.valid_flag == False:
        snow._logger.debug('No stations listed in config file for validation figure!')
        return

    rundirs = snow.run_dirs
    stns = snow.val_stns
    lbls = snow.val_lbls
    client = snow.val_client

    # get metadata from the data base from snotel sites
    if snow.basin == 'BRB':
        qry = ('SELECT tbl_metadata.* FROM tbl_metadata '
               + 'INNER JOIN tbl_stations ON tbl_metadata.primary_id = '
               + 'tbl_stations.station_id WHERE tbl_stations.client = '
               + ' "'"%s"'" HAVING network_name = "'"SNOTEL"'";'%client)
    else:
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
    swe_meas = pd.DataFrame(index = pd.date_range(datetime(snow.wy - 1,10,1),
                                                     snow.end_date,
                                                     freq='D'),columns = stns)
    swe_mod = pd.DataFrame(index = pd.date_range(datetime(snow.wy - 1,10,1),
                                                     snow.end_date,
                                                     freq='D'),columns = stns)

    tbl = 'tbl_level1'
    var = 'snow_water_equiv'
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

        # data = pd.read_sql(var_qry, cnx, index_col=bytes(bytearray(b'date_time')))
        data = pd.read_sql(var_qry, cnx)
        data = data.set_index('date_time')
        dind = pd.date_range(st_time,end_time,freq='D')
        swe_meas[stn] = data.reindex(dind)

    if snow.wy == 2017:
        if 'TIOC1' in swe_meas:
            swe_meas.TIOC1 = swe_meas.TIOC1 - 300
        if 'AGP' in swe_meas:
            swe_meas.AGP = swe_meas.AGP - 40
        if 'VLC' in swe_meas:
            swe_meas.VLC = swe_meas.VLC + 250
        if 'UBC' in swe_meas:
            swe_meas.UBC = swe_meas.UBC - 50

    sns.set_style('darkgrid')
    sns.set_context('notebook')

    plt.close(9)
    fig, axs = plt.subplots(num = 9,figsize = (10,10),nrows = 3,ncols = 2)
    axs = axs.flatten()

    # First need to combine all nc files...
    px = (1,1,1,0,0,0,-1,-1,-1)
    py = (1,0,-1,1,0,-1,1,0,-1)
    # stnpix = []
    # for n,m in zip(px,py):
    #     stnpix = stnpix + ['%s,%s'%(str(m),str(n))]

    for iters,stn in enumerate(stns):
        # iters = 0
        # stn = stns[iters]
        # swe_mod_pix = pd.DataFrame(index = pd.date_range(datetime(snow.wy - 1,10,1),
        #                                 snow.end_date, freq='D'), columns = stnpix)

        for n,m in zip(px,py):
            # n = 0
            # m = 0
            iswe = snow.offset

            for rname in rundirs:
                # rname = rundirs[1]

                ncpath = rname.split('output')[0]
                ncf = nc.Dataset(ncpath + 'snow.nc', 'r')
                nctvec = ncf.variables['time'][:]
                vswe = ncf.variables['specific_mass']
                ncxvec = ncf.variables['x'][:]
                ncyvec = ncf.variables['y'][:]
                ll = utm.from_latlon(meta_sno.ix[stn,'latitude'],meta_sno.ix[stn,'longitude'])
                # ll      = utm.from_latlon(37.641922,-119.055443)

                xind = np.where(abs(ncxvec-ll[0]) == min(abs(ncxvec-ll[0])))[0]
                yind = np.where(abs(ncyvec-ll[1]) == min(abs(ncyvec-ll[1])))[0]
                swe = pd.Series(vswe[:,yind+m,xind+n].flatten(),index=nctvec)

                try:
                    swe_mod.loc[iswe:(iswe + len(swe.values)),stn] = swe.values
                    # swe_mod_pix.loc[iswe:(iswe + len(swe.values)),'%s,%s'%(n,m)] = swe.values
                except:
                    sv = swe_mod[stn].values
                    lx = len(sv[iswe::])
                    swe_mod.loc[iswe:(iswe + lx),stn] = swe.values[0:lx]
                    # swe_mod_pix.loc[iswe:(iswe + lx),'%s,%s'%(n,m)] = swe.values[0:lx]

                ncf.close()
                iswe = iswe + len(swe.values)

            z = snow.dem[yind,xind]
            axs[iters].plot(swe_meas[stn],'k',label='measured')
            axs[iters].plot(swe_mod[stn],'b',linewidth = 0.75,label='model')
            axs[iters].set_title(lbls[iters])
            axs[iters].set_xlim((datetime(snow.wy - 1, 10, 1),snow.end_date))

        if iters == 1 or iters == 3 or iters == 5:
            axs[iters].yaxis.tick_right()

        if iters == 4 or iters == 5:
            for tick in axs[iters].get_xticklabels():
                tick.set_rotation(30)
        else:
            axs[iters].set_xticklabels('')

        # swe_mod_pix = swe_mod_pix.astype(float).round(decimals=1)
        # swe_mod_pix.to_csv('/home/markrobertson/mrworkspace/projects/ferix/tuol_wy2018_swe_mm_%s.csv'%(stn))
        # print('/home/markrobertson/mrworkspace/projects/ferix/tuol_wy2017_swe_mm_%s.csv'%(stn))

    # Plot
    maxm = np.nanmax(swe_meas)
    maxi = np.nanmax(swe_mod.max().values)

    if maxm > maxi:
        maxswe = maxm
    else:
        maxswe = maxi

    for iters in range(0,len(stns)):
        axs[iters].set_ylim((-0.1,maxswe + maxswe*0.05))

    axs[0].legend(['measured','modeled'],loc='upper left')
    axs[0].set_ylabel('SWE [mm]')

    plt.suptitle('Validation at Measured Sites')
    plt.tight_layout()
    plt.subplots_adjust(top=0.92)

    snow.swe_meas = swe_meas
    snow.swe_mod = swe_mod

    snow._logger.info('saving figure to %svalidation%s.png'%(snow.figs_path,snow.name_append))
    plt.savefig('%svalidation%s.png'%(snow.figs_path,snow.name_append))
