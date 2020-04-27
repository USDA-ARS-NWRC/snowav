from datetime import datetime
from matplotlib import pyplot as plt
import mysql.connector
import netCDF4 as nc
import numpy as np
import os
import pandas as pd
import seaborn as sns
import utm

from snowav.utils.stats import nashsutcliffe
import snowav.framework.figures


def stn_validate(rundirs, lbls, client, end_date, wy, snow_x, snow_y, stns, px,
                 py, login, figs_path, fig_name, dem, logger=None, factor=25.4,
                 nash_sut_flag=False, ncfile='snow.nc', nc_var='specific_mass',
                 index_col='date_time', elevlbl='ft', tbl='tbl_level1',
                 var='snow_water_equiv', dpi=200):
    """ SWE validation at snow pillow sites.

    Args
    ----------
    rundirs {list}: list of rundirs
    lbls {list}: stn labels
    client {str}: weather db client
    end_date {datetime}: period end date
    wy {int}: water year
    snow_x {}: netcdf x vector
    snow_y {}: netcdf y vector
    stns {list}: list of stations
    px {arr}: +- 1 array
    py {arr}: +- 1 array
    login {dict}: database login info
    figs_path {str}: base path for fig saving
    fig_name {str}: base str for fig name
    dem {arr}: dem
    logger {class}: snowav logger
    nash_sut_flag {bool}: include nash sutcliffe
    ncfile {str}: name of snow.nc
    nc_var {str}: name of nc variable
    indec_col {str}: weather db index
    dpi {int}: fig dpi
    elevlbl {str}: elevation label
    tbl {str}: weather db table
    factor {float}: conversion from mm
    logger {class}: snowav logger

    Returns
    ------
    flag {bool}: if there are errors stn_validate() figure will be removed
        from the report
    """

    logger.debug(' Beginning stn_validate()...')

    st_time = '{}-10-1 00:00:00'.format(str(wy-1))
    end_time = end_date.date().strftime("%Y-%-m-%-d")
    flag = True
    wy_start = datetime(wy-1, 10, 1)
    pr = range(0, 9)

    # container for +- 1 pixel
    stns_pixel = []
    for sta in stns:
        stns_pixel = stns_pixel + ['{}_'.format(sta) + str(s) for s in pr]

    # containers for measure and model results
    measure = pd.DataFrame(index=pd.date_range(wy_start, end_date, freq='D'),
                           columns=stns)
    model = pd.DataFrame(index=pd.date_range(wy_start, end_date, freq='D'),
                         columns=stns_pixel)

    # weather db query
    qry = ('SELECT tbl_metadata.* FROM tbl_metadata INNER JOIN tbl_stations ON '
           'tbl_metadata.id = tbl_stations.metadata_id WHERE '
           'tbl_stations.client="{}";'.format(client))

    cnx = mysql.connector.connect(user=login['user'],
                                  password=login['password'],
                                  host=login['host'],
                                  port=login['port'],
                                  database='weather_db')

    metadata = pd.read_sql(qry, cnx)

    if metadata.empty:
        flag = False
        return

    metadata.index = metadata['primary_id']
    metadata = metadata[~metadata.index.duplicated(keep='first')]

    # station results
    for iters, stn in enumerate(stns):
        cnx = mysql.connector.connect(user=login['user'],
                                      password=login['password'],
                                      host=login['host'],
                                      port=login['port'],
                                      database='weather_db')

        qry = ('SELECT weather_db.{0}.date_time, weather_db.{0}.{1} FROM '
               'weather_db.{0} WHERE weather_db.{0}.date_time between "{2}" '
               'and "{3}" AND weather_db.{0}.station_id IN ("{4}")'
               ';'.format(tbl, var, st_time, end_time, stn))

        # data = pd.read_sql(qry, cnx, index_col=bytes(bytearray(b'date_time')))
        data = pd.read_sql(qry, cnx, index_col=index_col)
        data.index.names = [index_col]
        dind = pd.date_range(st_time, end_time, freq='D')
        measure[stn] = data.reindex(dind)

    sns.set_style('darkgrid')
    sns.set_context('notebook')

    plt.close(9)

    if len(stns) == 2:
        fig, axs = plt.subplots(num=9, figsize=(6, 4), dpi=dpi,
                                nrows=1, ncols=2)
    if len(stns) in [3, 4]:
        fig, axs = plt.subplots(num=9, figsize=(6, 6), dpi=dpi,
                                nrows=2, ncols=2)
    if (len(stns) > 4) and (len(stns) <= 6):
        fig, axs = plt.subplots(num=9, figsize=(10, 10), dpi=dpi,
                                nrows=3, ncols=2)
    if (len(stns) > 6) and (len(stns) <= 9):
        fig, axs = plt.subplots(num=9, figsize=(10, 10), dpi=dpi,
                                nrows=3, ncols=3)
    if (len(stns) > 9) and (len(stns) <= 12):
        fig, axs = plt.subplots(num=9, figsize=(10, 10), dpi=dpi,
                                nrows=4, ncols=3)

    axs = axs.flatten()

    set_x_on = 12
    if len(stns) == 3:
        fig.delaxes(axs[len(stns)])
    if len(stns) in (5, 8):
        fig.delaxes(axs[len(stns)])
        if len(stns) == 5:
            set_x_on = 3
        if len(stns) == 8:
            set_x_on = 5
    if len(stns) == 10:
        fig.delaxes(axs[10])
        fig.delaxes(axs[11])
    if len(stns) == 11:
        fig.delaxes(axs[11])

    for rname in rundirs:
        logger.debug(' Loading pixel values in '
                     '{}...'.format(rname.split('runs')[-1]))

        snowfile = os.path.join(rname, ncfile)

        if not os.path.isfile(snowfile):
            raise Exception('invalid file --> {}'.format(snowfile))

        ncf = nc.Dataset(snowfile, 'r')
        t = nc.num2date(ncf.variables['time'][0], ncf.variables['time'].units)

        # Only load the rundirs that we need
        if t.date() <= end_date.date():

            for iters, stn in enumerate(stns):
                ll = utm.from_latlon(metadata.ix[stn, 'latitude'],
                                     metadata.ix[stn, 'longitude'])
                xind = np.where(abs(snow_x - ll[0]) ==
                                min(abs(snow_x - ll[0])))[0]
                yind = np.where(abs(snow_y - ll[1]) ==
                                min(abs(snow_y - ll[1])))[0]
                elev = dem[yind, xind]

                if ((xind == 0) or (xind >= len(snow_x) - 1) or
                        (yind == 0) or (yind >= len(snow_y) - 1)):
                    logger.info('Station {} in [validate] stations is '
                                'outside of domain, exiting stn_validate() '
                                'and not generating figure'.format(stn))
                    flag = False
                    return flag

                for ix, (n, m) in enumerate(zip(px, py)):
                    land_stn = stn + '_{}'.format(str(ix))
                    swe = ncf.variables[nc_var][:, yind + m, xind + n][0][0][0]
                    model.loc[t.date(), land_stn] = swe

                    if (n == 0 and m == 0) and (rname == rundirs[-1]):
                        lbl = 'measured'
                        lblm = 'modeled'
                    else:
                        lbl = '__nolabel__'
                        lblm = 'modeled'

                    axs[iters].plot(measure[stn] / factor,
                                    'k',
                                    label=lbl)
                    axs[iters].plot(model[land_stn] / factor,
                                    'b',
                                    linewidth=0.75,
                                    label=lblm)
                    axs[iters].set_title('{} [{} {}]'.format(lbls[iters],
                                                             str(int(elev)),
                                                             elevlbl))
                    axs[iters].set_xlim((datetime(wy - 1, 10, 1), end_date))

            ncf.close()

        else:
            ncf.close()

    # Nash-Sutcliffe on pixel 0,0
    if nash_sut_flag:
        for iters, sta in enumerate(stns):
            pz_sta = sta + '_5'
            stime = datetime(wy - 1, 10, 15)
            nsv = nashsutcliffe(measure[measure.index > stime][stn].values.tolist(),
                                model[model.index > stime][pz_sta].values.tolist())
            nstr = 'Nash-Sutcliffe:\n{}'.format(str(np.round(nsv, 2)))

            nsx = 0.05
            nsy = 0.85

            if iters == 0:
                nsy = 0.675
            if nsv > 0:
                axs[iters].text(nsx, nsy, nstr, horizontalalignment='left',
                                transform=axs[iters].transAxes, fontsize=8)

    if len(stns) == 2:
        axs[1].yaxis.tick_right()
        for n in (0, 1):
            for tick in axs[n].get_xticklabels():
                tick.set_rotation(30)

    if len(stns) == 3:
        axs[1].yaxis.tick_right()
        axs[0].set_xticklabels('')
        for n in (0, 1, 2):
            for tick in axs[n].get_xticklabels():
                tick.set_rotation(30)

    if len(stns) == 4:
        axs[1].yaxis.tick_right()
        axs[3].yaxis.tick_right()
        axs[0].set_xticklabels('')
        axs[1].set_xticklabels('')
        for n in (0, 1, 2, 3):
            for tick in axs[n].get_xticklabels():
                tick.set_rotation(30)

    if 4 < len(stns) <= 6:
        for n in (1, 3, 5):
            axs[n].yaxis.tick_right()
        for n in (4, 5):
            for tick in axs[n].get_xticklabels():
                tick.set_rotation(30)
        for n in (0, 1, 2, 3):
            if n != set_x_on:
                axs[n].set_xticklabels('')
            else:
                for tick in axs[n].get_xticklabels():
                    tick.set_rotation(30)

    if 6 < len(stns) <= 9:
        for n in (2, 5, 8):
            axs[n].yaxis.tick_right()

        for n in (6, 7, 8):
            for tick in axs[n].get_xticklabels():
                tick.set_rotation(30)
        for n in (0, 1, 2, 3, 4, 5):
            axs[n].set_xticklabels('')
        for n in (1, 4, 7):
            axs[n].set_yticklabels('')

    if 9 < len(stns) <= 12:
        for n in (2, 5, 8, 11):
            axs[n].yaxis.tick_right()
        for n in (9, 10, 11):
            for tick in axs[n].get_xticklabels():
                tick.set_rotation(30)
        for n in (0, 1, 2, 3, 4, 5, 6, 7, 8):
            axs[n].set_xticklabels('')

    # Plot
    measure = measure.replace('[]', np.nan)
    maxm = np.nanmax(measure.max().values / factor)
    maxi = np.nanmax(model.max().values / factor)

    if maxm > maxi:
        maxswe = maxm
    else:
        maxswe = maxi

    for iters in range(0, len(stns)):
        axs[iters].set_ylim((-0.1, maxswe + maxswe * 0.05))

    axs[0].legend(['measured', 'modeled'], loc='upper left')
    axs[0].set_ylabel('SWE [in]')

    plt.suptitle('SWE Validation at Snow Pillow Sites')

    plt.tight_layout()

    if len(stns) <= 4:
        plt.subplots_adjust(top=0.85, bottom=0.15)

    if len(stns) > 4:
        plt.subplots_adjust(top=0.92)

    fig_name = os.path.join(os.path.abspath(figs_path), fig_name)

    snowav.framework.figures.save_fig(fig, fig_name)

    if logger is not None:
        logger.info(' Saved: {}'.format(fig_name))

    return flag
