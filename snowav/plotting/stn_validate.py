
import numpy as np
from matplotlib import pyplot as plt
import seaborn as sns
import utm
import netCDF4 as nc
import mysql.connector
import pandas as pd
from datetime import datetime
import os
from snowav.utils.stats import nashsutcliffe
import snowav.framework.figures

def stn_validate(args, logger):
    '''
    SWE validation at snow pillow sites.

    Requires proper config options for weather database connection and figure
    creation, see CoreConfig.ini for more information.

    Args
    ----------
    args : dict
        dictionary with required inputs, see swi() figure for more information.
    logger : list

    Returns
    ------
    logger : list
    flag : bool
        False if there are errors, will be used to remove stn_validate() figure
        from the report
    '''

    logger.debug(' Beginning stn_validate()...')

    rundirs = args['dirs']
    stns = args['stns']
    lbls = args['lbls']
    client = args['client']
    factor = args['factor']
    end_date = args['end_date']
    tbl = args['tbl']
    var = args['var']
    st_time = '{}-10-1 00:00:00'.format(str(args['wy'] - 1))
    end_time = end_date.date().strftime("%Y-%-m-%-d")
    ncxvec = args['snow_x']
    ncyvec = args['snow_y']
    stns = args['stns']
    fig_name_short = 'validation_'
    flag = True
    px = args['px']
    py = args['py']

    stns_pixel = []
    for sta in stns:
        stns_pixel = stns_pixel + ['{}_'.format(sta) + str(s) for s in range(0,9)]

    measure = pd.DataFrame(index = pd.date_range(datetime(args['wy']-1,10,1),
                            end_date,freq='D'), columns = stns)

    model = pd.DataFrame(index = pd.date_range(datetime(args['wy']-1,10,1),
                         end_date,freq='D'), columns = stns_pixel)

    qry = ('SELECT tbl_metadata.* FROM tbl_metadata INNER JOIN tbl_stations ON '
          'tbl_metadata.id = tbl_stations.metadata_id WHERE tbl_stations.client'
          '="{}";'.format(client))

    cnx = mysql.connector.connect(user=args['user'], password=args['password'],
                                  host=args['host'], port=args['port'],
                                  database='weather_db')

    metadata = pd.read_sql(qry, cnx)

    if metadata.empty:
        flag = False
        return

    metadata.index = metadata['primary_id']
    metadata = metadata[~metadata.index.duplicated(keep='first')]

    # station results
    for iters,stn in enumerate(stns):
        cnx = mysql.connector.connect(user=args['user'],password=args['password'],
                                      host=args['host'], port=args['port'],
                                      database='weather_db')

        qry = ('SELECT weather_db.{0}.date_time, weather_db.{0}.{1} FROM '
               'weather_db.{0} WHERE weather_db.{0}.date_time between "{2}" '
               'and "{3}" AND weather_db.{0}.station_id IN '
               ' ("{4}") ;'.format(args['tbl'],args['var'],st_time,end_time,stn))

        # data = pd.read_sql(qry, cnx, index_col=bytes(bytearray(b'date_time')))
        data = pd.read_sql(qry, cnx, index_col='date_time')
        data.index.names=['date_time']
        dind = pd.date_range(st_time,end_time,freq='D')
        measure[stn] = data.reindex(dind)

    sns.set_style('darkgrid')
    sns.set_context('notebook')

    plt.close(9)

    if len(stns) <= 3:
        fig, axs = plt.subplots(num=9,figsize=(6,6),dpi=args['dpi'],nrows=3,ncols=1)

    if (len(stns) > 3) and (len(stns) <= 6):
        fig, axs = plt.subplots(num=9,figsize=(10,10),dpi=args['dpi'],nrows=3,ncols=2)

    if (len(stns) > 6) and (len(stns) <= 9):
        fig, axs = plt.subplots(num=9,figsize=(10,10),dpi=args['dpi'],nrows=3,ncols=3)

    if (len(stns) > 9) and (len(stns) <= 12):
        fig, axs = plt.subplots(num=9,figsize=(10,10),dpi=args['dpi'],nrows=4,ncols=3)

    axs = axs.flatten()

    set_x_on = 12
    if len(stns) in (2,5,8):
        fig.delaxes(axs[len(stns)])
        if len(stns) == 5:
            set_x_on = 3
        if len(stns) == 8:
            set_x_on = 5

    for rname in rundirs:
        logger.debug(' Loading pixel values in {}...'.format(rname.split('runs')[-1]))

        # could add if args['ncfile'] == 'em.nc' change to data/smrf/precip
        snowfile = os.path.join(rname, args['ncfile'])

        if not os.path.isfile(snowfile):
            raise Exception('In stn_validate(), {} not a valid file'.format(snowfile))

        ncf = nc.Dataset(snowfile, 'r')
        t = nc.num2date(ncf.variables['time'][0],ncf.variables['time'].units)

        # Only load the rundirs that we need
        if (t.date() <= end_date.date()):

            for iters,stn in enumerate(stns):
                ll = utm.from_latlon(metadata.ix[stn,'latitude'],metadata.ix[stn,'longitude'])
                xind = np.where(abs(ncxvec-ll[0]) == min(abs(ncxvec-ll[0])))[0]
                yind = np.where(abs(ncyvec-ll[1]) == min(abs(ncyvec-ll[1])))[0]

                if ((xind == 0) or (xind >= len(ncxvec)-1) or
                   (yind == 0) or (yind >= len(ncyvec)-1)):
                    logger.info('Station {} in [validate] stations is '
                                'outside of domain, exiting stn_validate() '
                                'and not generating figure'.format(stn))
                    flag = False
                    return '', flag

                for ix, (n,m) in enumerate(zip(px,py)):
                    land_stn = stn + '_{}'.format(str(ix))
                    swe = ncf.variables['specific_mass'][:,yind+m,xind+n][0][0][0]
                    model.loc[t.date(),land_stn] = swe

                    if (n == 0 and m == 0) and (rname == rundirs[-1]):
                        lbl = 'measured'
                        lblm = 'modeled'
                    else:
                        lbl = '__nolabel__'
                        lblm = 'modeled'

                    axs[iters].plot(measure[stn]/factor,'k',label=lbl)
                    axs[iters].plot(model[land_stn]/factor,'b',linewidth = 0.75,label=lblm)
                    axs[iters].set_title(lbls[iters])
                    axs[iters].set_xlim((datetime(args['wy']-1, 10, 1),end_date))

            ncf.close()

        else:
            ncf.close()

    # Nash-Sutcliffe on pixel 0,0
    if args['nash_sut_flag']:
        for iters,sta in enumerate(stns):
            pz_sta = sta + '_5'
            stime = datetime(args['wy']-1,10,15)
            nsv = nashsutcliffe(measure[measure.index>stime][stn].values.tolist(),
                                model[model.index>stime][pz_sta].values.tolist())
            # nstr = '{} Nash-Sutcliffe:\n{}'.format(sta,str(np.round(nsv,2)))
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
            if n != set_x_on:
                axs[n].set_xticklabels('')
            else:
                for tick in axs[n].get_xticklabels():
                    tick.set_rotation(30)

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
    measure = measure.replace('[]',np.nan)
    maxm = np.nanmax(measure.max().values/factor)
    maxi = np.nanmax(model.max().values/factor)

    if maxm > maxi:
        maxswe = maxm
    else:
        maxswe = maxi

    for iters in range(0,len(stns)):
        axs[iters].set_ylim((-0.1,maxswe + maxswe*0.05))

    axs[0].legend(['measured','modeled'],loc='upper left')
    axs[0].set_ylabel('SWE [in]')

    plt.suptitle('SWE Validation at Snow Pillow Sites')
    plt.subplots_adjust(top=0.92)

    fig_name = '{}{}{}.png'.format(args['figs_path'],fig_name_short,args['directory'])
    if logger is not None:
        logger.info(' saving {}'.format(fig_name))

    snowav.framework.figures.save_fig(fig, fig_name)

    return fig_name_short, flag
