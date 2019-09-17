
from snowav.utils.MidpointNormalize import MidpointNormalize
import numpy as np
from matplotlib import pyplot as plt
import matplotlib.colors as mcolors
from mpl_toolkits.axes_grid1 import make_axes_locatable
import seaborn as sns
import copy
import cmocean
import matplotlib.patches as mpatches
import pandas as pd
from datetime import timedelta
import netCDF4 as nc
from snowav.utils.wyhr import calculate_wyhr_from_date, calculate_date_from_wyhr
from snowav.database.database import collect, make_session
import snowav.framework.figures
from sys import exit

def flt_image_change(args, logger = None):
    '''
    Difference in SWE from one day prior to flight updates to the day of flight
    updates.

    The file specified in [plots] update_file is used to get dates of flight
    updates and the masks for which portions of the basin were updated. The
    actual SWE differences are between the snow.nc files on the dates of the
    flights and the day immediately prior. [plots] update_numbers is a 1-based
    list of subsetting flights that can be applied.

    Args
    ------
    args : dict
        dictionary with required inputs, see swi() figure for more information.

    logger : list
        snowav logger
    '''

    file = args['update_file']
    update_numbers = args['update_numbers']
    start_date = args['start_date']
    end_date = args['end_date']
    flight_outputs = args['flight_outputs']
    pre_flight_outputs = args['pre_flight_outputs']
    masks = args['masks']
    lims = args['lims']
    barcolors = args['barcolors']
    edges = args['edges']
    connector = args['connector']
    plotorder = args['plotorder']
    e = False
    flight_diff_fig_names = []
    flight_delta_vol_df = {}

    p = nc.Dataset(file, 'r')

    if update_numbers is None:
        times = p.variables['time'][:]
    else:
        if max(update_numbers) > len(p.variables['time']) - 1:
            raise Exception('Max flight update_numbers greater than available '
                            'flights')
        else:
            times = p.variables['time'][update_numbers]

    # remove flight indices that might be after the report period
    idx = []

    for i, n in enumerate(times):
        date = calculate_date_from_wyhr(int(n),args['wy'])
        if date > end_date:
            idx = np.append(idx,i)

    times = np.delete(times, idx)
    ix = np.argsort(times)
    times = times[ix]

    if update_numbers is not None:
        update_numbers = [int(x) for x in update_numbers]

    for i,time in enumerate(times):
        if update_numbers is not None:
            depth = p.variables['depth'][update_numbers[ix[i]],:,:]
        else:
            depth = p.variables['depth'][ix[i],:,:]

        delta_swe = flight_outputs['swe_z'][i][:] - pre_flight_outputs['swe_z'][i][:]
        delta_swe = delta_swe * args['depth_factor']

        # We will always make the difference between the previous day
        start_date = flight_outputs['dates'][i] - timedelta(hours = 24)
        end_date = flight_outputs['dates'][i]

        try:
            delta_swe_byelev = collect(connector, args['plotorder'], args['basins'],
                                       start_date, end_date,'swe_vol',
                                       args['run_name'], args['edges'],'difference')
        except:
            e = True

        try:
            end_swe_byelev = collect(connector, args['plotorder'], args['basins'],
                                     start_date, end_date, 'swe_vol',
                                     args['run_name'], args['edges'],'end')
        except:
            e = True

        if e:
            if logger is not None:
                logger.info(' Failed requesting database records ending on {} '
                            'for flight difference figure. This may mean that '
                            '[run] directory has not been processed with '
                            '[snowav] run_name: {} for the periods in '
                            '[plots] update_file. Try subsetting with [plots] '
                            'update_numbers or processing the full [run] '
                            'directory.'.format(end_date, args['run_name']))
                logger.info(' Flight figures being set to False...')

            # Returning empty lists will set flights to False
            return [], []

        # Make copy so that we can add nans for the plots
        delta_state = copy.deepcopy(delta_swe)
        vMin,vMax = np.nanpercentile(delta_state,[2,98])

        colorsbad = plt.cm.Set1_r(np.linspace(0., 1, 1))
        colors1 = cmocean.cm.matter_r(np.linspace(0., 1, 127))
        colors2 = plt.cm.Blues(np.linspace(0, 1, 128))
        colors = np.vstack((colorsbad,colors1, colors2))
        mymap = mcolors.LinearSegmentedColormap.from_list('my_colormap', colors)

        ixf = delta_state == 0
        delta_state[ixf] = -100000
        pmask = masks[plotorder[0]]['mask']
        ixo = pmask == 0
        delta_state[ixo] = np.nan
        cmap = copy.copy(mymap)
        cmap.set_bad('white',1.)

        sns.set_style('darkgrid')
        sns.set_context("notebook")

        plt.close(i)
        fig,(ax,ax1) = plt.subplots(num = i, figsize = args['figsize'],
                                    dpi = args['dpi'], nrows = 1, ncols = 2)

        norm = MidpointNormalize(midpoint=0, vmin=vMin-0.01, vmax=vMax+0.01)

        # this is primarily for awsm_test_cases, in which the flight update
        # nc file doesn't contain mask information. wy2019-forward flight nc
        # files *should* have a mask
        if hasattr(depth,'mask'):
            mask = np.ma.masked_array(depth.mask, ~depth.mask)
            if mask.shape != delta_state.shape:
                raise Exception('Dimensions {}: {} do not match snow.nc: {}'.format(
                                file, mask.shape, delta_state.shape))

            delta_state[mask] = np.nan
            h = ax.imshow(delta_state*mask, cmap = cmap, norm = norm)

        else:
            h = ax.imshow(delta_state, cmap = cmap, norm = norm)

        for name in masks:
            ax.contour(masks[name]['mask'],cmap = "Greys",linewidths = 1)

        h.axes.get_xaxis().set_ticks([])
        h.axes.get_yaxis().set_ticks([])
        divider = make_axes_locatable(ax)
        cax = divider.append_axes("right", size="5%", pad=0.2)
        cbar = plt.colorbar(h, cax = cax)
        cbar.set_label(r'$\Delta$ SWE [{}]'.format(args['depthlbl']))

        start_date = start_date + timedelta(hours=1)
        d = end_date + timedelta(hours=1)

        h.axes.set_title('Change in SWE from Snow Depth Update\n{}'
                         .format(start_date.date().strftime("%Y-%-m-%-d")))

        if args['dplcs'] == 0:
            tlbl = '{}: {} {}'.format(plotorder[0],
                                 str(int(delta_swe_byelev[plotorder[0]].sum())),
                                 args['vollbl'])
        else:
            tlbl = '{}: {} {}'.format(plotorder[0],
                                 str(np.round(delta_swe_byelev[plotorder[0]].sum(),
                                              args['dplcs'])),args['vollbl'])

        for iters,name in enumerate(lims.sumorder):
            if args['dplcs'] == 0:
                lbl = '{}: {} {}'.format(name,
                                    str(int(delta_swe_byelev[name].sum())),
                                    args['vollbl'])
            else:
                lbl = '{}: {} {}'.format(name,
                                    str(np.round(delta_swe_byelev[name].sum(),
                                    args['dplcs'])),args['vollbl'])

            if iters == 0:
                ax1.bar(range(0,len(args['edges'])),delta_swe_byelev[name],
                        color = barcolors[iters],
                        edgecolor = 'k',label = lbl)

            else:
                ax1.bar(range(0,len(args['edges'])),delta_swe_byelev[name],
                        color = barcolors[iters],
                        edgecolor = 'k',
                        label = lbl,
                        alpha = 0.5)

        end_swe_byelev = end_swe_byelev.fillna(0)
        datestr = flight_outputs['dates'][i].date().strftime("%Y%m%d")
        percent_delta_byelev = (delta_swe_byelev.sum(skipna=True)/end_swe_byelev.sum())*100
        delta_swe_byelev.to_csv('{}flight_{}_{}_taf.csv'.format(args['figs_path'],
                                                    args['directory'],datestr))
        percent_delta_byelev.to_csv('{}flight_{}_{}_percent.csv'.format(args['figs_path'],
                                                    args['directory'],datestr))

        flight_delta_vol_df[flight_outputs['dates'][i].date().strftime("%Y%m%d")] = delta_swe_byelev

        plt.tight_layout()
        xts = ax1.get_xticks()
        edges_lbl = []
        for x in xts[0:len(xts)-1]:
            edges_lbl.append(str(int(edges[int(x)])))

        ax1.set_xticklabels(str(x) for x in edges_lbl)
        for tick in ax1.get_xticklabels():
            tick.set_rotation(30)

        ylims = ax1.get_ylim()
        ymax = ylims[1] + abs(ylims[1] - ylims[0])
        ymin = ylims[0] - abs(ylims[1] - ylims[0])*0.1
        ax1.set_ylim((ymin, ymax))

        ax1.set_ylabel('{}'.format(args['vollbl']))
        ax1.set_xlabel('elevation [{}]'.format(args['elevlbl']))
        ax1.axes.set_title('Change in SWE')

        ax1.yaxis.set_label_position("right")
        ax1.tick_params(axis='x')
        ax1.tick_params(axis='y')
        ax1.yaxis.tick_right()

        if len(plotorder) > 1:
            ax1.legend(loc=(lims.legx,lims.legy))

        ax1.text(lims.btx,lims.bty,tlbl,horizontalalignment='center',
                 transform=ax1.transAxes,fontsize = 10)

        plt.tight_layout()
        fig.subplots_adjust(top=0.88)

        fig_name = '{}flight_{}_{}_diff.png'.format(args['figs_path'],args['directory'],datestr)
        flight_diff_fig_names += ['flight_{}_{}_diff.png'.format(args['directory'],datestr)]

        if logger is not None:
            logger.info(' saving {}'.format(fig_name))

        snowav.framework.figures.save_fig(fig, fig_name)

    p.close()

    return flight_diff_fig_names, flight_delta_vol_df
