import copy
import cmocean
from datetime import timedelta
from matplotlib import pyplot as plt
import matplotlib.colors as mcolors
from mpl_toolkits.axes_grid1 import make_axes_locatable
import netCDF4 as nc
import numpy as np
import os
import seaborn as sns

from snowav.utils.wyhr import calculate_date_from_wyhr
from snowav.utils.MidpointNormalize import MidpointNormalize
from snowav.database.database import collect
import snowav.framework.figures


def flt_image_change(file, update_numbers, end_date, flight_outputs,
                     pre_flight_outputs, masks, lims, barcolors, edges,
                     connector, plotorder, wy, depth_factor, basins,
                     run_name, figsize, depthlbl, elevlbl, vollbl, dplcs,
                     figspath, dpi=200, logger=None):
    """ Difference in SWE from one day prior to flight updates to the day of
    flight updates.

    The file specified in [plots] update_file is used to get dates of flight
    updates and the masks for which portions of the basin were updated. The
    actual SWE differences are between the snow.nc files on the dates of the
    flights and the day immediately prior. [plots] update_numbers is a 1-based
    list of subsetting flights that can be applied.

    Args
    ------
    file {str}: path to lidar_depths.nc
    update_numbers {list}: flight numbers to generate
    end_date {datetime}: end date
    flight_outputs {dict}: from config.py
    pre_flight_outputs {dict}: from config.py
    masks {dict}: from config.py
    lims {list}: percent of min/max for color limits
    barcolors {list}: colors
    edges {list}: elevation bands
    connector {str}: database connector
    plotorder {list}: basins
    wy {int}: water year
    depth_factor {float}: depth factor
    basins {dict}: from config.py
    run_name {str}: snowav run_name
    figsize {list}: figure dimensions
    depthlbl {str}: depth label
    elevlbl {str}: elevation label
    vollbl {str}: volume label
    dplcs {int}: decimal places
    figspath {str}: base path for figures
    dpi {int}: figure dpi
    """

    e = False
    basin_name = plotorder[0].lower().split(" ")[0]
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
        date = calculate_date_from_wyhr(int(n), wy)
        if date > end_date:
            idx = np.append(idx, i)

    times = np.delete(times, idx)
    ix = np.argsort(times)
    times = times[ix]

    if update_numbers is not None:
        update_numbers = [int(x) for x in update_numbers]

    for i, time in enumerate(times):
        if update_numbers is not None:
            depth = p.variables['depth'][update_numbers[ix[i]], :, :]
        else:
            depth = p.variables['depth'][ix[i], :, :]

        delta_swe = flight_outputs['swe_z'][i][:] - pre_flight_outputs['swe_z'][i][:]
        delta_swe = delta_swe * depth_factor

        # We will always make the difference between the previous day
        start_date = flight_outputs['dates'][i] - timedelta(hours=24)
        end_date = flight_outputs['dates'][i]

        try:
            delta_swe_byelev = collect(connector, plotorder, basins, start_date,
                                       end_date, 'swe_vol', run_name, edges,
                                       'difference')
        except:
            e = True

        try:
            end_swe_byelev = collect(connector, plotorder, basins, start_date,
                                     end_date, 'swe_vol', run_name, edges,
                                     'end')
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
                            'directory.'.format(end_date, run_name))
                logger.info(' Flight figures being set to False...')

            return [], []

        # Make copy so that we can add nans for the plots
        delta_state = copy.deepcopy(delta_swe)
        v_min, v_max = np.nanpercentile(delta_state, [lims[0], lims[1]])

        colorsbad = plt.cm.Set1_r(np.linspace(0., 1, 1))
        colors1 = cmocean.cm.matter_r(np.linspace(0., 1, 127))
        colors2 = plt.cm.Blues(np.linspace(0, 1, 128))
        colors = np.vstack((colorsbad, colors1, colors2))
        mymap = mcolors.LinearSegmentedColormap.from_list('my_colormap', colors)

        ixf = delta_state == 0
        delta_state[ixf] = -100000
        pmask = masks[plotorder[0]]['mask']
        ixo = pmask == 0
        delta_state[ixo] = np.nan
        cmap = copy.copy(mymap)
        cmap.set_bad('white', 1.)

        sns.set_style('darkgrid')
        sns.set_context("notebook")

        plt.close(i)
        fig, (ax, ax1) = plt.subplots(num=i, figsize=figsize, dpi=dpi, nrows=1,
                                      ncols=2)

        norm = MidpointNormalize(midpoint=0, vmin=v_min, vmax=v_max)

        # this is primarily for awsm_test_cases, in which the flight update
        # nc file doesn't contain mask information. wy2019-forward flight nc
        # files *should* have a mask
        if hasattr(depth, 'mask'):
            mask = np.ma.masked_array(depth.mask, ~depth.mask)
            if mask.shape != delta_state.shape:
                raise Exception('Dimensions {}: {} do not match snow.nc: '
                                '{}'.format(file, mask.shape, delta_state.shape))

            delta_state[mask] = np.nan
            h = ax.imshow(delta_state * mask, cmap=cmap, clim=(v_min, v_max),
                          norm=norm)

        else:
            h = ax.imshow(delta_state, cmap=cmap, clim=(v_min, v_max),
                          norm=norm)

        for name in masks:
            ax.contour(masks[name]['mask'], cmap="Greys", linewidths=1)

        h.axes.get_xaxis().set_ticks([])
        h.axes.get_yaxis().set_ticks([])
        divider = make_axes_locatable(ax)
        cax = divider.append_axes("right", size="5%", pad=0.2)
        cbar = plt.colorbar(h, cax=cax)
        cbar.set_label(r'$\Delta$ SWE [{}]'.format(depthlbl))


        end_date = start_date + timedelta(hours=1)
        d = end_date + timedelta(hours=1)

        h.axes.set_title('Change in SWE Depth\n{} to {}'
                         .format(start_date.date().strftime("%Y-%-m-%-d"),
                                 end_date.date().strftime("%Y-%-m-%-d")))

        if len(plotorder) == 1:
            porder = plotorder
        else:
            porder = plotorder[1::]

        for iters, name in enumerate(porder):
            if dplcs == 0:
                lbl = '{}: {} {}'.format(name,
                                         str(int(delta_swe_byelev[name].sum())),
                                         vollbl)
            else:
                lbl = '{}: {} {}'.format(name,
                                         str(np.round(delta_swe_byelev[name].sum(),
                                                      dplcs)), vollbl)

            if iters == 0:
                ax1.bar(range(0, len(edges)), delta_swe_byelev[name],
                        color=barcolors[iters],
                        edgecolor='k', label=lbl)

            else:
                ax1.bar(range(0, len(edges)), delta_swe_byelev[name],
                        color=barcolors[iters],
                        edgecolor='k',
                        label=lbl,
                        alpha=0.5)

        end_swe_byelev = end_swe_byelev.fillna(0)
        datestr = flight_outputs['dates'][i].date().strftime("%Y%m%d")
        percent_delta_byelev = (delta_swe_byelev.sum(skipna=True) / end_swe_byelev.sum()) * 100
        flight_delta_vol_df[flight_outputs['dates'][i].date().strftime("%Y%m%d")] = delta_swe_byelev

        fp = os.path.join(figspath,
                          '{}_flight_{}_delta_taf.csv'.format(basin_name, datestr))
        delta_swe_byelev.to_csv(fp)
        fp = os.path.join(figspath,
                          '{}_flight_{}_percent_change.csv'.format(basin_name, datestr))
        percent_delta_byelev.to_csv(fp)

        plt.tight_layout()
        xts = ax1.get_xticks()
        edges_lbl = []
        for x in xts[0:len(xts) - 1]:
            edges_lbl.append(str(int(edges[int(x)])))

        ax1.set_xticklabels(str(x) for x in edges_lbl)
        for tick in ax1.get_xticklabels():
            tick.set_rotation(30)

        ylims = ax1.get_ylim()
        ymax = ylims[1] + abs(ylims[1] - ylims[0])
        ymin = ylims[0] - abs(ylims[1] - ylims[0]) * 0.1
        ax1.set_ylim((ymin, ymax))

        ax1.set_ylabel(r'$\Delta$ SWE [{}]'.format(vollbl))
        ax1.set_xlabel('elevation [{}]'.format(elevlbl))

        ax1.yaxis.set_label_position("right")
        ax1.tick_params(axis='x')
        ax1.tick_params(axis='y')
        ax1.yaxis.tick_right()

        if len(plotorder) > 1:
            ax1.legend(loc=2, fontsize=10)

        plt.tight_layout()
        fig.subplots_adjust(top=0.88)

        fstr = '{}_flight_difference_{}.png'.format(basin_name, datestr)
        flight_diff_fig_names += [fstr]
        fig_path = os.path.join(figspath, fstr)

        snowav.framework.figures.save_fig(fig, fig_path)

        if logger is not None:
            logger.info(' Saved: {}'.format(fig_path))

    p.close()

    return flight_diff_fig_names, flight_delta_vol_df
