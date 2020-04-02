import cmocean
import copy
from datetime import datetime
from matplotlib import pyplot as plt
import matplotlib.colors
from mpl_toolkits.axes_grid1 import make_axes_locatable
from mpl_toolkits.axes_grid1.inset_locator import inset_axes
import numpy as np
import pandas as pd
import os
import seaborn as sns
import utm

import snowav.framework.figures


def point_values_csv(output, value, df, imgxy, filename, model_date, basin,
                     heading, logger):
    """Read in specified csv file of snow course validation locations, write
    model values at those points.

    Args
    ------
    output {arr}: model output array from date of interest
    value {str}: snowav database value
    df {DataFrame}: DataFrame from point_values_csv file
    imgxy {arr}: x, y coordinates for the image, from get_topo_stats.py
    filename {str}: file name
    model_date {str}: date that model values were pulled
    logger {class}: logger
    """

    if value == 'swe_z':
        header = 'SWE [in]'
        intro_col = 'Snow Water Equivalent (SWE) depth in inches'
    if value == 'depth':
        header = 'snow depth [in]'
        intro_col = 'Snow depth in inches'
    if value == 'density':
        header = 'density [kg/m^3]'
        intro_col = 'Snow density in kg/m^3'

    date_str = datetime.now().date().strftime("%Y-%m-%d")
    date_col = 'Date generated: {}'.format(date_str)
    basin_col = 'Model basin: {}'.format(basin)
    model_col = 'Model results from: {}'.format(model_date)
    med_col = '9 pixel median {}'.format(header)
    lat_col = 'latitude'
    lon_col = 'longitude'
    x_col = 'model x'
    y_col = 'model y'
    stn_col = 'name'
    out_cols = [stn_col, lat_col, lon_col, x_col, y_col, header, med_col]
    all_cols = [stn_col, x_col, y_col, lat_col, lon_col, header, med_col,
                heading]
    sic = 'The "{}" column refers the model value in the pixel that ' \
          'contains the "{}" "{}" point'.format(header, lat_col, lon_col)
    npmc = 'The "{}" column means the median value ' \
           'of the 9 closest pixels to the "{}" "{}" ' \
           'point'.format(med_col, lat_col, lon_col)

    file_header = ['USDA Agicultural Research Service Snowpack Summary Data',
                   intro_col,
                   'Data provided are model results from the iSnobal model',
                   basin_col,
                   model_col,
                   date_col,
                   'The "model x" column refers to the model column index',
                   'The "model y" column refers to the model row index',
                   sic,
                   npmc,
                   'Contact: Scott Havens <scott.havens@usda.gov>\n']

    pixel_swe = pd.DataFrame(columns=all_cols)

    # get closest model pixel
    for i in df.index:
        sta = df.loc[i, stn_col]
        lat = df.loc[i, lat_col]
        lon = df.loc[i, lon_col]

        ll = utm.from_latlon(lat, lon)
        xind = np.where(abs(imgxy[0] - ll[0]) == min(abs(imgxy[0] - ll[0])))[0]
        yind = np.where(abs(imgxy[1] - ll[1]) == min(abs(imgxy[1] - ll[1])))[0]

        row_lims = [0, len(output[:, 0]) - 1]
        col_lims = [0, len(output[0, :]) - 1]

        if xind in col_lims or yind in row_lims:
            logger.warn(' Skipping {} in point_values, latitude and/or '
                        'longitude may be outside of the model '
                        'domain'.format(sta))
        else:
            pixel_swe.loc[i, stn_col] = sta
            pixel_swe.loc[i, lat_col] = lat
            pixel_swe.loc[i, lon_col] = lon
            pixel_swe.loc[i, x_col] = int(xind)
            pixel_swe.loc[i, y_col] = int(yind)
            pixel_swe.loc[i, heading] = df.loc[i, heading]

            # Get the median of the nine closest pixels
            val = np.array([])
            for n in range(-1, 1):
                for m in range(-1, 1):
                    val = np.append(val, output[yind + n, xind + m])
                    if n == 0 and m == 0:
                        pixel_swe.loc[i, header] = \
                            np.round(float(output[yind + n, xind + m]),
                                     decimals=2)

            pixel_swe.loc[i, med_col] = \
                np.round(np.nanmedian(val), decimals=2)

    pixel_swe = pixel_swe.round(4)

    if os.path.isfile(filename):
        os.remove(filename)

    with open(filename, mode='w', encoding='utf-8') as f:
        f.write('\n'.join(file_header))

    pixel_swe.to_csv(filename, encoding='utf-8', mode='a', columns=out_cols,
                     index=False)
    logger.info(" Saved: {}".format(filename))

    return pixel_swe


def point_values_figures(output, value, df, dem, figs_path, veg_type,
                         model_date, course_date, settings, pixel, heading,
                         logger):
    """ Make point values figures.

    Args
    ------
    output {arr}: model output array from date of interest
    value {str}: snowav database value
    df {DataFrame}: DataFrame from point_values_csv function
    dem {arr}: dem
    figs_path {str}: base path for saving
    veg_type {arr}: veg type
    model_date {str}: date that model values were pulled
    course_date {str}: snow course data
    settings {list}: figure settings
    pixel {int}: pixel sixe
    heading {string}: df heading
    logger {class}: logger
    """

    if value == 'swe_z':
        header = 'SWE [in]'
    if value == 'depth':
        header = 'snow depth [in]'
    if value == 'density':
        header = 'density [kg/m^3]'

    x_col = 'model x'
    y_col = 'model y'
    stn_col = 'name'
    width = settings[0]
    height = settings[1]
    dpi = settings[2]
    nrows = settings[3]
    ncols = settings[4]
    font_small = settings[5]
    font_medium = settings[6]
    npix = settings[7]
    ss = settings[8]
    levels = settings[9]
    annot_x1 = settings[10]
    annot_y1 = settings[11]
    annot_x2 = settings[12]
    annot_y2 = settings[13]

    sns.set_style('white')
    plt.close(0)
    f, a = plt.subplots(num=0, figsize=(int(width), int(height)),
                        dpi=dpi, nrows=nrows, ncols=ncols)
    a = a.flatten()

    plt.close(1)
    f1, a1 = plt.subplots(num=1, figsize=(int(width), int(height)),
                          dpi=dpi, nrows=nrows, ncols=ncols)
    a1 = a1.flatten()

    plt.close(2)
    f2 = plt.figure(num=2, figsize=(width, height - 2), dpi=dpi)
    a2 = plt.gca()

    plt.close(3)
    f3 = plt.figure(num=3, figsize=(width, height - 2), dpi=400)
    a3 = plt.gca()

    lvls = list(np.arange(np.min(dem.flatten()), np.max(dem.flatten()), 500))

    snow_map = cmocean.cm.haline_r

    veg = []
    for n, idx in enumerate(df.index):

        ix = int(df.loc[idx][x_col])
        iy = int(df.loc[idx][y_col])

        veg_sub = veg_type[(iy - npix):(iy + npix + 1),
                  (ix - npix):(ix + npix + 1)].flatten()
        veg = np.append(veg_sub, veg)

        if n == 0:
            g_min = np.nanmin(output[(iy - npix):(iy + npix),
                              (ix - npix):(ix + npix)])
            g_max = np.nanmax(output[(iy - npix):(iy + npix),
                              (ix - npix):(ix + npix)])

        else:
            vmin = np.nanmin(output[(iy - npix):(iy + npix),
                             (ix - npix):(ix + npix)])
            vmax = np.nanmax(output[(iy - npix):(iy + npix),
                             (ix - npix):(ix + npix)])

            if vmin < g_min:
                g_min = vmin
            if vmax > g_max:
                g_max = vmax

    ymin = copy.deepcopy(g_min)
    veg = np.unique(veg)
    vl = list(veg)
    norm = matplotlib.colors.BoundaryNorm(vl, len(vl))

    h3 = a3.imshow(veg_type, cmap='Vega20', norm=norm, alpha=0.75)
    a3.contour(dem, colors='k', levels=lvls, linewidths=0.15)

    place = np.linspace(g_min, g_max, 256)

    for idx, v in enumerate(df.name.unique()):
        mflag = False

        if idx < len(df[stn_col].unique()) + 1:

            # this plots repeating values
            if len(df[df[stn_col] == v][stn_col].values) > 1:
                mflag = True
                vstr = mstr = ''

                for idx3 in range(0, len(df[df[stn_col] == v][stn_col].values)):
                    course = df.loc[idx + idx3][heading]
                    ix = int(df[df[stn_col] == v][x_col].values[idx3])
                    iy = int(df[df[stn_col] == v][y_col].values[idx3])
                    ixs = int(np.where(abs(place - course) ==
                                       min(abs(place - course)))[0])
                    color = snow_map(ixs)

                    if idx3 == 0:
                        vstr = '{}'.format(round(course, 1))
                        mstr = '{}'.format(output[iy, ix].round(1))
                    else:
                        vstr = vstr + ', {}'.format(round(course, 1))
                        mstr = mstr + ', {}'.format(output[iy, ix].round(1))

                    a[idx].scatter(x=ix, y=iy, marker='^', s=ss, c=color,
                                   edgecolors='k', linewidths=0.5)

                cstr = 'validation: {}\nmodel: {}'.format(vstr, mstr)
                a[idx].text(annot_x1,
                            annot_y1,
                            cstr,
                            horizontalalignment='left',
                            transform=a[idx].transAxes,
                            fontsize=font_small - 1,
                            color='w')

            ix = int(df[df[stn_col] == v][x_col].values[0])
            iy = int(df[df[stn_col] == v][y_col].values[0])
            swe = df[df[stn_col] == v][header].values[0]
            course = df[df[stn_col] == v][heading].values[0]

            s = dem[(iy - npix):(iy + npix + 1),
                    (ix - npix):(ix + npix + 1)]
            snow_sec = output[(iy - npix):(iy + npix + 1),
                              (ix - npix):(ix + npix + 1)]

            l1 = [int(x) * snow_sec.shape[0] +
                  npix - 1 for x in range(npix - 1, npix + 2)] + \
                 [int(x) * snow_sec.shape[0] +
                  npix for x in range(npix - 1, npix + 2)] + \
                 [int(x) * snow_sec.shape[0] +
                  npix + 1 for x in range(npix - 1, npix + 2)]

            if idx == 0:
                sl = r'$\pm${} pixels'.format(npix)
                bl1 = r'$\pm$1 pixel'.format(npix)
            else:
                sl = '__nolabel__'
                bl1 = '__nolabel__'

            b1 = a2.boxplot(snow_sec.flatten(), positions=[idx - 0.22],
                            widths=0.18, showfliers=False)
            b2 = a2.boxplot(snow_sec.flatten()[l1], positions=[idx + 0.22],
                            widths=0.18, showfliers=False)
            a2.plot(0, ymin, 'dimgray', label=sl)
            a2.plot(0, ymin, 'forestgreen', label=bl1)
            if idx == 0:
                astr = 'val-mod [in], %:\n{}\n' \
                       '{}%'.format(round((course - swe), 1),
                                    int((swe / course) * 100))
            else:
                astr = '{}\n{}%'.format(round((course - swe), 1),
                                        int((swe / course) * 100))

            a2.annotate(astr, (idx - 0.25, ymin - 2), fontsize=font_small,
                        color='k')

            for prop in ['boxes', 'whiskers', 'fliers', 'caps']:
                plt.setp(b1[prop], color='dimgray')
                plt.setp(b2[prop], color='forestgreen')

            lvls = list(np.arange(np.min(s.flatten()),
                                  np.max(s.flatten()), levels))
            ixs = int(np.where(abs(place - course) ==
                               min(abs(place - course)))[0])

            color = snow_map(ixs)
            a[idx].contour(dem, colors='k', levels=lvls, linewidths=0.25)
            h = a[idx].imshow(output, cmap=snow_map, clim=(g_min, g_max))

            if not mflag:
                a[idx].scatter(x=ix, y=iy, marker='^', s=ss, c=color,
                               edgecolors='k', linewidths=0.5)
                cstr = 'validation: {}\nmodel: ' \
                       '{}'.format(round(course, 1),
                                   output[iy, ix].round(1))
                a[idx].text(annot_x1, annot_y1, cstr,
                            horizontalalignment='left',
                            transform=a[idx].transAxes,
                            fontsize=font_small,
                            color='w')

            a[idx].get_xaxis().set_ticks([])
            a[idx].get_yaxis().set_ticks([])
            a[idx].set_ylim((iy + npix + 0.5), (iy - (npix + 0.5)))
            a[idx].set_xlim((ix - (npix + 0.5), ix + npix + 0.5))
            a[idx].set_title(v, fontsize=font_small)

            # legend of sorts
            veg_subp = veg_type[(iy - npix):(iy + npix + 1),
                       (ix - npix):(ix + npix + 1)]
            vegmap = plt.cm.get_cmap('Vega20')

            a1[idx].contour(dem, colors='k', levels=lvls, linewidths=0.25)
            a1[idx].imshow(veg_type, cmap='Vega20', norm=norm)
            a1[idx].scatter(x=ix, y=iy, marker='^', s=ss, c='k',
                            edgecolors='k', linewidths=0.5)
            a1[idx].get_xaxis().set_ticks([])
            a1[idx].get_yaxis().set_ticks([])

            a1[idx].set_ylim((iy + npix + 0.5), (iy - (npix + 0.5)))
            a1[idx].set_xlim((ix - (npix + 0.5), ix + npix + 0.5))
            a1[idx].set_title(v, fontsize=font_small)

            a3.scatter(x=ix, y=iy, marker='s', s=10, color="None",
                       edgecolors='k', linewidths=0.5)
            a3.annotate(v,
                        (ix + 5, iy - 5),
                        fontsize=font_small - 1,
                        color='k')

            for vs in np.unique(veg_subp):
                ixf = np.where(veg_subp == vs)
                fx = ixf[0][0]
                fy = ixf[1][0]
                nx = ix + fy - npix
                ny = iy + fx - npix
                a1[idx].scatter(x=nx, y=ny, marker='X', s=3, c='k',
                                edgecolors='k', linewidths=0.5)
                a1[idx].annotate(str(vs), (nx, ny), fontsize=font_small - 3,
                                 color='w')

    # tack on another subplot for legend and colorbar
    a[idx + 1].imshow(output, cmap=snow_map, alpha=0)
    a[idx + 1].set_title('', fontsize=font_small)
    a[idx + 1].get_xaxis().set_ticks([])
    a[idx + 1].get_yaxis().set_ticks([])
    a[idx + 1].axis('off')
    cbaxes = inset_axes(a[idx + 1], width="90%", height="10%", loc=8)
    cbar = f.colorbar(h,
                      cax=cbaxes,
                      ticks=[int(g_min + 1),
                             int(((g_max - g_min) / 2) + g_min),
                             int(g_max)],
                      orientation='horizontal')
    cbar.ax.tick_params(labelsize=font_medium)
    cbar.set_label(header, fontsize=font_medium)
    cbar.ax.xaxis.set_label_position('top')
    cstr = (r'Validation pixel $\pm${} pixels'
            '\nmodel: {}\nvalidation: {}\ncontour interval: {} '
            'ft\npixel width: {} m, {} '
            'ft'.format(npix, model_date, course_date,
                        levels, str(pixel), str(int(pixel * 3.28))))
    a[idx + 1].text(annot_x2, annot_y2, cstr, horizontalalignment='left',
                    transform=a[idx + 1].transAxes, fontsize=font_medium)

    a1[idx + 1].imshow(veg_type,
                       cmap=vegmap,
                       clim=(np.min(veg), np.max(veg)),
                       alpha=0)
    a1[idx + 1].set_title('', fontsize=font_small)
    a1[idx + 1].get_xaxis().set_ticks([])
    a1[idx + 1].get_yaxis().set_ticks([])
    a1[idx + 1].axis('off')
    cstr = ('model: {}\nvalidation: {}\ncontour interval: {} '
            'ft\npixel width: {} m, {} '
            'ft'.format(model_date, course_date,
                        levels, str(pixel), str(int(pixel * 3.28))))
    a1[idx + 1].text(annot_x2, annot_y2, cstr, horizontalalignment='left',
                     transform=a1[idx + 1].transAxes, fontsize=font_medium)

    for n in range(idx + 2, nrows * ncols):
        f.delaxes(a[n])
    for n in range(idx + 2, nrows * ncols):
        f1.delaxes(a1[n])

    a2.set_ylabel(header)
    a2.set_xticks(list(range(0, len(df[stn_col].unique()))))
    a2.set_xticklabels(df[stn_col].unique(), rotation=90,
                       fontsize=font_small - 1)
    a2.set_xlim((-1, len(df[stn_col].unique()) + 1))
    a2.legend(loc='upper left')
    a2.grid(linewidth=0.25)
    a2.set_title('Validation and Model Pixel Values, {}'.format(model_date))

    a3.get_xaxis().set_ticks([])
    a3.get_yaxis().set_ticks([])
    divider3 = make_axes_locatable(a3)
    cax3 = divider3.append_axes("right", size="4%", pad=0.2)
    cbar3 = plt.colorbar(h3, cax=cax3, ticks=list(veg))
    cbar3.ax.tick_params()
    cbar3.set_label('vegetation type classification')

    f.tight_layout()
    f1.tight_layout()
    f2.tight_layout()
    f3.tight_layout()

    fig_name_short = 'validation_list_{}'.format(value)
    fig_name = '{}{}.png'.format(figs_path, fig_name_short)
    logger.info(' Saved: {}'.format(fig_name))
    snowav.framework.figures.save_fig(f2, fig_name)

    # currently only passing the swe_z image
    if value == 'swe_z':
        fig_name_short = 'validation_map_{}'.format(value)
        fig_name = '{}{}.png'.format(figs_path, fig_name_short)
        logger.info(' Saved: {}'.format(fig_name))
        snowav.framework.figures.save_fig(f, fig_name)

    fig_name_short = 'validation_veg_map'
    fig_name = '{}{}.png'.format(figs_path, fig_name_short)
    if not os.path.isfile(os.path.abspath(fig_name)):
        logger.info(' Saved: {}'.format(fig_name))
        snowav.framework.figures.save_fig(f1, fig_name)

    fig_name_short = 'validation_locations'
    fig_name = '{}{}.png'.format(figs_path, fig_name_short)
    if not os.path.isfile(os.path.abspath(fig_name)):
        logger.info(' Saved: {}'.format(fig_name))
        snowav.framework.figures.save_fig(f3, fig_name)
