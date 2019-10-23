
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
import cmocean
import snowav.framework.figures
import matplotlib.gridspec as gridspec
from mpl_toolkits.axes_grid1.inset_locator import inset_axes

def point_values(output, value, df, imgxy, filename, dem, figs_path, veg_type,
                 heading, model_date, course_date, settings, pixel, logger):
    '''
    Read in specified csv file of snow course validation locations, write out
    model values at those points, and create figures.

    Args
    ------
    output : array
        array of output values
    value : str
    df : DataFrame
        from point_values_csv file:
        name        latitude    longitude
        Bishop Pass 37.1        -118.557
        Bench Lake  36.958      -118.445

    imgxy : array
    x and y coordinates for the image, generated from get_topo_stats.py
    filename : str
    dem : array
    figs_path : str
    veg_type : array
    heading : str
    model_date : str
    course_date : str
    settings : list
    pixel : int
    logger

    '''

    if value == 'swe_z':
        header = cbarlabel = 'swe [in]'

    if value == 'density':
        header = 'density [kg/m^3]'
        cbarlabel = r'density [$kg/m^3$]'

    if value == 'depth':
        header = cbarlabel = 'depth [in]'

    if heading is not None:
        cols = ['station','x','y','model x','model y', header,
                '9 pixel median', heading]
    else:
        cols = ['station','x','y','model x','model y', header,
                '9 pixel median']

    stns = df['name']
    pixel_swe = pd.DataFrame(columns = cols)
    # ps = pd.DataFrame(columns= cols)

    # get closest model pixel
    for i in df.index:
        sta = df.loc[i,'name']
        lat = df.loc[i,'latitude']
        lon = df.loc[i,'longitude']

        ll = utm.from_latlon(lat,lon)
        xind = np.where(abs(imgxy[0]-ll[0]) == min(abs(imgxy[0]-ll[0])))[0]
        yind = np.where(abs(imgxy[1]-ll[1]) == min(abs(imgxy[1]-ll[1])))[0]

        if xind in [0, len(output[0,:])-1] or yind in [0, len(output[:,0])-1]:
            logger.warn('Latitude and/or longitude for {} in [validate] '
                        'point_values_csv are outside of the model domain, '
                        'exiting point_values'.format(sta,filename))
            return

        pixel_swe.loc[i,'station'] = sta
        pixel_swe.loc[i,'x'] = lat
        pixel_swe.loc[i,'y'] = lon
        pixel_swe.loc[i,'model x'] = int(xind)
        pixel_swe.loc[i,'model y'] = int(yind)

        if heading is not None:
            pixel_swe.loc[i,heading] = df.loc[i,heading]

        # Get the median of the nine closest pixels
        val = np.array([])
        for n in range(-1,1):
            for m in range(-1,1):
                val = np.append(val,output[yind+n,xind+m])
                if n == 0 and m == 0:
                    pixel_swe.loc[i,header] = np.round(float(output[yind+n,xind+m]),decimals = 2)

        pixel_swe.loc[i,'9 pixel median'] = np.round(np.nanmedian(val),decimals = 2)

    pixel_swe.to_csv(filename)

    ###########################################################################
    #   figures
    ###########################################################################

    if heading is not None:

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
        f, a = plt.subplots(num=0,figsize=(int(width),int(height)),
                            dpi=dpi,nrows=nrows,ncols=ncols)
        a = a.flatten()

        plt.close(1)
        f1, a1 = plt.subplots(num=1,figsize=(int(width),int(height)),
                            dpi=dpi,nrows=nrows,ncols=ncols)
        a1 = a1.flatten()

        plt.close(2)
        f2 = plt.figure(num=2,figsize=(width,height-2),dpi=dpi)
        a2 = plt.gca()

        plt.close(3)
        f3 = plt.figure(num=3,figsize=(width,height-2),dpi=dpi)
        a3 = plt.gca()

        lvls = list(np.arange(np.min(dem.flatten()), np.max(dem.flatten()), 500))

        a3.imshow(veg_type, cmap = plt.cm.get_cmap('Vega20'), alpha = 0.75)
        a3.contour(dem, colors = 'k', levels=lvls, linewidths = 0.15)

        vegmap = plt.cm.get_cmap('brg', len(np.unique(veg_type)))

        swe_values = pixel_swe[header].values

        cvals = plt.cm.YlGn
        snow_map = cmocean.cm.haline_r

        # get min/max for all sub-domains
        pixel_swe = pd.DataFrame(pixel_swe)
        veg = []
        for idx in pixel_swe.index:
            name = pixel_swe.loc[idx,'station']
            ix = int(pixel_swe.loc[idx]['model x'])
            iy = int(pixel_swe.loc[idx]['model y'])

            veg_sub = veg_type[(iy-npix):(iy+npix+1),(ix-npix):(ix+npix+1)].flatten()
            veg = np.append(veg_sub, veg)

            if idx == 0:
                g_min = np.nanmin(output[(iy-npix):(iy+npix),(ix-npix):(ix+npix)])
                g_max = np.nanmax(output[(iy-npix):(iy+npix),(ix-npix):(ix+npix)])

            else:
                vmin = np.nanmin(output[(iy-npix):(iy+npix),(ix-npix):(ix+npix)])
                vmax = np.nanmax(output[(iy-npix):(iy+npix),(ix-npix):(ix+npix)])

                if vmin < g_min:
                    g_min = vmin

                if vmax > g_max:
                    g_max = vmax

        g_min = 0
        veg = np.unique(veg)

        place = np.linspace(g_min,g_max,256)
        stns = stns.tolist() + [' ']

        for idx, v in enumerate(pixel_swe.station.unique()):
            mflag = False

            if idx < len(pixel_swe.station.unique()) + 1:

                # this plots repeating values
                if len(pixel_swe[pixel_swe['station'] == v]['station'].values) > 1:
                    mflag = True
                    vstr = mstr = ''

                    for idx3 in range(0,len(pixel_swe[pixel_swe['station'] == v]['station'].values)):
                        course = pixel_swe.loc[idx + idx3][heading]
                        ix = int(pixel_swe[pixel_swe['station'] == v]['model x'].values[idx3])
                        iy = int(pixel_swe[pixel_swe['station'] == v]['model y'].values[idx3])
                        ixs = int(np.where(abs(place-course) == min(abs(place-course)))[0])
                        color = snow_map(ixs)

                        if idx3 == 0:
                            vstr = '{}'.format(round(course,1))
                            mstr = '{}'.format(output[iy,ix].round(1))
                        else:
                            vstr = vstr + ', {}'.format(round(course,1))
                            mstr = mstr + ', {}'.format(output[iy,ix].round(1))

                        a[idx].scatter(x = ix, y = iy, marker = '^', s=ss, c = color,
                                       edgecolors = 'k', linewidths = 0.5)

                    cstr = 'validation: {}\nmodel: {}'.format(vstr, mstr)
                    a[idx].text(annot_x1, annot_y1, cstr, horizontalalignment='left',
                                transform=a[idx].transAxes,fontsize = font_small - 1,
                                color = 'w')

                ix = int(pixel_swe[pixel_swe['station'] == v]['model x'].values[0])
                iy = int(pixel_swe[pixel_swe['station'] == v]['model y'].values[0])
                swe = pixel_swe[pixel_swe['station'] == v][header].values[0]
                course = pixel_swe[pixel_swe['station'] == v][heading].values[0]

                s = dem[(iy-npix):(iy+npix+1),(ix-npix):(ix+npix+1)]
                snow_sec = output[(iy-npix):(iy+npix+1),(ix-npix):(ix+npix+1)]
                l1 = [int(x)*snow_sec.shape[0]+npix-1 for x in range(npix-1,npix+2)] + \
                     [int(x)*snow_sec.shape[0]+npix for x in range(npix-1,npix+2)] + \
                     [int(x)*snow_sec.shape[0]+npix+1 for x in range(npix-1,npix+2)]

                if idx == 0:
                    ml = 'model, pixel containing point'
                    vl = 'validation'
                    sl = r'model, $\pm${} pixels'.format(npix)

                else:
                    ml = '__nolabel__'
                    vl = '__nolabel__'
                    sl = '__nolabel__'

                for idx2, v2 in enumerate(list(snow_sec.flatten())):
                    clr = 'b'
                    mkr = '.'

                    if idx == 0 and idx2 == 0:
                        sl = r'model, $\pm${} pixels'.format(npix)
                    else:
                        sl = '__nolabel__'

                    if idx2 in l1:
                        clr = 'g'
                        mkr = 'X'
                        mkrs = 10

                        if idx == 1 and idx2 == l1[-1]:
                            sl = 'model, $\pm$1 pixel'

                    a2.plot(idx, v2, '{}{}'.format(clr, mkr), markersize = 7, label = sl)

                    sl = '__nolabel__'

                a2.plot(idx, course, 'kP', markersize = 10, label = vl)
                a2.plot(idx, swe, 'rX', markersize = 10, label = ml)

                lvls = list(np.arange(np.min(s.flatten()), np.max(s.flatten()), levels))
                ixs = int(np.where(abs(place-course) == min(abs(place-course)))[0])

                color = snow_map(ixs)

                a[idx].contour(dem, colors = 'k', levels=lvls, linewidths = 0.25)
                h = a[idx].imshow(output, cmap = snow_map, clim = (g_min, g_max))

                if not mflag:
                    a[idx].scatter(x = ix, y = iy, marker = '^', s=ss, c = color,
                                   edgecolors = 'k', linewidths = 0.5)
                    cstr = 'validation: {}\nmodel: {}'.format(round(course,1),
                                                             output[iy,ix].round(1))
                    a[idx].text(annot_x1, annot_y1,cstr,horizontalalignment='left',
                                transform=a[idx].transAxes,fontsize = font_small,
                                color = 'w')

                a[idx].get_xaxis().set_ticks([])
                a[idx].get_yaxis().set_ticks([])
                a[idx].set_ylim((iy - (npix + 0.5), iy + npix + 0.5))
                a[idx].set_xlim((ix - (npix + 0.5), ix + npix + 0.5))
                a[idx].set_title(v, fontsize = font_small)

                # legend of sorts
                veg_subp = veg_type[(iy-npix):(iy+npix+1),(ix-npix):(ix+npix+1)]
                vegmap = plt.cm.get_cmap('Vega20')

                a1[idx].contour(dem, colors = 'k', levels=lvls, linewidths = 0.25)
                h1 = a1[idx].imshow(veg_type, cmap = vegmap, clim =
                                    (min(np.unique(veg_type[(iy-npix):(iy+npix+1),(ix-npix):(ix+npix+1)].flatten())),
                                     max(np.unique(veg_type[(iy-npix):(iy+npix+1),(ix-npix):(ix+npix+1)].flatten()))))
                a1[idx].scatter(x = ix, y = iy, marker = '^', s=ss, c = 'k',
                               edgecolors = 'k', linewidths = 0.5)
                a1[idx].get_xaxis().set_ticks([])
                a1[idx].get_yaxis().set_ticks([])
                a1[idx].set_ylim((iy - (npix + 0.5), iy + npix + 0.5))
                a1[idx].set_xlim((ix - (npix + 0.5), ix + npix + 0.5))
                a1[idx].set_title(v, fontsize = font_small)

                a3.scatter(x = ix, y = iy, marker = '*', s=ss-10, c = 'k',
                           edgecolors = 'k', linewidths = 0.5)
                a3.annotate(v, (ix+5, iy-5), fontsize = font_small-1, color = 'k')

                for vs in np.unique(veg_subp):
                    ixf = np.where(veg_subp == vs)
                    fx = ixf[0][0]
                    fy = ixf[1][0]
                    nx = ix + fy - npix
                    ny = iy + fx - npix
                    a1[idx].scatter(x = nx, y = ny, marker = 'X', s=3, c = 'k',
                                   edgecolors = 'k', linewidths = 0.5)
                    a1[idx].annotate(str(vs), (nx, ny), fontsize = font_small - 3, color = 'w')


        # tack on another subplot for legend and colorbar
        a[idx+1].imshow(output, cmap = snow_map, alpha = 0)
        a[idx+1].set_title('', fontsize = font_small)
        a[idx+1].get_xaxis().set_ticks([])
        a[idx+1].get_yaxis().set_ticks([])
        a[idx+1].axis('off')
        cbaxes = inset_axes(a[idx+1], width="90%", height="10%", loc=8)
        cbar = f.colorbar(h, cax=cbaxes, ticks = [0, int(g_max/2), int(g_max)],
                          orientation='horizontal')
        cbar.ax.tick_params(labelsize=font_medium)
        cbar.set_label(cbarlabel, fontsize = font_medium)
        cbar.ax.xaxis.set_label_position('top')
        cstr = ('Values in pixel containing coordinates\n'
              'model: {}\nvalidation: {}\ncontour interval: {} '
              'ft\npixel width: {} m, {} ft'.format(model_date, course_date,
              levels, str(pixel), str(int(pixel*3.28))))
        a[idx+1].text(annot_x2, annot_y2,cstr,horizontalalignment='left',
                    transform=a[idx+1].transAxes,fontsize = font_medium)

        a1[idx+1].imshow(veg_type, cmap = vegmap, clim = (np.min(veg),np.max(veg)), alpha = 0)
        a1[idx+1].set_title('', fontsize = font_small)
        a1[idx+1].get_xaxis().set_ticks([])
        a1[idx+1].get_yaxis().set_ticks([])
        a1[idx+1].axis('off')
        cstr = ('model: {}\nvalidation: {}\ncontour interval: {} '
              'ft\npixel width: {} m, {} ft'.format(model_date, course_date,
              levels, str(pixel), str(int(pixel*3.28))))
        a1[idx+1].text(annot_x2, annot_y2,cstr,horizontalalignment='left',
                    transform=a1[idx+1].transAxes,fontsize = font_medium)

        for n in range(idx+2,nrows*ncols):
            f.delaxes(a[n])

        for n in range(idx+2,nrows*ncols):
            f1.delaxes(a1[n])

        a2.set_ylabel(cbarlabel)
        a2.set_xticks(list(range(0,len(stns))))
        a2.set_xticklabels(pixel_swe.station.unique(), rotation = 90, fontsize = font_small - 1)
        a2.set_xlim((-0.5,len(pixel_swe.station.unique())))
        a2.legend(loc = 2)
        a2.grid(linewidth = 0.25)
        a2.set_title('Validation and Model Pixel Values, {}'.format(model_date))

        a3.get_xaxis().set_ticks([])
        a3.get_yaxis().set_ticks([])

        f.tight_layout()
        f1.tight_layout()
        f2.tight_layout()
        f3.tight_layout()

        fig_name_short = 'validation_list_{}'.format(value)
        fig_name = '{}{}.png'.format(figs_path,fig_name_short)
        logger.info(' saving {}'.format(fig_name))
        snowav.framework.figures.save_fig(f2, fig_name)

        fig_name_short = 'validation_map_{}'.format(value)
        fig_name = '{}{}.png'.format(figs_path,fig_name_short)
        logger.info(' saving {}'.format(fig_name))
        snowav.framework.figures.save_fig(f, fig_name)

        fig_name_short = 'validation_veg_map'
        fig_name = '{}{}.png'.format(figs_path,fig_name_short)
        logger.info(' saving {}'.format(fig_name))
        snowav.framework.figures.save_fig(f1, fig_name)

        fig_name_short = 'validation_locations'
        fig_name = '{}{}.png'.format(figs_path,fig_name_short)
        logger.info(' saving {}'.format(fig_name))
        snowav.framework.figures.save_fig(f3, fig_name)


        return
