
from matplotlib import pyplot as plt
import matplotlib.dates as mdates
import seaborn as sns
import snowav.framework.figures
import numpy as np
from datetime import datetime

def diagnostics(args, logger):
    '''
    Simple model output diagnostics for the water year and report
    period, plots water year and report period mean precip [in], mean swe [in],
    snow line, and density.

    Args
    ----------
    args : dict
        dictionary with required inputs, see swi() figure for more information.

    logger : object

    '''
    if args['print']:
        print('diagnostics() figure args:\n')
        for name in args.keys():
            if name not in ['masks','image','swe']:
                print(name, ': ', args[name])

    basins = args['dbasins']
    depthlbl = args['depthlbl']
    elevlbl = args['elevlbl']
    start_date = args['start_date']
    barcolors = args['barcolors']
    barcolors.insert(0,'black')

    ix = args['snow_line'] == 0
    args['snow_line'][ix] = np.nan
    ix = args['density'] == 0
    args['density'][ix] = np.nan

    sns.set_style('darkgrid')
    sns.set_context('notebook')
    dfmt = mdates.DateFormatter('%b-%d')

    plt.close(0)
    f, a = plt.subplots(num = 0, figsize = (8,9), facecolor = 'white',
                        dpi = args['dpi'], nrows = 5, ncols = 2,
                        linewidth = 0.75)

    a = a.flatten()

    for iters, name in enumerate(basins):
        clr = barcolors[iters]

        # wy precip
        a[0].plot(args['precip'][name], label = name, color = clr, linewidth = 1)
        a[0].set_ylabel('precip [{}]'.format(depthlbl))

        # report period precip
        a[1].plot(args['precip_per'][name], label = name, color = clr, linewidth = 1)
        a[1].set_ylabel('precip [{}]'.format(depthlbl))
        m = args['precip_per'].max().max()
        if m < 0.25:
            m = 0.25
        a[1].set_ylim((-0.01, m + m*0.1))

        # wy SWE depth
        a[2].plot(args['swe'][name], label = name, color = clr, linewidth = 1)
        a[2].set_ylabel('swe [{}]'.format(depthlbl))

        # report period SWE depth
        a[3].plot(args['swe_per'][name], label = name, color = clr, linewidth = 1)
        a[3].set_ylabel('$\Delta$ swe [{}]'.format(depthlbl))

        # wy density
        a[4].plot(args['density'][name], label = name, color = clr, linewidth = 1)
        a[4].set_ylabel(r'density [$kg/m^3$]')
        a[4].set_ylim((100,600))

        # report period density
        a[5].plot(args['density_per'][name], label = name, color = clr, linewidth = 1)
        a[5].set_ylabel(r'$\Delta$ density [$kg/m^3$]')

        # wy snow line
        a[6].plot(args['snow_line'][name], label = name, color = clr, linewidth = 1)
        a[6].set_ylabel('snow line [{}]'.format(elevlbl))

        # report period snow line
        a[7].plot(args['snow_line_per'][name], label = name, color = clr, linewidth = 1)
        a[7].set_ylabel('$\Delta$ snow line [{}]'.format(elevlbl))

        # wy evap_z
        a[8].plot(args['evap_z'][name].cumsum(), label = name, color = clr, linewidth = 1)
        a[8].set_ylabel('evap_z [{}]'.format(depthlbl))

        # report period evap_z
        a[9].plot(args['evap_z_per'][name].cumsum(), label = name, color = clr, linewidth = 1)
        a[9].set_ylabel('$\Delta$ evap_z [{}]'.format(depthlbl))

    for n in range(0,len(a)):

        if args['flt_flag']:
            for i,d in enumerate(args['flight_dates']):
                if n == 0 and i == 0:
                    lb = 'flight'.format(args['wy'])
                else:
                    lb = '__nolabel__'

                a[n].axvline(x=d, linestyle=':', linewidth=0.75, color='r', label=lb)

        if n in [1,3,5,7,9]:
            a[n].yaxis.tick_right()
            a[n].yaxis.set_label_position('right')
            a[n].set_xlim((args['start_date'], args['end_date']))

        else:
            a[n].set_xlim((datetime(args['wy']-1, 10, 1), args['end_date']))

            if n == 0:
                lbl = 'report period'
            else:
                lbl = '__nolabel__'

            a[n].axvline(x=start_date, linestyle=':', linewidth=0.65, color='k',
                         label = lbl)

        if n < 8:
            a[n].set_xticklabels('')

        else:
            a[n].xaxis.set_major_formatter(dfmt)
            for tick in a[n].get_xticklabels():
                tick.set_rotation(45)

    a[0].set_title('Water Year')
    a[0].legend(loc='upper left', fontsize = 8)
    a[1].set_title('Report Period')
    plt.tight_layout()

    del barcolors[0]

    fig_name_short = 'diagnostics_'
    fig_name = '{}{}{}.png'.format(args['figs_path'],fig_name_short,args['directory'])
    if logger is not None:
        logger.info(' saving {}'.format(fig_name))
    snowav.framework.figures.save_fig(f, fig_name)

    return fig_name_short
