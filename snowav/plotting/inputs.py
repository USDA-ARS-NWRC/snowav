
import seaborn as sns
from matplotlib import pyplot as plt
import snowav.framework.figures
import matplotlib.dates as mdates
from datetime import datetime

def inputs(args, logger):
    '''
    Simple summary of smrf inputs.

    Args
    ----------
    args : dict
        dictionary with required inputs, see swi() figure for more information.

    logger : list
        snowav logger

    '''

    dfmt = mdates.DateFormatter('%b-%d')
    sns.set_style('whitegrid')

    values = args['inputs']
    barcolors = args['barcolors']
    lines = ['-',':','--','-.','-',':','--','-.','-',':','--','-.']

    plt.close(0)
    f, a = plt.subplots(num=0,figsize=(8,6),dpi=args['dpi'],nrows=3,ncols=3)
    a = a.flatten()

    p = []
    for v in values[args['var_list'][0]][args['inputs_basins'][0]]:
        if 'percentile' in v:
            p.append(v)

    for idx, val in enumerate(args['var_list']):

        for ixn, name in enumerate(args['inputs_basins']):

            for ixs, stat in enumerate(args['inputs_methods']):

                clr = barcolors[ixn]
                lstyle = lines[ixs]
                vlstyle = lines[ixs + 1]

                lbl = '__nolabel__'

                if idx == 0 and stat == args['inputs_methods'][0]:
                    lbl = name

                if idx == 1 and name == args['inputs_basins'][0]:
                    lbl = stat

                if 'percentile' not in stat:
                    if idx == 0:
                        a[idx].plot(values[val][name][stat]['value'],
                                    label = lbl,
                                    color = clr,
                                    linestyle = lstyle,
                                    linewidth = 0.5)

                    else:
                        a[idx].plot(values[val][name][stat]['value'],
                                    label = lbl,
                                    color = clr,
                                    linestyle = lstyle,
                                    linewidth = 0.5)

                else:
                    x = values[val][name][stat]['value'].index
                    y1 = values[val][name][p[0]]['value']
                    y2 = values[val][name][p[1]]['value']

                    if not y1.empty or not y2.empty:
                        a[idx].fill_between(x, y1, y2, color=clr, alpha='0.075')

                a[idx].set_xlim((datetime(args['wy']-1, 10, 1), args['end_date']))
                a[idx].set_title(val, fontsize = 6)
                a[idx].xaxis.set_major_formatter(dfmt)
                for tick in a[idx].get_xticklabels():
                    tick.set_rotation(45)
                    tick.set_fontsize(6)

                for tick in a[idx].get_yticklabels():
                    tick.set_fontsize(6)

    for n in range(idx+1,len(a)):
        f.delaxes(a[n])

    a[0].legend(loc=3,fontsize=6)
    a[1].legend(loc=3,fontsize=6)

    plt.tight_layout()

    fig_name_short = 'inputs_'
    fig_name = '{}{}{}.png'.format(args['figs_path'],fig_name_short,
                                   args['directory'])
    if logger is not None:
        logger.info(' saving {}'.format(fig_name))
    snowav.framework.figures.save_fig(f, fig_name)

    for idx, val in enumerate(values.keys()):
        a[idx].set_xlim((args['start_date'], args['end_date']))

    fig_name_short = 'inputs_period_'
    fig_name = '{}{}{}.png'.format(args['figs_path'],fig_name_short,
                                   args['directory'])
    if logger is not None:
        logger.info(' saving {}'.format(fig_name))
    snowav.framework.figures.save_fig(f, fig_name)
