
import seaborn as sns
from matplotlib import pyplot as plt
import snowav.framework.figures
import matplotlib.dates as mdates

def inputs(args, logger):
    '''
    Simple summary plots of smrf inputs.

    Args
    ----------
    args : dict
        dictionary with required inputs, see swi() figure for more information.

    logger : list
        snowav logger

    '''
    dfmt = mdates.DateFormatter('%b-%d')
    sns.set_style('whitegrid')

    # variables = ['air_temp','precip']
    values = args['inputs']
    barcolors = args['barcolors']
    lines = ['-',':','--','-.','-',':','--','-.']
    # print(values)

    plt.close(0)
    f, a = plt.subplots(num=0,figsize=(8,6),dpi=args['dpi'],nrows=3,ncols=3)
    a = a.flatten()

    for idx, val in enumerate(values.keys()):
        # air_temp
        for ixn, name in enumerate(values[val].keys()):
            # Cherry
            for ixs, stat in enumerate(values[val][name]):
                # mean

                clr = barcolors[ixn]
                lstyle = lines[ixs]
                vlstyle = lines[ixs + 1]


                if idx == 0 and stat == [*values[val][name].keys()][0]:
                    lbl = name
                elif idx == 1 and name == [*values[val].keys()][0]:
                    lbl = stat
                else:
                    lbl = '__nolabel__'

                if 'percentile' in stat:
                    x = values[val][name][stat]['value'].index
                    y1 = values[val][name]['nanpercentile_25']['value']
                    y2 = values[val][name]['nanpercentile_75']['value']
                    a[idx].fill_between(x, y1, y2, color=clr, label = lbl, alpha='0.1')

                else:
                    a[idx].plot(values[val][name][stat]['value'],
                                label = lbl,
                                color = clr,
                                linestyle = lstyle,
                                linewidth = 0.75)

                a[idx].set_title(val, fontsize = 6)
                a[idx].xaxis.set_major_formatter(dfmt)
                for tick in a[idx].get_xticklabels():
                    tick.set_rotation(45)
                    tick.set_fontsize(6)

                for tick in a[idx].get_yticklabels():
                    tick.set_fontsize(6)

    # a[0].legend(loc=2, fontsize=6)
    # a[1].legend(loc=2, fontsize=6)
    plt.tight_layout()

    fig_name_short = 'inputs_'
    fig_name = '{}{}{}.png'.format(args['figs_path'],fig_name_short,args['directory'])
    if logger is not None:
        logger.info(' saving {}'.format(fig_name))
    snowav.framework.figures.save_fig(f, fig_name)

    plt.show()
