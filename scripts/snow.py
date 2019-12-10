
from sys import exit
import argparse
import os
import coloredlogs
import logging
import pandas as pd
from snowav.framework.framework import snowav
from snowav.framework import framework
from snowav.framework.process_day import process
from snowav.plotting.swe_volume import swe_volume
from snowav.plotting.image_change import image_change
from snowav.plotting.plotlims import plotlims as plotlims
from snowav.utils.utilities import masks

def run():

    parser = argparse.ArgumentParser(description='Process AWSM model results, '
                                     'put results on database, create figures, '
                                     'and make pdf report. For other features, '
                                     'options, and behavior see README.md '
                                     'and CoreConfig.ini.')

    parser.add_argument('-f', '--config_file', dest='config_file', type=str,
                        help='Path to snowav configuration file.')

    parser.add_argument('-end_date', '--end-date', dest='end_date', type=str,
                        help='End date that will override what is given in the '
                        'config file, intended for airflow application.')

    parser.add_argument('-t', '--topo_path', dest='topo_path', type=str,
                        help='Path to topo.nc file.')

    parser.add_argument('-d', '--nc_path', dest='nc_path', type=str,
                        help='Path to single snow.nc file for swe '
                        'volume calculation and figure without database '
                        'processing.')

    parser.add_argument('-v', '--value', dest='value', type=str,
                        default='swe_vol',
                        choices=['swe_z','swe_vol','swe_avail','swe_unavail'],
                        help='SWE value to calculate for single snow.nc file.')

    parser.add_argument('-s', dest='show', action='store_true',
                        help='Flag to display single snow.nc calculation.')

    parser.add_argument('-no-s', dest='show', action='store_false')

    parser.set_defaults(show=False)

    parser.add_argument('-A', '--snow_a', dest='snow_a', type=str,
                        help='Path to "A" snow.nc file for simple difference, '
                        'must also include "-B", and must not use "-d".')

    parser.add_argument('-B', '--snow_b', dest='snow_b', type=str,
                        help='Path to "B" snow.nc file for simple difference, '
                        'must also include "-A", and must not use "-d".')

    parser.add_argument('-p', '--figs_path', dest='figs_path', type=str,
                        default=os.getcwd() + '/',
                        help='Path to save the figure, default is ./snowav/')

    args = parser.parse_args()

    ###########################################################################
    #  standard snowav processing                                             #
    ###########################################################################

    if args.config_file is not None:
        if not os.path.isfile(args.config_file):
            raise Exception('Config file {} does not '
                            'exist'.format(args.config_file))

        if args.end_date is not None:
            try:
                end_date = pd.to_datetime(args.end_date)
            except:
                raise Exception('pandas failed parsing {} to datetime'.format(
                    args.end_date))

        else:
            end_date = None

        snowav(config_file = args.config_file, end_date = end_date)
        exit()

    ###########################################################################
    #  single snow.nc file, or simple difference, no database                 #
    ###########################################################################
    if args.topo_path is None:
        raise Exception('Must use -t <topo.nc> argument for command line tools')

    if not os.path.isdir(args.figs_path):
        raise Exception('{} not a valid directory'.format(args.figs_path))

    if not os.path.isfile(args.topo_path):
        raise Exception('{} not a valid topo file'.format(args.topo_path))

    # Must use either single snow.nc, or difference between two snow.nc
    if args.nc_path is not None:
        nc_path = [os.path.abspath(args.nc_path)]
        diff = False
        if not os.path.isfile(nc_path[0]):
            raise Exception('{} not a valid file name'.format(nc_path[0]))

    elif (args.snow_a is not None) and (args.snow_b is not None):
        nc_path = [os.path.abspath(args.snow_a)] + [os.path.abspath(args.snow_b)]
        diff = True

        for path in nc_path:
            if not os.path.isfile(path):
                raise Exception('{} not a valid file name'.format(path))

    else:
        parser.error('Must use either -d, or both -A and -B arguments')

    log = logging.getLogger(__name__)
    level = 'INFO'
    levelname = 'INFO'
    coloredlogs.install(fmt='%(levelname)-5s %(message)s', level=level,
                        level_styles={'info':{'color':'green'}},logger=log)

    results = process(nc_path, os.path.abspath(args.topo_path), args.value, log)

    barcolors = ['xkcd:cobalt', 'xkcd:mustard green', 'xkcd:lichen',
                 'xkcd:pale green', 'xkcd:blue green', 'xkcd:bluish purple',
                 'xkcd:lightish purple', 'xkcd:deep magenta']

    fargs = {'show': args.show,
             'print':False,
             'masks':results['masks'],
             'barcolors':barcolors,
             'plotorder':results['plotorder'],
             'lims':plotlims(results['plotorder']),
             'vollbl':'TAF',
             'elevlbl':'ft',
             'depthlbl':'in',
             'edges':results['edges'],
             'labels':results['labels'],
             'percent_max':99.5,
             'percent_min':0.05,
             'figsize': (10,5),
             'dpi':250,
             'dplcs':1,
             'xlims':(0,len(results['edges'])),
             'figs_path':os.path.abspath(args.figs_path),
             'directory':'cli'}

    # difference between two snow.nc files
    if diff:
        if len(nc_path[0]) > 20:
            title = nc_path[1][-20::] + ' - ' + nc_path[0][-20::]
        else:
            title = nc_path

        df = results['df'][nc_path[1]] - results['df'][nc_path[0]]
        fargs['df'] = df
        fargs['image'] = results['outputs']['swe_z'][1] \
                         - results['outputs']['swe_z'][0]
        fargs['title'] = title

        log.info('Difference generated from:\n      {} '
                 'subtract\n      {}'.format(nc_path[1],nc_path[0]))
        log.info('Saved figure in {}'.format(os.path.abspath(args.figs_path)))

        image_change(fargs, None)
        exit()

    # snow.nc swe volume
    if not diff:
        if len(args.nc_path) > 35:
            title = '...' + args.nc_path[-20::]
        else:
            title = args.nc_path

        fargs['image'] = results['outputs']['swe_z'][0]/25.4
        fargs['title'] = '{}, {}'.format(args.value,title)
        fargs['df'] = results['df'][os.path.abspath(args.nc_path)]

        swe_volume(fargs, None)
        log.info('Saved figure in {}'.format(args.figs_path))

def can_i_snowav(config_file):
    '''
    Function that wraps snowav run in try/except for test case.

    Args
    ------
    config_file : str
        path to config file
    '''

    try:
        snowav(config_file = config_file)
        return True

    except Exception as e:
        print(e)
        return False

if __name__ == '__main__':
    run()
