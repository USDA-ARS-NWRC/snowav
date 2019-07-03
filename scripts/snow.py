
from sys import exit
import snowav
import argparse
import os
from snowav.framework.framework import snowav
from snowav.framework import framework
from snowav.framework.process_day import process
from snowav.plotting.swe_volume import swe_volume
from snowav.plotting.swe_difference import swe_difference
from snowav.plotting.plotlims import plotlims as plotlims
from snowav.utils.utilities import masks

def run():

    '''
    For the standard snowav run, pass a config.ini file:
        $ snowav -f <config.ini>


    Otherwise, must pass either:
        $ snowav -d <snow.nc> -t <topo.nc>
    or
        $ snowav -A <snow.nc> -B <snow.nc> -t <topo.nc>

    '''

    parser = argparse.ArgumentParser(description='Process AWSM model results, '
                                     'put results on database, create figures, '
                                     'and make pdf report. For other features, '
                                     'options, and behavior see CoreConfig.ini '
                                     'and README.md.')

    parser.add_argument('-f', '--config_file', dest='config_file', type=str,
                        help='Path to snowav configuration file for standard '
                        'processing runs.')

    parser.add_argument('-t', '--topo_path', dest='topo_path', type=str,
                        help='Path to topo.nc file associated with the snow.nc '
                        'file.')

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
    #                       standard snowav processing                        #
    ###########################################################################

    if args.config_file is not None:
        if not os.path.isfile(args.config_file):
            raise Exception('Config file {} does not '
                            'exist'.format(args.config_file))

        snowav(config_file = args.config_file)
        exit()

    ###########################################################################
    #                       single file                                       #
    ###########################################################################
    if args.topo_path is None:
        raise Exception('Must supply topo.nc file to use command line tools')

    if not os.path.isdir(args.figs_path):
        raise Exception('{} not a valid directory'.format(args.figs_path))

    if not os.path.isfile(args.topo_path):
        raise Exception('{} not a valid topo file'.format(args.topo_path))

    # Must use either single day, or difference between two days
    if args.nc_path is not None:
        nc_path = [args.nc_path]
        diff = False
        if not os.path.isfile(nc_path[0]):
            raise Exception('{} not a valid file name'.format(nc_path[0]))

    elif (args.snow_a is not None) and (args.snow_b is not None):
        nc_path = [args.snow_a] + [args.snow_b]
        diff = True

        for path in nc_path:
            if not os.path.isfile(path):
                raise Exception('{} not a valid file name'.format(path))

    else:
        parser.error('Must use either -d, or both -A and -B arguments')

    results = process(nc_path, args.topo_path, args.value)

    barcolors = ['xkcd:cobalt', 'xkcd:mustard green', 'xkcd:lichen',
                 'xkcd:pale green', 'xkcd:blue green', 'xkcd:bluish purple',
                 'xkcd:lightish purple', 'xkcd:deep magenta']

    fargs = {'print':False,
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
             'dplcs':2,
             'xlims':(0,len(results['edges'])),
             'figs_path':args.figs_path,
             'directory':''}

    if diff:
        swe_difference(day=day)
        print('Saved figure in ')
        exit()

    fargs['image'] = results['outputs']['swe_z'][0]
    fargs['title'] = 'title'
    fargs['df'] = results['df']

    swe_volume(fargs, None)

    print('Saved figure in {}'.format(args.figs_path))

if __name__ == '__main__':
    run()
