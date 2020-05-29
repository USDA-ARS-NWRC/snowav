import argparse
import coloredlogs
import logging
import os

from snowav.framework.framework import Snowav
from snowav.framework.process_day import process
from snowav.framework.point_values import run_point_values
from snowav.plotting.swe_volume import swe_volume
from snowav.plotting.image_change import image_change
from snowav.plotting.plotlims import plotlims as plotlims


def run():
    parser = argparse.ArgumentParser(description="Process AWSM model "
                                                 "results.")

    parser.add_argument('-f', '--snowav_config', dest='snowav_config',
                        type=str, help='Path to snowav configuration file.')

    parser.add_argument('--end_date', '--end-date', dest='end_date', type=str,
                        help='End date that will override what is given in the'
                             ' config file, intended for airflow application.')

    parser.add_argument('-t', '--topo_path', dest='topo_path', type=str,
                        help='Path to topo.nc file.')

    parser.add_argument('-d', '--nc_path', dest='nc_path', type=str,
                        help='Path to single snow.nc file for swe '
                             'volume calculation and figure without database '
                             'processing.')

    parser.add_argument('-v', '--value', dest='value', type=str,
                        default='swe_vol',
                        choices=['swe_z', 'swe_vol', 'swe_avail', 'swe_unavail'],
                        help='SWE value to calculate for single snow.nc file.')

    parser.add_argument('-A', '--snow_a', dest='snow_a', type=str,
                        help='Path to "A" snow.nc file for simple difference, '
                             'must also include "-B", and must not use "-d".')

    parser.add_argument('-B', '--snow_b', dest='snow_b', type=str,
                        help='Path to "B" snow.nc file for simple difference, '
                             'must also include "-A", and must not use "-d".')

    parser.add_argument('-p', '--figs_path', dest='figs_path', type=str,
                        default=os.getcwd() + '/',
                        help='Path to save the figure, default is ./snowav/')

    parser.add_argument('-pc', '--point_values_config',
                        dest='point_values_config', type=str,
                        help='Path to point values config.ini file.')

    parser.add_argument('-pm', '--point_values_master',
                        dest='point_values_master', type=str,
                        help='Path to PointValuesCoreConfig.ini master config '
                             'file.')

    parser.add_argument('-pb', '--blank', dest='blank', action='store_true',
                        default=False, help='Make a blank example csv file '
                                            'for [basin] locations_csv')

    args = parser.parse_args()
    snowav_main(config_file=args.snowav_config, topo_path=args.topo_path,
                nc_path=args.nc_path, value=args.value, snow_a=args.snow_a,
                snow_b=args.snow_b, figs_path=args.figs_path,
                point_values_config=args.point_values_config,
                point_values_master=args.point_values_master,
                blank=args.blank)


def snowav_main(config_file=None, topo_path=None, nc_path=None, value=None,
                snow_a=None, snow_b=None, figs_path=None,
                point_values_config=None, point_values_master=None,
                blank=None):
    """ Command line function for snowav.

    Runs the standard snowav processing of AWSM files if config_file is
    present.

    """

    # standard snowav processing
    if config_file is not None:
        if not os.path.isfile(config_file):
            raise Exception('Invalid config file {} '.format(config_file))
        else:
            end_date = None

        Snowav(config_file=config_file, end_date=end_date)

    # point values processing
    point_values_args = [point_values_config, blank]

    if sum([x is not None for x in point_values_args]) > 0:
        run_point_values(point_values_config, point_values_master, blank)

    # single snow.nc file
    single_args = [topo_path, figs_path]
    if sum([x is not None for x in single_args]) > 1:
        barcolors = ['xkcd:cobalt',
                     'xkcd:mustard green',
                     'xkcd:lichen',
                     'xkcd:pale green',
                     'xkcd:blue green',
                     'xkcd:bluish purple',
                     'xkcd:lightish purple',
                     'xkcd:deep magenta',
                     'xkcd:burgundy',
                     'red']
        diff = False

        if not os.path.isdir(figs_path):
            raise Exception('{} not a valid directory'.format(figs_path))

        if not os.path.isfile(topo_path):
            raise Exception('{} not a valid topo file'.format(topo_path))

        # single snow.nc
        if nc_path is not None:
            nc_path = [os.path.abspath(nc_path)]
            if not os.path.isfile(nc_path[0]):
                raise Exception('{} not a valid file name'.format(nc_path[0]))

        # difference between two snow.nc
        if snow_a is not None and snow_b is not None:
            diff = True
            nc_path = [os.path.abspath(snow_a)] + [os.path.abspath(snow_b)]

            for path in nc_path:
                if not os.path.isfile(path):
                    raise Exception('{} not a valid file name'.format(path))

        log = create_log()
        results = process(nc_path, os.path.abspath(topo_path), value, log)

        # difference between two snow.nc files
        if diff:
            if len(nc_path[0]) > 20:
                title = nc_path[1][-20::] + ' - ' + nc_path[0][-20::]
            else:
                title = nc_path

            df = results['df'][nc_path[1]] - results['df'][nc_path[0]]
            image = (results['outputs']['swe_z'][1] - results['outputs']['swe_z'][0]) * 0.03937
            fig_name = os.path.join(os.path.abspath(figs_path), 'volume_change.png')

            log.info(' Difference generated from:\n      {} '
                     'subtract\n      {}'.format(nc_path[1], nc_path[0]))

            image_change(results['masks'], image, df, results['plotorder'],
                         plotlims(results['plotorder']), results['edges'],
                         results['labels'], barcolors, [0.5, 99.5],
                         'TAF', 'ft', 'in', (0, len(results['edges'])), [10, 5],
                         title, 1, os.path.abspath(figs_path), fig_name)

            log.info(' Saved: {}'.format(os.path.abspath(figs_path)))

        else:
            if len(nc_path[0]) > 35:
                title = '...' + nc_path[0][-20::]
            else:
                title = nc_path

            fig_name = os.path.join(os.path.abspath(figs_path), 'swe_volume.png')
            image = results['outputs']['swe_z'][0] / 25.4
            title = '{}, {}'.format(value, title)
            df = results['df'][os.path.abspath(nc_path[0])]

            swe_ylims = swe_volume(results['masks'], image, df,
                                   results['plotorder'],
                                   plotlims(results['plotorder']),
                                   results['edges'], results['labels'],
                                   barcolors, [0.5, 99.5], title, 'in', 'TAF',
                                   'ft', (0, len(results['edges'])), 1,
                                   os.path.abspath(figs_path), fig_name,
                                   [10, 5], dpi=200, show=True)

            log.info(' Saved: {}'.format(os.path.join(figs_path, fig_name)))


def create_log(log_level='INFO'):
    """ Create logger. """

    level_styles = {'info': {'color': 'white'},
                    'notice': {'color': 'magenta'},
                    'verbose': {'color': 'blue'},
                    'success': {'color': 'green', 'bold': True},
                    'spam': {'color': 'green', 'faint': True},
                    'critical': {'color': 'red', 'bold': True},
                    'error': {'color': 'red'},
                    'debug': {'color': 'green'},
                    'warning': {'color': 'yellow'}}

    field_styles = {'hostname': {'color': 'magenta'},
                    'programname': {'color': 'cyan'},
                    'name': {'color': 'white'},
                    'levelname': {'color': 'white', 'bold': True},
                    'asctime': {'color': 'green'}}

    log_level = log_level.upper()

    numeric_level = getattr(logging, log_level, None)
    if not isinstance(numeric_level, int):
        raise ValueError('Invalid log level: %s' % log_level)

    fmt = '%(levelname)s:%(name)s:%(message)s'
    logging.basicConfig(level=numeric_level)
    coloredlogs.install(level=numeric_level,
                        fmt=fmt,
                        level_styles=level_styles,
                        field_styles=field_styles)

    log = logging.getLogger(__name__)

    return log


def can_i_snowav(config_file):
    """ Function that wraps snowav for test case.

    Args
    ------
    config_file : str
        path to config file
    """

    success = True
    try:
        Snowav(config_file=config_file)
    except Exception as e:
        print(e)
        success = False

    return success


if __name__ == '__main__':
    run()
