
from snowav.framework.process_day import Day
import argparse
from snowav.plotting.swe_volume import swe_volume
from snowav.plotting.swe_difference import swe_difference
from snowav.database.tables import Watershed, Basin, Watersheds, Basins
from sys import exit
from snowav.plotting.plotlims import plotlims as plotlims
from datetime import datetime
from sys import exit
import os

def main():
    '''
    Command line program to process and plot either SWE volume by elevation
    band for a single snow.nc file, or the difference between two snow.nc files.

    '''

    parser = argparse.ArgumentParser(description="Simple processing and "
                                     "plotting tool for single snow.nc file, "
                                     "or difference between two snow.nc files.")

    parser.add_argument('-d', '--nc_path', dest='nc_path', type=str,
                        help='Path to snow.nc file for simple swe volume '
                        'calculation and figure.')

    parser.add_argument('-wy', '--wy', dest='wy', type=int,
                        help='Water year, will be assigned as current if '
                        'empty.')

    parser.add_argument('-v', '--value', dest='value', type=str,
                        default='swe_z',
                        choices=['swe_z','swe_vol','swe_avail','swe_unavail'],
                        help='Value to plot.')

    parser.add_argument('-s', dest='show', action='store_true',
                        help='Optional, boolean to show figures, '
                        'default is False.')

    parser.add_argument('-no-s', dest='show', action='store_false')

    parser.set_defaults(show=False)

    parser.add_argument('-A', '--snow_a', dest='snow_a', type=str,
                        help='Path to "A" snow.nc file for simple difference, '
                        'must also include "-B", and must not use "-d".')

    parser.add_argument('-B', '--snow_b', dest='snow_b', type=str,
                        help='Path to "B" snow.nc file for simple difference, '
                        'must also include "-A", and must not use "-d".')

    parser.add_argument('-b', '--basin', dest='basin', type=str,
                        choices=['sanjoaquin','kings','lakes',
                                 'tuolumne','kaweah', 'merced'],
                        help='Basin, options are: sanjoaquin, kings, lakes, '
                        'tuolumne, brb, kaweah, merced.')

    parser.add_argument('-p', '--figs_path', dest='figs_path', type=str,
                        default=os.getcwd() + '/',
                        help='Path to save the figure, default is snowav/.')

    args = parser.parse_args()

    if not os.path.isdir(args.figs_path):
        raise OSError('{} not a valid directory'.format(args.figs_path))

    # Must use either single day, or difference between two days
    if args.nc_path is not None:
        nc_path = args.nc_path
        diff = False
        if not os.path.isfile(nc_path):
            raise OSError('{} not a valid file name'.format(nc_path))

    elif (args.snow_a is not None) and (args.snow_b is not None):
        nc_path = [args.snow_a] + [args.snow_b]
        diff = True

        for path in nc_path:
            if not os.path.isfile(path):
                raise OSError('{} not a valid file name'.format(path))

    else:
        parser.error('Must use either -d, or both -A and -B arguments')

    lims = plotlims(args.basin, plotorder=None)
    day = Day(lims.database_name,
              nc_path,
              args.value,
              args.figs_path,
              show=args.show)
    Day.process_day(day)

    if diff:
        swe_difference(day=day)
    else:
        swe_volume(day=day)

if __name__ == '__main__':
    main()
