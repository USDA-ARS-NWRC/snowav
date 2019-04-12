
from snowav.framework.process_day import day
import argparse
from snowav.plotting.swe_volume import swe_volume
from snowav.database.tables import Watershed, Basin, Watersheds, Basins
from sys import exit


def main():
    '''
    python scripts/day.py -b SJ -d /data/blizzard/sanjoaquin/ops/wy2019/ops/runs/run20190301/snow.nc -p /home/markrobertson/wkspace/


    '''

    parser = argparse.ArgumentParser(description="Examine and auto populate")

    parser.add_argument('-d', '--nc_path', dest='nc_path', type=str,
                        help='Path to snow.nc file.')

    parser.add_argument('-b', '--basin', dest='basin', type=str,
                        help='Basin.')

    parser.add_argument('-p', '--figs_path', dest='figs_path', type=str,
                        help='Path to save the figure.')

    args = parser.parse_args()

    if args.basin not in Watersheds.watersheds.keys():
        print('Given basin "{}", must be in {}'.format(args.basin,basins))
        exit()

    day.process_day(day, args.nc_path, basin=args.basin)

    if args.figs_path is not None:
        day.figs_path = args.figs_path

    else:
        day.figs_path = './'

    swe_volume(day=day)

if __name__ == '__main__':
    main()
