
from snowav.framework.process_day import process_day
import argparse
from snowav.plotting.swe_volume import swe_volume


def main():


    parser = argparse.ArgumentParser(description="Examine and auto populate")

    parser.add_argument('-d', '--nc_path', dest='nc_path', type=str,
                        help='Path to snow.nc file.')

    parser.add_argument('-b', '--basin', dest='basin', type=str,
                        help='Basin.')

    parser.add_argument('-p', '--figs_path', dest='figs_path', type=str,
                        help='Path to save the figure.')

    args = parser.parse_args()
    print(args)
    print(args.nc_path)
    # do a check on args.basin for in tables

    day = process_day(args.nc_path, args.basin)

    if args.figs_path is not None:
        day.figs_path = args.figs_path
        
    else:
        day.figs_path = '.'

    swe_volume(day=day)

if __name__ == '__main__':
    main()
