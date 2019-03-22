
from snowav.framework.process_day import process_day
import argparse
from snowav.plotting.swe_volume import swe_volume


def main():


    parser = argparse.ArgumentParser(description="Examine and auto populate")

    parser.add_argument('-d', '--nc_path', dest='nc_path', type=str,
                        help='Path to snow.nc file.')

    parser.add_argument('-b', '--basin', dest='basin', type=str,
                        help='Basin.')

    args = parser.parse_args()

    # do a check on args.basin for in tables

    process_day(args.nc_path, args.basin)

    swe_volume()

if __name__ == '__main__':
    main()
