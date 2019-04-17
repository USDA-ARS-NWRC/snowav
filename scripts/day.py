
from snowav.framework.process_day import day
import argparse
from snowav.plotting.swe_volume import swe_volume
from snowav.plotting.swe_difference import swe_difference
from snowav.database.tables import Watershed, Basin, Watersheds, Basins
from sys import exit


def main():
    '''

    Default processing settings in snowav/database/tables.py

    SWE volume by elevation for snow.nc file:
    $ snowav_day -b "San Joaquin River Basin" -d snow.nc -p wkspace/

    Difference in volume for two snow.nc files:
    $ snowav_day -b "San Joaquin River Basin" -A snow.nc -B snow.nc -p ~/wkspace/


    '''

    parser = argparse.ArgumentParser(description="Examine and auto populate")

    parser.add_argument('-d', '--nc_path', dest='nc_path', type=str,
                        help='Path to snow.nc file.')

    parser.add_argument('-v', '--value', dest='value', type=str,
                        help='Value to plot, options are: coldcont, evap_z, '
                             'rain_z, density, depth, precip_z, precip_vol, '
                             'swi_z, swi_vol, swe_z, swe_vol, swe_avail, '
                             'swe_unavail. Default is swe_z.')

    parser.add_argument('-A', '--snow_a', dest='snow_a', type=str,
                        help='Path to "A" snow.nc file.')

    parser.add_argument('-B', '--snow_b', dest='snow_b', type=str,
                        help='Path to "B" snow.nc file.')

    parser.add_argument('-b', '--basin', dest='basin', type=str,
                        help='Basin.')

    parser.add_argument('-p', '--figs_path', dest='figs_path', type=str,
                        help='Path to save the figure.')

    args = parser.parse_args()

    if args.basin not in Watersheds.watersheds.keys():
        print('Given basin "{}", must be one of '
              '{}'.format(args.basin,Watersheds.watersheds.keys()))
        exit()

    var_options = ['coldcont','evap_z','rain_z','density','depth',
            'precip_z','precip_vol','swi_z','swi_vol','swe_z',
            'swe_vol','swe_avail','swe_unavail']

    if args.value is None:
        args.value = 'swe_z'

    if args.value not in var_options:
        print('Given value "{}", must be one of '
              '{}'.format(args.value,var_option))
        exit()

    if (args.snow_a is not None) and (args.snow_b is not None):
        nc_path = [args.snow_a] + [args.snow_b]
    else:
        nc_path = args.nc_path

    day.process_day(day, nc_path, args.value, basin=args.basin)

    if args.figs_path is not None:
        day.figs_path = args.figs_path

    else:
        day.figs_path = './'

    if (args.snow_a is not None) and (args.snow_b is not None):
        swe_difference(day=day)

    elif args.value in ['swe_vol','swe_avail','swe_unavail','precip_vol','swi_vol']:
        swe_volume(day=day)

if __name__ == '__main__':
    main()
