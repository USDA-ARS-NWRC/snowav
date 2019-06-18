

import sys
import snowav
import argparse
import os
from snowav.framework.framework import snowav

def run():

    parser = argparse.ArgumentParser(description='Process AWSM model results, '
                                     'put results on database, create figures, '
                                     'and make pdf report. SNOWAV expects '
                                     'results in awsm_daily format and topo.nc '
                                     'files created with basin_setup.')

    parser.add_argument('-f', '--config_file', dest='config_file', type=str,
                        help='Path to snowav configuration file.')

    args = parser.parse_args()

    if os.path.isfile(args.config_file):
        snowav(config_file = args.config_file)

    else:
        raise IOError('Config file {} does not exist'.format(args.config_file))

if __name__ == '__main__':
    run()
