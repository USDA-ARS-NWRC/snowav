
import sys
import snowav

'''
- basin_total() needs cleaned up; pull from database once all runs available;
    some dataframes are built on the order they are pulled from db rather than
    chronological order
- output csv option
- SWI versus streamflow (-> streamflow on db?)
- SWE v distance to reservoir?

'''

def run(config_file):

    snow = snowav.framework.framework.SNOWAV(config_file = config_file)

if __name__ == '__main__':
    config_file = sys.argv[1]
    run(config_file)
