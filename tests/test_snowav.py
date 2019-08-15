
"""
Initial unittest for snowav processing, using wy2019 results in the Lakes Basin.

tests/ includes 'gold' data in tests/lakes/results/results.db, with
run_name: gold on the database. The end-to-end snowav run generated here places
run_name: test results and compares.

These tests in the Lakes should be expanded, including tests for:
    config options and parsing
    elevation bin calculations
    report creation
    query
    flights

Tests could also be expanded by including a Tuolumne dataset and multi-day
snow.nc files.

"""

import unittest
import shutil
import os
from snowav.database.database import collect
from snowav.utils.utilities import calculate, masks
from scripts.snow import can_i_snowav
from datetime import datetime
from snowav.utils.OutputReader import iSnobalReader

def check_utils_masks():
    ''' Test mask creation from topo.nc file. '''

    result = True
    dem = os.path.abspath('./tests/lakes/topo/topo.nc')
    out = masks(dem, False)

    if out['nrows'] != 168:
        result = False

    if out['ncols'] != 156:
        result = False

    if out['plotorder'] != ['Lakes Basin']:
        result = False

    return result

def check_utils_calculate():
    ''' Test snowav volume and mean depth calculation. '''

    result = True
    path = os.path.abspath('./tests/lakes/gold/runs/run20190402/')
    dem = os.path.abspath('./tests/lakes/topo/topo.nc')
    filetype = 'netcdf'
    wy = 2019
    pixel = 50
    out = masks(dem, False)
    mask = out['masks']['Lakes Basin']['mask']

    output = iSnobalReader(path, filetype, snowbands = [0,1,2],
                           embands = [6,7,8,9], wy = wy)

    array = output.snow_data[0][0,:,:]
    out = calculate(array, pixel, mask, 'mean', 'snow_depth')

    if out != 116.411:
        result = False

    array = output.snow_data[2][0,:,:]
    out = calculate(array, pixel, mask, 'sum', 'volume')

    if out != 26.742:
        result = False

    return result

def compare_database_swe_vol():
    ''' Compare 'gold' and 'test' swe_vol. '''

    connector = 'sqlite:///'+os.path.abspath('./tests/lakes/results/results.db')
    plotorder = ['Lakes Basin']
    basins = {'Lakes Basin': {'watershed_id': 1, 'basin_id': 1}}
    start_date = datetime(2019,4,1,23,0,0)
    end_date = datetime(2019,4,2,23,0,0)
    value = 'swe_vol'
    run_name_gold = 'gold'
    run_name_test = 'test'
    edges = [8000,9000,10000,11000]

    gold = collect(connector, plotorder, basins, start_date, end_date, value,
                 run_name_gold, edges, 'end')
    test = collect(connector, plotorder, basins, start_date, end_date, value,
                 run_name_test, edges, 'end')

    test_result = gold - test

    if test_result.sum().any() != 0:
        result = False

    else:
        result = True

    return result

def compare_database_swi_vol():
    ''' Compare 'gold' and 'test' swi_vol. '''

    connector = 'sqlite:///'+os.path.abspath('./tests/lakes/results/results.db')
    plotorder = ['Lakes Basin']
    basins = {'Lakes Basin': {'watershed_id': 1, 'basin_id': 1}}
    start_date = datetime(2019,4,1,23,0,0)
    end_date = datetime(2019,4,2,23,0,0)
    value = 'swi_vol'
    run_name_gold = 'gold'
    run_name_test = 'test'
    edges = [8000,9000,10000,11000]

    gold = collect(connector, plotorder, basins, start_date, end_date, value,
                 run_name_gold, edges, 'sum')
    test = collect(connector, plotorder, basins, start_date, end_date, value,
                 run_name_test, edges, 'sum')

    test_result = gold - test

    if test_result.sum().any() != 0:
        result = False

    else:
        result = True

    return result

def check_figure_creation():
    """ Simple check if .png figures were created. """

    result = True
    swe = './tests/lakes/results/lakes_test_20190401_20190402/swe_volume_lakes_test.png'
    swi = './tests/lakes/results/lakes_test_20190401_20190402/swi_lakes_test.png'

    if not os.path.isfile(os.path.abspath(swe)):
        result = False
    else:
        os.remove(os.path.abspath(swe))

    if not os.path.isfile(os.path.abspath(swi)):
        result = False
    else:
        os.remove(os.path.abspath(swi))

    return result


class TestStandardLakes(unittest.TestCase):
    ''' Test snowav processing using wy2019 Lakes basin for comparison.'''

    @classmethod
    def setUpClass(self):
        ''' End-to-end snowav run '''

        config_file = os.path.abspath('tests/lakes/snowav_lakes_test.ini')
        flag = can_i_snowav(config_file)

    def test_utils_masks(self):
        ''' Utils masks '''

        a = check_utils_masks()
        assert(a)

    def test_utils_calculate(self):
    	''' Utils calculate '''

    	a = check_utils_calculate()
    	assert(a)

    def test_database_swe_vol(self):
    	''' Standard and current swe_vol DataFrames '''

    	a = compare_database_swe_vol()
    	assert(a)

    def test_database_swi_vol(self):
    	""" Standard and current swi_vol DataFrames """

    	a = compare_database_swi_vol()
    	assert(a)

    def test_figure_creation(self):
    	""" Output png files for swe and swi figures"""

    	a = check_figure_creation()
    	assert(a)

if __name__ == '__main__':
    unittest.main()
