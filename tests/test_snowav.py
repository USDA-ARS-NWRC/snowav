
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

def check_gold_results():
    ''' '''
    connector = 'sqlite:///'+os.path.abspath('./tests/lakes/results/results.db')
    plotorder = ['Lakes Basin']
    basins = {'Lakes Basin': {'watershed_id': 1, 'basin_id': 1}}
    start_date = datetime(2019,4,1,23,0,0)
    end_date = datetime(2019,4,2,23,0,0)
    value = 'swe_vol'
    run_name_gold = 'gold'
    edges = [8000,9000,10000,11000]
    gold_values = [3.075, 13.63, 9.543, 0.494]

    gold = collect(connector, plotorder, basins, start_date, end_date, value,
                 run_name_gold, edges, 'end')

    result = True
    for ix, edge in enumerate(edges):
        if (gold.iloc[ix,0] - gold_values[ix]) != 0.0:
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


    # swe = './tests/lakes/results/lakes_test_20190401_20190402/swe_volume_lakes_tes.png'
    # swi = './tests/lakes/results/lakes_test_20190401_20190402/swi_lakes_test.png'
    # cold = './tests/lakes/results/lakes_test_20190401_20190402/cold_content_lakes_test.png'
    # swe_change = './tests/lakes/results/lakes_test_20190401_20190402/swe_change_lakes_test.png'
    # swe_volume = './tests/lakes/results/lakes_test_20190401_20190402/swe_volume_lakes_test.png'
    # precip_depth = './tests/lakes/results/lakes_test_20190401_20190402/precip_depth_lakes_test.png'
    # diagnostics = './tests/lakes/results/lakes_test_20190401_20190402/diagnostics_lakes_test.png'
    # report = './tests/lakes/results/lakes_test_20190401_20190402/SnowpackSummary20190403.pdf'

def check_swi_figure():
    """ Simple check if .png figures were created. """
    result = True
    fig = './tests/lakes/results/lakes_test_20190401_20190402/swi_lakes_test.png'
    if not os.path.isfile(os.path.abspath(fig)):
        result = False

    return result

def check_swe_figure():
    """ Simple check if .png figures were created. """
    result = True
    fig = './tests/lakes/results/lakes_test_20190401_20190402/swe_volume_lakes_test.png'
    if not os.path.isfile(os.path.abspath(fig)):
        result = False

    return result

def check_cold_content_figure():
    """ Simple check if .png figures were created. """
    result = True
    fig = './tests/lakes/results/lakes_test_20190401_20190402/cold_content_lakes_test.png'
    if not os.path.isfile(os.path.abspath(fig)):
        result = False

    return result

def check_inputs_figure():
    """ Simple check if .png figures were created. """
    result = True
    fig = './tests/lakes/results/lakes_test_20190401_20190402/inputs_lakes_test.png'
    fig2 = './tests/lakes/results/lakes_test_20190401_20190402/inputs_period_lakes_test.png'

    if not os.path.isfile(os.path.abspath(fig)):
        result = False

    if not os.path.isfile(os.path.abspath(fig2)):
        result = False

    return result

def check_swe_change_figure():
    """ Simple check if .png figures were created. """
    result = True
    fig = './tests/lakes/results/lakes_test_20190401_20190402/swe_change_lakes_test.png'
    if not os.path.isfile(os.path.abspath(fig)):
        result = False

    return result

def check_diagnostics_figure():
    """ Simple check if .png figures were created. """
    result = True
    fig = './tests/lakes/results/lakes_test_20190401_20190402/diagnostics_lakes_test.png'
    if not os.path.isfile(os.path.abspath(fig)):
        result = False

    return result

def check_precip_depth_figure():
    """ Simple check if .png figures were created. """
    result = True
    fig = './tests/lakes/results/lakes_test_20190401_20190402/precip_depth_lakes_test.png'
    if not os.path.isfile(os.path.abspath(fig)):
        result = False

    return result

def check_report():
    """ Simple check if .png figures were created. """
    result = True
    fig = './tests/lakes/results/lakes_test_20190401_20190402/SnowpackSummary20190403.pdf'
    if not os.path.isfile(os.path.abspath(fig)):
        result = False

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

    def test_gold_results(self):
    	''' Check that gold results are on database '''

    	a = check_gold_results()
    	assert(a)

    def test_database_swe_vol(self):
    	''' Standard and current swe_vol DataFrames '''

    	a = compare_database_swe_vol()
    	assert(a)

    def test_database_swi_vol(self):
    	""" Standard and current swi_vol DataFrames """

    	a = compare_database_swi_vol()
    	assert(a)

    def test_swe_figure(self):
    	""" Output swe figure .png"""

    	a = check_swe_figure()
    	assert(a)

    def test_inputs_figure(self):
    	""" Output swe figure .png"""

    	a = check_inputs_figure()
    	assert(a)

    def test_swe_change_figure(self):
    	""" Output swe change figure .png"""

    	a = check_swe_change_figure()
    	assert(a)

    def test_swi_figure(self):
    	""" Output swi figure .png"""

    	a = check_swi_figure()
    	assert(a)

    def test_cold_content_figure(self):
    	""" Output cold content figure .png"""

    	a = check_cold_content_figure()
    	assert(a)

    def test_diagnostics_figure(self):
    	""" Output diagnostics figure .png"""

    	a = check_diagnostics_figure()
    	assert(a)

    def test_report(self):
    	""" Output report .pdf """

    	a = check_report()
    	assert(a)

    @classmethod
    def tearDownClass(self):
        """ Remove figures and report """
        basin_total = './tests/lakes/results/lakes_test_20190401_20190402/basin_total_lakes_test.png'
        swe_change = './tests/lakes/results/lakes_test_20190401_20190402/swe_change_lakes_test.png'
        swe_vol = './tests/lakes/results/lakes_test_20190401_20190402/swe_volume_lakes_test.png'
        swi_vol = './tests/lakes/results/lakes_test_20190401_20190402/swi_lakes_test.png'
        precip_depth = './tests/lakes/results/lakes_test_20190401_20190402/precip_depth_lakes_test.png'
        cold_content = './tests/lakes/results/lakes_test_20190401_20190402/cold_content_lakes_test.png'
        diagnostics = './tests/lakes/results/lakes_test_20190401_20190402/diagnostics_lakes_test.png'
        inputs = './tests/lakes/results/lakes_test_20190401_20190402/inputs_lakes_test.png'
        inputs_period = './tests/lakes/results/lakes_test_20190401_20190402/inputs_period_lakes_test.png'
        report = './tests/lakes/results/lakes_test_20190401_20190402/SnowpackSummary20190403.pdf'
        config = './tests/lakes/results/lakes_test_20190401_20190402/lakes_test_20190401_20190402.ini'

        os.remove(os.path.abspath(basin_total))
        os.remove(os.path.abspath(swe_vol))
        os.remove(os.path.abspath(swi_vol))
        os.remove(os.path.abspath(swe_change))
        os.remove(os.path.abspath(precip_depth))
        os.remove(os.path.abspath(cold_content))
        os.remove(os.path.abspath(diagnostics))
        os.remove(os.path.abspath(inputs))
        os.remove(os.path.abspath(inputs_period))
        os.remove(os.path.abspath(report))
        os.remove(os.path.abspath(config))        

if __name__ == '__main__':
    unittest.main()
