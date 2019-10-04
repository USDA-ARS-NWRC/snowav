
"""
Initial unittest for snowav processing, using wy2019 results in the Lakes Basin.

Tests:
- end-to-end snowav run
- simple topo.nc read and masks creation
- volume and depth calculations for SWE volume and mean depth
- database 'gold' and current values for swe_z, swe_vol, swe_unavail, precip_z,
    swi_z, swi_vol, and density
- standard figure .png creation
- report .pdf creation

Tests could be expanded by including a Tuolumne dataset and multi-day
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
import subprocess
import matplotlib
matplotlib.use('agg')

topo_path = './tests/lakes/topo/topo.nc'
db_path = './tests/lakes/results/gold.db'
tbl_path = './tests/lakes/gold/lakes_report_table.txt'
run_name_gold = 'lakes_wy2019_gold'
run_name_test = 'test'
start_date = datetime(2019,4,1,23,0,0)
end_date = datetime(2019,4,2,23,0,0)
plotorder = ['Lakes']
connector = 'sqlite:///'+os.path.abspath(db_path)
edges = [8000,9000,10000,11000]

# This basins dictionary matches the current snowav database.
# Note: this will differ from new sqlite databases created in >v0.10.0
basins = {'Lakes': {'watershed_id': 4, 'basin_id': 16}}

def check_utils_masks():
    ''' Test mask creation from topo.nc file. '''

    result = True
    dem = os.path.abspath(topo_path)
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
    dem = os.path.abspath(topo_path)
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
    ''' Check database 'gold' values '''

    value = 'swe_vol'
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

    value = 'swe_vol'

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

    value = 'swi_vol'

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

def compare_database_density():
    ''' Compare 'gold' and 'test' density. '''

    value = 'density'

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

def compare_database_precip():
    ''' Compare 'gold' and 'test' precip. '''

    value = 'precip_z'

    gold = collect(connector, plotorder, basins, start_date, end_date, value,
                 run_name_gold, 'total', 'sum')
    test = collect(connector, plotorder, basins, start_date, end_date, value,
                 run_name_test, 'total', 'sum')

    test_result = gold - test

    if test_result.sum().any() != 0:
        result = False
    else:
        result = True

    return result

def compare_database_unavail():
    ''' Compare 'gold' and 'test' swe_unavail. '''

    value = 'swe_unavail'

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

def check_density_figure():
    """ Simple check if .png figures were created. """
    result = True
    fig = './tests/lakes/results/lakes_test_20190401_20190402/density_lakes_test.png'
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

def check_cli_process():
    ''' Check command line snow.nc processing '''

    result = True
    cmd = 'snowav -t {} -d ./tests/lakes/gold/runs/run20190402/snow.nc -v swe_vol -no-s'.format(topo_path)
    process = subprocess.Popen(cmd.split(), stdout=subprocess.PIPE)
    output, error = process.communicate()

    if error is not None:
        result = False

    return result

class TestStandardLakes(unittest.TestCase):
    ''' Test snowav processing using wy2019 Lakes basin for comparison.'''

    @classmethod
    def setUpClass(self):
        ''' End-to-end snowav run '''

        config_file = os.path.abspath('tests/lakes/snowav_lakes_test.ini')
        flag = can_i_snowav(config_file)

    def test_cli_process(self):
        ''' Check command line process \n'''

        a = check_cli_process()
        assert(a)

    def test_utils_masks(self):
        ''' Check masks dictionary utility '''

        a = check_utils_masks()
        assert(a)

    def test_utils_calculate(self):
    	''' Check calculate utility '''

    	a = check_utils_calculate()
    	assert(a)

    def test_gold_results(self):
    	''' Check that gold results are on database '''

    	a = check_gold_results()
    	assert(a)

    def test_database_swe_vol(self):
    	''' Gold and current swe_vol DataFrames '''

    	a = compare_database_swe_vol()
    	assert(a)

    def test_database_swi_vol(self):
    	""" Gold and current swi_vol DataFrames """

    	a = compare_database_swi_vol()
    	assert(a)

    def test_database_density(self):
    	""" Gold and current density DataFrames """

    	a = compare_database_density()
    	assert(a)

    def test_database_precip(self):
    	""" Gold and current precip DataFrames """

    	a = compare_database_precip()
    	assert(a)

    def test_database_unavail(self):
    	""" Gold and current swe_unavail DataFrames """

    	a = compare_database_unavail()
    	assert(a)

    def test_density_figure(self):
    	""" Output density figure .png"""

    	a = check_density_figure()
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

        base = './tests/lakes/results/lakes_test_20190401_20190402/'

        basin_total = '{}basin_total_lakes_test.png'.format(base)
        swe_change = '{}swe_change_lakes_test.png'.format(base)
        swe_vol = '{}swe_volume_lakes_test.png'.format(base)
        swi_vol = '{}swi_lakes_test.png'.format(base)
        precip_depth = '{}precip_depth_lakes_test.png'.format(base)
        density = '{}density_lakes_test.png'.format(base)
        cold_content = '{}cold_content_lakes_test.png'.format(base)
        diagnostics = '{}diagnostics_lakes_test.png'.format(base)
        inputs = '{}inputs_lakes_test.png'.format(base)
        inputs_period = '{}inputs_period_lakes_test.png'.format(base)
        report = '{}SnowpackSummary20190403.pdf'.format(base)
        config = '{}lakes_test_20190401_20190402.ini'.format(base)

        os.remove(os.path.abspath(basin_total))
        os.remove(os.path.abspath(swe_vol))
        os.remove(os.path.abspath(swi_vol))
        os.remove(os.path.abspath(density))
        os.remove(os.path.abspath(swe_change))
        os.remove(os.path.abspath(precip_depth))
        os.remove(os.path.abspath(cold_content))
        os.remove(os.path.abspath(diagnostics))
        os.remove(os.path.abspath(inputs))
        os.remove(os.path.abspath(inputs_period))
        os.remove(os.path.abspath(report))
        os.remove(os.path.abspath(config))

        if os.path.isfile(os.path.abspath('./swe_volume_cli.png')):
            os.remove(os.path.abspath('./swe_volume_cli.png'))

        if os.path.isfile('{}swe_vol_test_TAF.csv'.format(base)):
            os.remove('{}swe_vol_test_TAF.csv'.format(base))

        if os.path.isfile('{}swi_vol_test_TAF.csv'.format(base)):
            os.remove('{}swi_vol_test_TAF.csv'.format(base))

if __name__ == '__main__':
    unittest.main()
