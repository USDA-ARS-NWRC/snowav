
"""
Test snowav processing.

"""

import unittest
import shutil
import os
import snowav
from snowav.database.database import collect
from scripts.snow import can_i_snowav

def compare_database_dataframes(run_name):
    """

    """
    # define arguments for collect()

    # df = collect( - test)
    # df = collect( - gold)

    # compare dataframes
    print(x)

    return True

def check_figure_creation():
    print(x)

# check_report_creation

# check_config_options

# check_utils_calculate

# check_process

# check_query



class TestStandardLakes(unittest.TestCase):
    """
    Test snowav processing using wy2019 Lakes basin.

    """

    @classmethod
    def setUpClass(self):
        """

        """

        # set up and run test snowav processing
        config_file = os.path.abspath('tests/lakes/snowav_lakes_test.ini')
        flag = can_i_snowav(config_file)

    def test_database_results(self):
    	"""
    	Query the test database for the gold standard results and compare to
        results for the current snowav run.

    	"""

    	a = compare_database_dataframes('test')
    	assert(a)

    def test_figure_creation(self):
    	"""
    	Query the test database for the gold standard results and compare to
        results for the current snowav run.

    	"""

    	a = check_figure_creation()
    	assert(a)

if __name__ == '__main__':
    unittest.main()
