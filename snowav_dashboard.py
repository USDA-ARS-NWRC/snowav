
import sys
sys.path.append('/home/markrobertson/mrworkspace/code/SNOWAV/')
import SNOWAV

config_file     = '/home/markrobertson/mrworkspace/code/SNOWAV/config/snowav_brb_wy2017.txt'
snow            = SNOWAV.snowav(config_file)

SNOWAV.snowav.calc(snow)
SNOWAV.snowav.current_image(snow)
SNOWAV.snowav.accumulated(snow)