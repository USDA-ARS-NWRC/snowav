
import sys
sys.path.append('/home/markrobertson/mrworkspace/code/SNOWAV/')
import SNOWAV

config_file     = '/home/markrobertson/mrworkspace/code/SNOWAV/config/snowav_brb_wy2017.txt'
snow            = SNOWAV.snowav(config_file)
    
# Make all the calculations
SNOWAV.snowav.process(snow)

# Plots
SNOWAV.snowav.current_image(snow)
SNOWAV.snowav.accumulated(snow)
SNOWAV.snowav.image_change(snow)

# Report
SNOWAV.report(snow)