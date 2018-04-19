

import sys
sys.path.append('/home/markrobertson/mrworkspace/code/SNOWAV/')
import snowav
from snowav.report import report
from snowav.plotting.accumulated import accumulated
from snowav.plotting.image_change import image_change
from snowav.plotting.current_image import current_image
from snowav.plotting.state_by_elev import state_by_elev
from snowav.plotting.basin_total import basin_total
from snowav.plotting.stn_validate import stn_validate

config_file = '/mnt/volumes/wkspace/config/snowav/snowav_brb_wy2018_ops.txt'
  
snow = snowav.plotting.framework.SNOWAV(config_file) 

# Make all the calculations
snow.process()

# Plots
accumulated(snow)
current_image(snow)
state_by_elev(snow)
image_change(snow)
basin_total(snow)
stn_validate(snow) 

# Report
report.report(snow)




