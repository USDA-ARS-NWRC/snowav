
import sys
sys.path.append('/home/markrobertson/mrworkspace/code/SNOWAV/')
import SNOWAV

config_file     = '/home/markrobertson/mrworkspace/code/SNOWAV/config/snowav_tuol_wy2018.txt'
snow            = SNOWAV.snowav(config_file)
    
# Make all the calculations
SNOWAV.snowav.process(snow)

# Plots
#SNOWAV.snowav.accumulated(snow)
#SNOWAV.snowav.current_image(snow)
#SNOWAV.snowav.state_by_elev(snow)
#SNOWAV.snowav.image_change(snow)
#SNOWAV.snowav.basin_total(snow)
#SNOWAV.snowav.stn_validate(snow) # just for BRB right now
#SNOWAV.snowav.write_summary(snow, 'state_byelev')

# Flight comparison example
#SNOWAV.snowav.process(snow,snow.psnowFile_flt,snow.csnowFile_flt)
#SNOWAV.snowav.image_change(snow,'_flt')

# Report
SNOWAV.report(snow,'_update')