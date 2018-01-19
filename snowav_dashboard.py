
import sys
sys.path.append('/home/markrobertson/mrworkspace/code/SNOWAV/')
import SNOWAV

config_file     = '/home/markrobertson/mrworkspace/code/SNOWAV/config/snowav_tuol_wy2018.txt'
snow            = SNOWAV.snowav(config_file)
    
# Make all the calculations
SNOWAV.snowav.process(snow)

# To make a report for the BRB, or with no flight updates, run these
SNOWAV.snowav.accumulated(snow,'sub')
SNOWAV.snowav.current_image(snow)
SNOWAV.snowav.state_by_elev(snow)
SNOWAV.snowav.image_change(snow)
SNOWAV.snowav.basin_total(snow)
SNOWAV.snowav.stn_validate(snow) # just for BRB right now
SNOWAV.report(snow)

# To make a report with flight updates, add this figure
#     [Report]  -> tex_file = tuol_report_flt.tex
#     [Outputs] -> psnowFile_flt and csnowFile_flt
SNOWAV.snowav.process(snow,snow.psnowFile_flt,snow.csnowFile_flt)
SNOWAV.snowav.image_change(snow,'_flt')
SNOWAV.report(snow,'_update')

# Summaries for previous years
# SNOWAV.snowav.write_summary(snow, 'state_byelev')

