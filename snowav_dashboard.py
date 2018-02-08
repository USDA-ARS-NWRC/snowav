
import sys
sys.path.append('/home/markrobertson/mrworkspace/code/SNOWAV/')
import SNOWAV

config_file     = '/home/markrobertson/mrworkspace/code/SNOWAV/config/snowav_sj_wy2018_ops.txt'
snow            = SNOWAV.snowav(config_file)
    
# Make all the calculations
SNOWAV.snowav.process(snow)

# To make a report for the BRB, or with no flight updates, run these
SNOWAV.snowav.accumulated(snow)
SNOWAV.snowav.current_image(snow)
SNOWAV.snowav.state_by_elev(snow)
SNOWAV.snowav.image_change(snow)
SNOWAV.snowav.basin_total(snow)
rundirs = ['/mnt/data/snowdrift/brb/ops/wy2018/runs/run20171001_20180107/',
           '/mnt/data/snowdrift/brb/ops/wy2018/runs/run20180108_20180117/',
           '/mnt/data/snowdrift/brb/ops/wy2018/ops/runs/run20180117_20180204/']
# SNOWAV.snowav.stn_validate(snow,rundirs) # just for BRB right now
# SNOWAV.snowav.basin_detail(snow)
SNOWAV.report(snow)

# To make a report with flight updates, add this figure
#     [Report]  -> tex_file = tuol_report_flt.tex
#     [Outputs] -> psnowFile_flt and csnowFile_flt
SNOWAV.snowav.process(snow,snow.psnowFile_flt,snow.csnowFile_flt)
SNOWAV.snowav.image_change(snow,'_flt')
SNOWAV.report(snow,'_debug')

# Summaries for previous years
# SNOWAV.snowav.write_summary(snow, 'state_summary')

