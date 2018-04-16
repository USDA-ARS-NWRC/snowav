

import SNOWAV
config_file = '/mnt/volumes/wkspace/config/snowav/snowav_brb_wy2018_ops.txt'
snow = SNOWAV.snowav(config_file)
    
# Make all the calculations
SNOWAV.snowav.process(snow)

# Figures
SNOWAV.accumulated(snow)
SNOWAV.current_image(snow)
SNOWAV.state_by_elev(snow)
SNOWAV.image_change(snow)
SNOWAV.basin_total(snow)
SNOWAV.pixel_swe(snow)
SNOWAV.density(snow)
SNOWAV.stn_validate(snow) 

# Report
SNOWAV.report(snow)

# Summaries for previous years
# SNOWAV.snowav.write_summary(snow, 'state_summary')
# SNOWAV.snowav.write_summary(snow, 'accum_summary')



