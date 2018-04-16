
import SNOWAV

config_file = '/mnt/volumes/wkspace/config/snowav/snowav_brb_wy2018_ops.txt'
snow = SNOWAV.snowav(config_file)
    
# Make all the calculations
SNOWAV.snowav.process(snow)

# Figures
SNOWAV.snowav.accumulated(snow,'sub')
SNOWAV.snowav.current_image(snow)
SNOWAV.snowav.state_by_elev(snow)
SNOWAV.snowav.image_change(snow)
SNOWAV.snowav.basin_total(snow)
SNOWAV.snowav.pixel_swe(snow)
# SNOWAV.snowav.density(snow)
# SNOWAV.snowav.stn_validate(snow,31+19) # 31+19 or 0

SNOWAV.report(snow)

# Summaries for previous years
# SNOWAV.snowav.write_summary(snow, 'state_summary')
# SNOWAV.snowav.write_summary(snow, 'accum_summary')



