

import sys
import snowav

def run(config_file):

    # Initialize snowav, read in config file
    # config_file = '/mnt/volumes/wkspace/config/snowav/snowav_tuol_wy2018_devel.txt'
    snow = snowav.plotting.framework.SNOWAV(config_file = config_file)

    if not hasattr(snow,'error'):

        # Make all the calculations
        snow.process()

        # Save processed variables to netcdf if desired
        if snow.nc_flag == True:
            snowav.methods.output_nc.output_nc(snow)

        # Plots
        snowav.plotting.accumulated.accumulated(snow)
        snowav.plotting.current_image.current_image(snow)
        snowav.plotting.state_by_elev.state_by_elev(snow)
        snowav.plotting.image_change.image_change(snow)
        snowav.plotting.swe_change.swe_change(snow)
        snowav.plotting.basin_total.basin_total(snow)
        snowav.plotting.pixel_swe.pixel_swe(snow)
        snowav.plotting.density.density(snow)
        snowav.plotting.water_balance.water_balance(snow)
        snowav.plotting.stn_validate.stn_validate(snow)
        
        # If options exist in config file
        if hasattr(snow,'flt_flag'):
            snowav.plotting.flt_image_change.flt_image_change(snow)
        
        # snowav.plotting.basin_detail.basin_detail(snow)
        snowav.plotting.write_summary.write_summary(snow,'accum_summary')
        snowav.plotting.write_summary.write_summary(snow,'state_summary')

        # Generate report if desired
        if snow.report_flag == True:
            snowav.report.report.report(snow)

if __name__ == '__main__':
    config_file = sys.argv[1]
    run(config_file)
