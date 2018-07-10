
import sys
import snowav

def run(config_file):

    snow = snowav.framework.framework.SNOWAV(config_file = config_file)
    snow.process()

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

    for name in snow.summary:
        snowav.plotting.write_summary.write_summary(snow,name)

    if snow.flt_flag is True:
        snowav.plotting.flt_image_change.flt_image_change(snow)

    if snow.report_flag is True:
        snowav.report.report.report(snow)

if __name__ == '__main__':
    config_file = sys.argv[1]
    run(config_file)
