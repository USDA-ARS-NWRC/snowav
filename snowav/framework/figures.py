
import numpy as np
from datetime import datetime
from snowav.plotting.swi import swi
from snowav.plotting.basin_total import basin_total
from snowav.plotting.cold_content import cold_content
from snowav.plotting.compare_runs import compare_runs
from snowav.plotting.current_image import current_image
from snowav.plotting.density import density
from snowav.plotting.flt_image_change import flt_image_change
from snowav.plotting.image_change import image_change
from snowav.plotting.precip_depth import precip_depth
from snowav.plotting.stn_validate import stn_validate
from snowav.plotting.swe_change import swe_change
from snowav.plotting.swe_volume import swe_volume
from snowav.plotting.plotlims import plotlims
from snowav.database.database import collect, make_session
from snowav.plotting.plotlims import plotlims as plotlims

def figures(self):
    '''
    Set up and call snowav figures. See CoreConfig.ini and README.md for more on
    config options and use.

    Note: swe_volume() must be called before cold_content() if you want to use
    the same ylims for each.

    '''

    args = {'report_start':self.report_start.date().strftime("%Y-%-m-%-d"),
            'report_date':self.report_date.date().strftime("%Y-%-m-%-d"),
            'print':self.print_args_dict,
            'run_name':self.run_name,
            'start_date':self.start_date,
            'end_date':self.end_date,
            'directory':self.directory,
            'figs_path':self.figs_path,
            'edges':self.edges,
            'plotorder':self.plotorder,
            'labels':self.labels,
            'lims':plotlims(self.plotorder),
            'masks':self.masks,
            'figsize':self.figsize,
            'dpi':self.dpi,
            'depthlbl':self.depthlbl,
            'vollbl':self.vollbl,
            'elevlbl':self.elevlbl,
            'dplcs':self.dplcs,
            'barcolors':self.barcolors,
            'xlims':self.xlims,
            'depth_clip':self.depth_clip,
            'percent_min':self.percent_min,
            'percent_max':self.percent_max,
            'basins':self.basins}

    fig_names = {}
    connector = self.connector

    # For each figure, first collect 2D array image, by-elevation
    # DataFrame, and set figure-specific inputs

    if self.swi_flag:
        image = np.zeros_like(self.outputs['swi_z'][0])
        for n in range(self.ixs,self.ixe):
            image = image + self.outputs['swi_z'][n]*self.depth_factor

        df = collect(connector, args['plotorder'], args['basins'],
                     args['start_date'], args['end_date'], 'swi_vol',
                     args['run_name'], args['edges'], 'sum')

        args['df'] = df
        args['image'] = image
        args['title'] = 'Accumulated SWI\n{} to {}'.format(
                        self.start_date.date().strftime("%Y-%-m-%-d"),
                        self.end_date.date().strftime("%Y-%-m-%-d"))

        fig_names['swi'] = swi(args, self._logger)

    if self.image_change_flag:
        image = self.outputs['swe_z'][self.ixs] - self.outputs['swe_z'][self.ixe]

        start = collect(connector, args['plotorder'], args['basins'],
                        args['start_date'], args['start_date'],'swe_vol',
                        args['run_name'], args['edges'],'end')
        end = collect(connector, args['plotorder'], args['basins'],
                      args['start_date'], args['end_date'], 'swe_vol',
                      args['run_name'], args['edges'],'end')

        df = end - start

        args['df'] = df
        args['image'] = image*self.depth_factor
        args['title'] = 'Change in SWE Depth\n{} to {}'.format(
                        args['report_start'], args['report_date'])

        fig_names['image_change'] = image_change(args, self._logger)

    if self.swe_volume_flag:
        df = collect(connector, args['plotorder'], args['basins'],
                     args['start_date'], args['end_date'],'swe_vol',
                     args['run_name'], args['edges'],'end')

        image = self.outputs['swe_z'][self.ixe]*self.depth_factor

        args['df'] = df
        args['image'] = image
        args['title'] = 'SWE {}'.format(args['report_date'])

        fig_names['swe_volume'], args['ylims'] = swe_volume(args, self._logger)

    if self.cold_content_flag:
        df = collect(connector, args['plotorder'], args['basins'],
                     args['start_date'], args['end_date'], 'swe_unavail',
                     args['run_name'], args['edges'],'end')

        swe = self.outputs['swe_z'][self.ixe]
        image = self.outputs['coldcont'][self.ixe]*0.000001

        args['df'] = df
        args['swe'] = swe
        args['image'] = image
        args['title'] = 'Cold Content {}'.format(args['report_date'])

        fig_names['cold_content'] = cold_content(args, self._logger)

    if self.density_flag:
        image = self.outputs['density'][self.ixe]

        # self.density is assigned in process()
        args['density'] = self.density
        args['image'] = image
        args['title'] = 'Density {}'.format(args['report_date'])

        fig_names['density'] = density(args, self._logger)

    if self.stn_validate_flag:
        stn_validate(self)

    if self.compare_runs_flag:
        compare_runs(self)

    if self.precip_depth_flag:
        precip_depth(self)

    # # Write out current model SWE values at snow course locations
    # if self.write_stn_csv_flag is True:
    #     point_values(self.outputs['swe_z'][-1],
    #                  self.stns_csv,
    #                  (self.snow_x, self.snow_y),
    #                  '{}model_pixel_swe_{}.csv'.format(self.figs_path,
    #                  self.end_date.date().strftime("%Y%m%d")))

    if self.precip_validate_flag:
        precip_validate(self)

    if self.basin_detail_flag:
        basin_detail(self)

    # Make flight difference figure in options in config file
    if self.flt_flag:
        flt_image_change(self)

    if self.basin_total_flag:

        '''
        clean up this and basin_total() still...

        '''
        wy_start = datetime(self.wy-1,10,1)
        swi_summary = collect(connector, args['plotorder'], args['basins'],
                              wy_start,args['end_date'],'swi_vol',
                              args['run_name'],'total','daily')
        df_swe = collect(connector, args['plotorder'], args['basins'],
                              wy_start,args['end_date'],'swe_vol',
                              args['run_name'],'total','daily')
        df_swi = swi_summary.cumsum()

        args['swi_summary'] = df_swi
        args['swe_summary'] = df_swe
        args['flag'] = False
        args['wy'] = self.wy

        fig_names['basin_total'] = basin_total(args, self._logger)

    if self.forecast_flag:

        if self.image_change_flag:
            image_change(self, forecast=self.for_run_name)

        if self.swi_flag:

            if forecast is None:
                run_name = snow.run_name
                outputs = copy.deepcopy(snow.outputs)
                ixs = snow.ixs
                ixe = snow.ixe
                start_date = snow.start_date
                end_date = snow.end_date
                directory = snow.directory
                title = 'Accumulated SWI\n{} to {}'.format(
                                            snow.report_start.date().strftime("%Y-%-m-%-d"),
                                            snow.report_date.date().strftime("%Y-%-m-%-d"))
            swi(args)

        if self.swe_volume_flag:
            '''
                if day is not None:
                    figs_path = day.figs_path
                    name_append = 'day'
                    date_stamp = day.date.strftime("%Y-%-m-%-d %H:%M") + ' (UTC)'

            '''

            args['directory'] = self.directory + '_forecast'
            args['run_name'] = self.for_run_name
            args['start_date'] = self.for_start_date
            args['end_date'] = self.for_end_date
            args['title'] = 'Forecast SWE \n {}'.format(self.for_end_date.date().strftime("%Y-%-m-%-d"))

            swe = collect(self,args['plotorder'],args['start_date'],
                          args['end_date'],'swe_vol',args['run_name'],
                          args['edges'],'end')

            image = self.for_outputs['swe_z'][self.for_ixe]*self.depth_factor

            args['df'] = swe
            args['image'] = image
            args['title'] = 'SWE {}'.format(args['report_end'])

            name, ylims = swe_volume(args, self._logger)

            if self.basin_total_flag:
                args['flag'] = True

                basin_total(self, forecast=self.for_run_name)

            if self.precip_depth_flag:
                precip_depth(self, forecast=self.for_run_name)


def save_fig(fig, paths):
    '''
    Save figures.

    Args
    ----------
    fig : object
        matplotlib figure object
    paths : str
        list of paths to save figure

    '''

    if type(paths) != list:
        paths = [paths]

    for path in paths:
        fig.savefig(path)
