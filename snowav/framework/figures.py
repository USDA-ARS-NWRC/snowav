
import os
import numpy as np
import pandas as pd
import copy
from datetime import datetime, timedelta
from tablizer.tablizer import get_existing_records
from snowav.plotting.swi import swi
from snowav.plotting.basin_total import basin_total
from snowav.plotting.cold_content import cold_content
from snowav.plotting.compare_runs import compare_runs
from snowav.plotting.density import density
from snowav.plotting.flt_image_change import flt_image_change
from snowav.plotting.image_change import image_change
from snowav.plotting.precip_depth import precip_depth
from snowav.plotting.stn_validate import stn_validate
from snowav.plotting.write_properties import write_properties
from snowav.plotting.swe_volume import swe_volume
from snowav.plotting.inputs import inputs
from snowav.inflow.inflow import inflow
from snowav.plotting.diagnostics import diagnostics
from snowav.plotting.plotlims import plotlims
from snowav.plotting.point_values import point_values
from snowav.database.database import collect
from snowav.plotting.plotlims import plotlims as plotlims
import matplotlib as mpl

if os.environ.get('DISPLAY','') == '':
    mpl.use('Agg')

def figures(self):
    '''
    Set up and call snowav figures. See CoreConfig.ini and README.md for more on
    config options and use.

    Notes:
    - swe_volume() must be called before cold_content() if you want to use
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
            'basins':self.basins,
            'wy':self.wy,
            'flag':False,
            'flt_flag':self.flt_flag}

    if self.flt_flag:
        args['flight_dates'] = self.flight_diff_dates

    fig_names = {}
    connector = self.connector

    ##########################################################################
    #       For each figure, collect 2D array image, by-elevation            #
    #       DataFrame, and set any figure-specific args inputs               #
    ##########################################################################
    if self.flt_flag:
        args['depth_factor'] = self.depth_factor
        args['update_file'] = self.update_file
        args['update_numbers'] = self.update_numbers
        args['flight_outputs'] = self.flight_outputs
        args['pre_flight_outputs'] = self.pre_flight_outputs
        args['connector'] = connector

        self.flight_diff_fig_names, self.flight_delta_vol_df = flt_image_change(args, self._logger)

        if self.flight_diff_fig_names == []:
            self.flt_flag = False

    if self.swi_flag:
        image = np.zeros_like(self.outputs['swi_z'][0])
        for n in range(self.ixs,self.ixe+1):
            image = image + self.outputs['swi_z'][n]*self.depth_factor

        df = collect(connector, args['plotorder'], args['basins'],
                     args['start_date'], args['end_date'], 'swi_vol',
                     args['run_name'], args['edges'], 'sum')

        args['df'] = df
        args['image'] = image
        args['title'] = 'Accumulated SWI\n{} to {}'.format(
                        args['report_start'],
                        args['report_date'])

        fig_names['swi'] = swi(args, self._logger)

    if self.image_change_flag:
        image = self.outputs['swe_z'][self.ixe] - self.outputs['swe_z'][self.ixs]

        start = collect(connector, args['plotorder'], args['basins'],
                        args['start_date'], args['start_date'],'swe_vol',
                        args['run_name'], args['edges'],'end')
        end = collect(connector, args['plotorder'], args['basins'],
                      args['start_date'], args['end_date'], 'swe_vol',
                      args['run_name'], args['edges'],'end')

        df = end - start

        args['df'] = df
        args['image'] = image*self.depth_factor
        args['title'] = 'Change in SWE\n{} to {}'.format(
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

    if self.basin_total_flag:
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
        args['forecast_flag'] = self.forecast_flag
        args['flt_flag'] = self.flt_flag

        fig_names['basin_total'] = basin_total(args, self._logger)

    if self.precip_depth_flag:
        swi_image = np.zeros_like(self.outputs['swi_z'][0])
        for n in range(self.ixs,self.ixe+1):
            swi_image = swi_image + self.outputs['swi_z'][n]*self.depth_factor

        swi_df = collect(connector, args['plotorder'], args['basins'],
                         args['start_date'], args['end_date'], 'swi_z',
                         args['run_name'], args['edges'], 'sum')
        precip_df = collect(connector, args['plotorder'], args['basins'],
                            args['start_date'], args['end_date'], 'precip_z',
                            args['run_name'], args['edges'], 'sum')
        rain_df = collect(connector, args['plotorder'], args['basins'],
                          args['start_date'], args['end_date'], 'rain_z',
                          args['run_name'], args['edges'], 'sum')

        args['swi_image'] = swi_image
        args['precip_image'] = self.precip_total*self.depth_factor
        args['rain_image'] = self.rain_total*self.depth_factor
        args['swi_df'] = swi_df
        args['precip_df'] = precip_df
        args['rain_df'] = rain_df
        args['title'] = 'Depth of SWI, Precipitation, and Rain\n{} to {}'.format(
                        args['report_start'],args['report_date'])

        fig_names['precip_depth'] = precip_depth(args, self._logger)

    if self.diagnostics_flag:
        wy_start = datetime(self.wy-1,10,1)
        precip = collect(connector, args['plotorder'], args['basins'],
                         wy_start, args['end_date'], 'precip_z',
                         args['run_name'], args['edges'], 'daily')
        precip_per = collect(connector, args['plotorder'], args['basins'],
                         args['start_date'], args['end_date'], 'precip_z',
                         args['run_name'], args['edges'], 'daily')
        swe = collect(connector, args['plotorder'], args['basins'],
                         wy_start, args['end_date'], 'swe_z',
                         args['run_name'], args['edges'], 'daily')
        swe_per = collect(connector, args['plotorder'], args['basins'],
                         args['start_date'], args['end_date'], 'swe_z',
                         args['run_name'], args['edges'], 'daily')
        rho = collect(connector, args['plotorder'], args['basins'],
                         wy_start, args['end_date'], 'density',
                         args['run_name'], args['edges'], 'daily')
        rho_per = collect(connector, args['plotorder'], args['basins'],
                         args['start_date'], args['end_date'], 'density',
                         args['run_name'], args['edges'], 'daily')
        snow_line = collect(connector, args['plotorder'], args['basins'],
                         wy_start, args['end_date'], 'snow_line',
                         args['run_name'], args['edges'], 'daily')
        snow_line_per = collect(connector, args['plotorder'], args['basins'],
                         args['start_date'], args['end_date'], 'snow_line',
                         args['run_name'], args['edges'], 'daily')
        evap_z = collect(connector, args['plotorder'], args['basins'],
                         wy_start, args['end_date'], 'evap_z',
                         args['run_name'], args['edges'], 'daily')
        evap_z_per = collect(connector, args['plotorder'], args['basins'],
                         args['start_date'], args['end_date'], 'evap_z',
                         args['run_name'], args['edges'], 'daily')

        snow_line_per = snow_line_per.fillna(0)
        first_row = snow_line_per.iloc[[0]].values[0]
        snow_line_per = snow_line_per.apply(lambda row: row - first_row, axis=1)
        args['snow_line'] = snow_line
        args['snow_line_per'] = snow_line_per

        swe = swe.fillna(0)
        swe_per = swe_per.fillna(0)
        first_row = swe_per.iloc[[0]].values[0]
        swe_per = swe_per.apply(lambda row: row - first_row, axis=1)

        evap_z = evap_z.fillna(0)
        evap_z_per = evap_z_per.fillna(0)
        first_row = evap_z_per.iloc[[0]].values[0]
        evap_z_per = evap_z_per.apply(lambda row: row - first_row, axis=1)

        rho = rho.fillna(0)
        rho_per = rho_per.fillna(0)
        first_row = rho_per.iloc[[0]].values[0]
        rho_per = rho_per.apply(lambda row: row - first_row, axis=1)

        precip = precip.fillna(0)
        precip_per = precip_per.fillna(0)
        precip = precip.cumsum()
        precip_per = precip_per.cumsum()
        first_row = precip_per.iloc[[0]].values[0]
        precip_per = precip_per.apply(lambda row: row - first_row, axis=1)

        if self.diag_basins is None:
            args['dbasins'] = copy.deepcopy(self.plotorder)
        else:
            args['dbasins'] = self.diag_basins

        args['precip'] = precip
        args['precip_per'] = precip_per
        args['swe'] = swe
        args['swe_per'] = swe_per
        args['evap_z'] = evap_z
        args['evap_z_per'] = evap_z_per
        args['density'] = rho
        args['density_per'] = rho_per
        args['elevlbl'] = self.elevlbl

        diagnostics(args, self._logger)

    if self.stn_validate_flag:
        args['dirs'] = self.all_dirs
        args['stns'] = self.val_stns
        args['lbls'] = self.val_lbls
        args['client'] = self.val_client
        args['factor'] = 25.4
        args['user'] = self.wxdb_user
        args['password'] = self.wxdb_password
        args['host'] = self.wxdb_host
        args['port'] = self.wxdb_port
        args['snow_x'] = self.snow_x
        args['snow_y'] = self.snow_y
        args['stns'] = self.val_stns
        args['nash_sut_flag'] = self.nash_sut_flag

        # could change these for precip, px = 0, 'em.nc'
        args['tbl'] = 'tbl_level1'
        args['var'] = 'snow_water_equiv'
        args['px'] = (1,1,1,0,0,0,-1,-1,-1)
        args['py'] = (1,0,-1,1,0,-1,1,0,-1)
        args['ncfile'] = 'snow.nc'

        fig_names['valid'], flag = stn_validate(args, self._logger)

        if not flag:
            self.stn_validate_flag = False

    if self.point_values_flag and self.point_values:

        # check that self.point_values_date falls within options
        pv_date = self.point_values_date

        if pv_date is None:
            self._logger.info(' Value in [validate] point_values_date being '
                              'assigned to {}'.format(
                              self.end_date.date().strftime("%Y-%m-%d")))
            pv_date = self.end_date

        if pv_date < self.start_date or pv_date > self.end_date:
            self._logger.info(' Value in [validate] point_values_date outside '
                              'of range in [run] start_date - end_date, '
                              'point_values_date being assigned to {}'.format(
                              self.end_date.date().strftime("%Y-%m-%d")))
            pv_date = self.end_date
            idx = -1

        else:
            # get index for self.outputs for that date
            x = np.abs([date - pv_date for date in self.outputs['dates']])
            idx = x.argmin(0)

        model_date = self.outputs['dates'][idx]
        model_date = model_date.date().strftime('%Y-%m-%d')
        course_date = self.point_values_date.date().strftime('%Y-%m-%d')
        nsubplots = (self.point_values_settings[3]*self.point_values_settings[4]-1)-1

        for idxp, value in enumerate(self.point_values_properties):
            pflag = True
            df = pd.read_csv(self.point_values_csv[idxp])

            if self.point_values_heading[idxp] is not None:
                check_headings = [self.point_values_heading[idxp],'name',
                                  'latitude','longitude']
            else:
                check_headings = ['name', 'latitude', 'longitude']

            for head in check_headings:
                if head not in df.columns.tolist():
                    self._logger.warn(' Config option [validate] '
                        'point_values_heading: {} not in headings {} in '
                        '{}, setting point_values: False'.format(head,
                        df.columns.tolist(), self.point_values_csv[idxp]))
                    pflag = False

            if not pflag:
                continue

            if len(df.name.unique()) > nsubplots:
                self._logger.warn(' Number of subplots that will be generated in '
                                  'point_values() may not fit well with settings '
                                  'in point_values_settings, consider changing '
                                  'nrows and/or ncols in [validate] '
                                  'point_values_settings')

            fig_name = '{}model_pixel_{}_{}.csv'.format(self.figs_path,
                        value, self.end_date.date().strftime("%Y%m%d"))

            if value == 'swe_z':
                factor = self.depth_factor

            if value == 'density':
                factor = 1

            if value == 'depth':
                factor = 39.37

            array = self.outputs[value][idx]*factor

            point_values(array, value, df, (self.snow_x, self.snow_y), fig_name,
                         self.dem, self.figs_path, self.veg_type,
                         self.point_values_heading[idxp], model_date, course_date,
                         self.point_values_settings, self.pixel, self._logger)

    if self.compare_runs_flag:
        args['variables'] = ['swe_vol','swi_vol']

        if self.flt_flag:
            args['flag'] = True
        else:
            args['flag'] = False

        dict = {}
        for var in args['variables']:
            dict[var] = {}
            for wy, run in zip(self.compare_run_wys, self.compare_run_names):
                wy_start = datetime(wy-1,10,1)
                df = collect(connector, args['plotorder'][0], args['basins'],
                             wy_start,args['end_date'],var,run,'total','daily')

                if wy != self.wy:
                    adj = self.wy - wy
                    df.index = df.index + timedelta(days = 365*adj)

                if var == 'swi_vol':
                    df = df.cumsum()

                dict[var][run] = df

        args['dict'] = dict

        compare_runs(args, self._logger)

    if self.inflow_flag:
        wy_start = datetime(self.wy-1,10,1)
        swi_summary = collect(connector, args['plotorder'], args['basins'],
                              wy_start,args['end_date'],'swi_vol',
                              args['run_name'],'total','daily')
        df_swi = swi_summary.cumsum()

        args['swi_summary'] = df_swi

        if self.inflow_data is None:
            raw = pd.read_csv(self.summary_csv, skiprows = 1,
                              parse_dates=[0], index_col = 0)
            args['inflow_summary'] = raw[self.basin_headings]

        else:
            args['inflow_summary'] = pd.read_csv(self.summary_csv,
                                                 parse_dates=[0], index_col = 0)


        args['inflow_headings'] = self.inflow_headings
        args['basin_headings'] = self.basin_headings

        inflow(args, self._logger)

    if self.write_properties is not None:
        args['connector'] = self.connector
        args['wy_start'] = datetime(self.wy-1,10,1)

        write_properties(args, self.write_properties, self._logger)

    if self.inputs_fig_flag:

        if self.mysql is not None:
            dbs = 'sql'
        else:
            dbs = 'sqlite'

        df = get_existing_records(connector, dbs)
        df = df.set_index('date_time')
        df.sort_index(inplace=True)

        ivalue = {}
        p = []

        for var in self.plots_inputs_variables:
            ivalue[var] = {}

            for basin in self.inputs_basins:
                bid = args['basins'][basin]['basin_id']
                ivalue[var][basin] = {}

                for func in self.inputs_methods:

                    if 'percentile' in func:
                        nfunc = '{}_{}'.format(func,str(self.inputs_percentiles[0]))

                        if ((var == self.plots_inputs_variables[0]) and
                             (basin == self.inputs_basins[0])):
                            p.append(nfunc)

                        ivalue[var][basin][nfunc] =  df[(df['function'] == nfunc) &
                                       (df['variable'] == var) &
                                       (df['basin_id'] == int(bid)) &
                                       (df['run_name'] == args['run_name'])]

                        nfunc = '{}_{}'.format(func,str(self.inputs_percentiles[1]))

                        if ((var == self.plots_inputs_variables[0]) and
                             (basin == self.inputs_basins[0])):
                            p.append(nfunc)

                        ivalue[var][basin][nfunc] =  df[(df['function'] == nfunc) &
                                       (df['variable'] == var) &
                                       (df['basin_id'] == int(bid)) &
                                       (df['run_name'] == args['run_name'])]
                    else:
                        ivalue[var][basin][func] =  df[(df['function'] == func) &
                                       (df['variable'] == var) &
                                       (df['basin_id'] == int(bid)) &
                                       (df['run_name'] == args['run_name'])]

                        if ((var == self.plots_inputs_variables[0]) and
                             (basin == self.inputs_basins[0])):
                            p.append(func)

        args['inputs'] = ivalue
        args['inputs_methods'] = p
        args['var_list'] = self.plots_inputs_variables
        args['inputs_basins'] = self.inputs_basins

        inputs(args, self._logger)

    if self.forecast_flag:
        print('Forecast figures in progress...')
        # if self.image_change_flag:
        #     image_change(self, forecast=self.for_run_name)
        #
        # if self.swi_flag:
        #
        #     if forecast is None:
        #         run_name = snow.run_name
        #         outputs = copy.deepcopy(snow.outputs)
        #         ixs = snow.ixs
        #         ixe = snow.ixe
        #         start_date = snow.start_date
        #         end_date = snow.end_date
        #         directory = snow.directory
        #         title = 'Accumulated SWI\n{} to {}'.format(
        #                                     snow.report_start.date().strftime("%Y-%-m-%-d"),
        #                                     snow.report_date.date().strftime("%Y-%-m-%-d"))
        #     swi(args)
        #
        # if self.swe_volume_flag:
        #     '''
        #         if day is not None:
        #             figs_path = day.figs_path
        #             name_append = 'day'
        #             date_stamp = day.date.strftime("%Y-%-m-%-d %H:%M") + ' (UTC)'
        #
        #     '''
        #
        #     args['directory'] = self.directory + '_forecast'
        #     args['run_name'] = self.for_run_name
        #     args['start_date'] = self.for_start_date
        #     args['end_date'] = self.for_end_date
        #     args['title'] = 'Forecast SWE \n {}'.format(self.for_end_date.date().strftime("%Y-%-m-%-d"))
        #
        #     swe = collect(self,args['plotorder'],args['start_date'],
        #                   args['end_date'],'swe_vol',args['run_name'],
        #                   args['edges'],'end')
        #
        #     image = self.for_outputs['swe_z'][self.for_ixe]*self.depth_factor
        #
        #     args['df'] = swe
        #     args['image'] = image
        #     args['title'] = 'SWE {}'.format(args['report_end'])
        #
        #     name, ylims = swe_volume(args, self._logger)
        #
        #     if self.basin_total_flag:
        #         '''
        #             for iter,d in enumerate(v['date_time'].values):
        #                 swe_summary.loc[d,bid] = v['value'].values[iter]
        #                 swi_summary.loc[d,bid] = v2['value'].values[iter]
        #
        #         swi_summary.sort_index(inplace=True)
        #
        #         # as a starting spot, add actual run
        #         swi_summary.iloc[0,:] = swi_summary.iloc[0,:] + swi_end_val.iloc[-1,:].values
        #         swi_summary = swi_summary.cumsum()
        #         '''
        #         # forecast True
        #         args['flag'] = True
        #
        #         basin_total(self, forecast=self.for_run_name)
        #
        #     if self.precip_depth_flag:
        #         precip_depth(self, forecast=self.for_run_name)

    self.fig_names = fig_names


def save_fig(fig, paths):
    '''
    Args
    ----------
    fig : object
        matplotlib figure object
    paths : list
        list of paths to save figure

    '''

    if type(paths) != list:
        paths = [paths]

    for path in paths:
        fig.savefig(path)
