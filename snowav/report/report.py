
from jinja2 import FileSystemLoader
from latex.jinja2 import make_env
from latex import build_pdf
import pandas as pd
import numpy as np
from datetime import datetime
import copy
from snowav import database
from snowav.database.tables import Basins, Watersheds
import pandas as pd
import os
from datetime import timedelta
import snowav
from snowav.database.database import collect


def report(self):
    '''
    This makes the pdf report.

    The latex formats for figures, tables, and summary paragraphs are created
    by templates in snowav/report/. Data is pulled from the snowav database.

    To add a figure to the report:
        - create the plot (of course)
        - add to variables dict in this file, with the same naming convention
            variables['PDEP_FIG'] = 'precip_depth%s.png'%(self.name_append)
        - add figure template to snowav_report.tex
            \VAR{PDEP_FIG_TPL}
        - add to section_dict
        - add template to snowav/report/figs/
        - add to snowav/config/CoreConfig.py [report] exclude_figs
        - if the figure may not exist (such as flt_diff, or those that require
          process() to run), address that with some form of exception before
          creating the pdf

    '''

    bid = Basins.basins[self.plotorder[0]]['basin_id']
    wy_start = datetime(self.wy-1,10,1)
    start_date = self.start_date
    end_date = self.end_date
    run_name = self.run_name
    edges = self.edges
    plotorder = self.plotorder
    dpts = int(self.dplcs)

    # Initialize variables to pass to latex file
    variables = {}

    dbval = collect(self,plotorder[0],wy_start,end_date,'swi_vol',run_name,'total','sum')
    variables['TOTAL_SWI'] = dbval.sum().values.round(decimals=dpts)[0]

    dbval = collect(self,plotorder[0],start_date,end_date,'swi_vol',run_name,'total','sum')
    variables['PER_SWI'] = dbval.sum().values.round(decimals=dpts)[0]

    dbval = collect(self,plotorder[0],start_date,end_date,'swe_vol',run_name,'total','end')
    variables['TOTAL_SWE'] = dbval.sum().values.round(decimals=dpts)[0]

    dbval = collect(self,plotorder[0],start_date,end_date,'swe_avail',run_name,'total','end')
    variables['TOTAL_SAV'] = dbval.sum().values.round(decimals=dpts)[0]

    dbval = collect(self,plotorder[0],start_date,end_date,'swe_z',run_name,'total','end')
    variables['TOTAL_PM'] = dbval.sum().values.round(decimals=dpts)[0]

    dbval = collect(self,plotorder[0],wy_start,end_date,'precip_z',run_name,'total','sum')
    variables['TOTALPRE_PM'] = dbval.sum().values.round(decimals=dpts)[0]

    s = collect(self,plotorder[0],start_date,start_date,'swe_vol',run_name,'total','end')
    e = collect(self,plotorder[0],start_date,end_date,'swe_vol',run_name,'total','end')
    start_swe = s.sum().values.round(decimals=dpts)[0]
    end_swe = e.sum().values.round(decimals=dpts)[0]
    diff = end_swe - start_swe
    variables['TOTAL_SDEL'] = diff
    variables['TOTAL_PVOL'] = '100'

    if float(end_swe) - float(start_swe) > 0:
        variables['SIGN'] = r'$+$'
    if float(end_swe) - float(start_swe) == 0.0:
        variables['SIGN'] = ''
    if float(end_swe) - float(start_swe) < 0.0:
        variables['SIGN'] = r'-'

    dbval = collect(self,plotorder[0],wy_start,end_date,'rain_z',run_name,'total','sum')
    total_rain = dbval.sum().values.round(decimals=dpts)[0]

    if variables['TOTALPRE_PM'] != 0:
        variables['TOTAL_RAT'] = str(int((total_rain/variables['TOTALPRE_PM'])*100))
    else:
        variables['TOTAL_RAT'] = '0'

    report_time = datetime.now().strftime("%Y-%-m-%-d %H:%M")
    # if hasattr(self,'orig_date'):
    #     report_time = self.orig_date + ', revised ' + report_time

    numsubs = range(1,len(self.plotorder))

    for n,sub in zip(numsubs,self.plotorder[1:]):
        SWIIND = 'SUB' + str(n) + '_SWI'
        PERSWIIND = 'SUB' + str(n) + '_PERSWI'
        SWEIND = 'SUB' + str(n) + '_SWE'
        AVSWEIND = 'SUB' + str(n) + '_SAV'
        SWEDEL = 'SUB' + str(n) + '_SDEL'
        PM = 'SUB' + str(n) + '_PM'
        PREPM = 'SUB' + str(n) + 'PRE_PM'
        RAIN = 'SUB' + str(n) + 'RAI'
        RATIO = 'SUB' + str(n) + '_RAT'
        PVOL = 'SUB' + str(n) + '_PVOL'

        dbval = collect(self,sub,wy_start,end_date,'swi_vol',run_name,'total','sum')
        variables[SWIIND] = dbval.sum().values.round(decimals=dpts)[0]

        dbval = collect(self,sub,start_date,end_date,'swi_vol',run_name,'total','sum')
        variables[PERSWIIND] = dbval.sum().values.round(decimals=dpts)[0]

        dbval = collect(self,sub,start_date,end_date,'swe_vol',run_name,'total','end')
        variables[SWEIND] = dbval.sum().values.round(decimals=dpts)[0]

        dbval = collect(self,sub,start_date,end_date,'swe_avail',run_name,'total','end')
        variables[AVSWEIND] = dbval.sum().values.round(decimals=dpts)[0]

        dbval = collect(self,sub,start_date,start_date,'swe_vol',run_name,'total','end')
        start_swe = dbval.sum().values.round(decimals=dpts)[0]
        dbval = collect(self,sub,start_date,end_date,'swe_vol',run_name,'total','end')
        end_swe = dbval.sum().values.round(decimals=dpts)[0]
        variables[SWEDEL] = end_swe - start_swe

        dbval = collect(self,sub,start_date,end_date,'swe_z',run_name,'total','end')
        variables[PM] = dbval.sum().values.round(decimals=dpts)[0]

        dbval = collect(self,sub,wy_start,end_date,'precip_z',run_name,'total','sum')
        variables[PREPM] = dbval.sum().values.round(decimals=dpts)[0]

        dbval = collect(self,sub,wy_start,end_date,'rain_z',run_name,'total','sum')
        variables[RAIN] = dbval.sum().values.round(decimals=dpts)[0]

        if (end_swe > 0) and (variables['TOTAL_SWE'] > 0):
            variables[PVOL] = end_swe/variables['TOTAL_SWE']*100
        else:
            variables[PVOL] = '-'

        if variables[PREPM] != 0.0:
            variables[RATIO] = str(int((variables[RAIN]/variables[PREPM])*100))
        else:
            variables[RATIO] = '0'

    # Upper case variables are used in the LaTex file,
    # lower case versions are assigned here

    # untested - if report title contains comma?
    if isinstance(self.rep_title, list):
        title = self.rep_title[0]
        for s in self.rep_title[1::]:
            title = title + ', ' + s
        variables['REPORT_TITLE'] = title

    else:
        if self.flt_flag:
            variables['REPORT_TITLE'] = self.rep_title + r' \\ ASO Updates'
            fst = ' Model snow depths were updated with ASO snow depths on '

            for i,d in enumerate(self.flight_diff_datestr):
                if len(self.flight_diff_datestr) == 1:
                    fst = fst + d + '.'

                if ((len(self.flight_diff_datestr) > 1) and
                    (i < len(self.flight_diff_datestr) - 1)):
                    fst = fst + d + ', '

                if ((len(self.flight_diff_datestr) > 1) and
                    (i == len(self.flight_diff_datestr) - 1)):
                    fst = fst + 'and ' + d + '.'

            variables['FLTSENT'] = fst

        else:
            variables['REPORT_TITLE'] = self.rep_title
            variables['FLTSENT'] = (' Model snow depths will be updated with '
                                    'ASO snow depths as they become available.')

    variables['REPORT_TIME'] = report_time
    variables['WATERYEAR'] = str(self.wy)
    variables['UNITS'] = self.vollbl
    variables['VOLLBL'] = self.vollbl
    variables['DEPLBL'] = self.depthlbl
    variables['START_DATE'] = self.report_start.date().strftime("%B %-d")
    variables['END_DATE'] = self.report_date.date().strftime("%B %-d")
    variables['SWE_IN'] = variables['TOTAL_PM']
    variables['SWI_IN'] = variables['TOTAL_SWI']
    variables['FIG_PATH'] = self.figs_path
    variables['SWI_FIG'] = 'swi_{}.png'.format(self.name_append)
    variables['CHANGES_FIG'] = 'swe_change_{}.png'.format(self.name_append)
    variables['CHANGES_DEP_FIG'] = 'swe_change_depth_{}.png'.format(self.name_append)
    variables['TOTALS_FIG'] = 'basin_total_{}.png'.format(self.name_append)
    variables['MULTITOTSWE_FIG'] = 'compare_swe_vol_{}.png'.format(self.name_append)
    variables['HYP_FIG'] = 'hypsometry_{}.png'.format(self.name_append)
    variables['MEAN_FIG'] = 'mean_swe_depth_{}.png'.format(self.name_append)
    variables['DETAIL_FIG'] = 'mean_detail_{}.png'.format(self.name_append)
    variables['DENSITY_FIG'] = 'density_{}.png'.format(self.name_append)
    variables['DENSITY_SUB_FIG'] = 'density_sub_{}.png'.format(self.name_append)
    variables['DENSITY_SWE_FIG'] = 'density_swe_{}.png'.format(self.name_append)
    variables['PDEP_FIG'] = 'precip_depth_{}.png'.format(self.name_append)
    variables['VALID_FIG'] = 'validation_{}.png'.format(self.name_append)
    variables['COLD_FIG'] = 'cold_content_{}.png'.format(self.name_append)
    variables['SWE_FIG'] = 'swe_volume_{}.png'.format(self.name_append)
    variables['INFLOW_FIG'] = 'inflow_{}.png'.format(self.name_append)
    variables['VERSION'] = snowav.__version__

    if (self.update_file is not None) and self.flt_flag :
        for name in self.flight_diff_fig_names:
            variables['DFLT_FIG'] = name

    if self.forecast_flag is True:
        variables['FORE_START_DATE'] = self.for_start_date.date().strftime("%B %-d")
        variables['FORE_DATE'] = self.for_end_date.date().strftime("%B %-d")
        variables['SWEFORECAST_FIG'] = 'swe_volume_{}.png'.format(self.name_append + '_forecast')
        variables['SWIFORECAST_FIG'] = 'swi_{}.png'.format(self.name_append + '_forecast')
        variables['CHANGESFORECAST_FIG'] = 'swe_change_{}.png'.format(self.name_append + '_forecast')
        variables['TOTALSFORECAST_FIG'] = 'basin_total_{}.png'.format(self.name_append + '_forecast')
        variables['PDEPFORECAST_FIG'] = 'precip_depth_{}.png'.format(self.name_append + '_forecast')

    if self.subs_fig is not None:
        variables['SUBBASINS_FIG'] = '{}'.format(self.subs_fig)

    else:
        if self.exclude_figs != None:
            self.exclude_figs = self.exclude_figs + ['SUBBASINS']
        else:
            self.exclude_figs = ['SUBBASINS']

    # Logos
    variables['ARSLOGO'] = self.figs_tpl_path + 'ARS.jpg'
    variables['ASOLOGO'] = self.figs_tpl_path + 'ASO.jpg'
    variables['USDALOGO'] = self.figs_tpl_path + 'USDA.png'
    variables['JPLLOGO'] = self.figs_tpl_path + 'JPL.jpg'
    variables['CDWRLOGO'] = self.figs_tpl_path + 'CDWR.png'
    variables['USBRLOGO'] = self.figs_tpl_path + 'USBR.jpg'
    variables['NRCSLOGO'] = self.figs_tpl_path + 'NRCS.jpg'
    variables['KRWALOGO'] = self.figs_tpl_path + 'KRWA.jpg'
    variables['FRIANTLOGO'] = self.figs_tpl_path + 'FRIANT.jpg'
    variables['AWSMLOGO'] = self.figs_tpl_path + 'logo.png'

    dfind = [str(i) for i in self.edges]
    dfindt = [str(i) for i in self.edges] + ['total']
    colstr = 'l' + 'r'*len(self.plotorder)

    if self.basin in ['KINGS', 'SJ', 'KAWEAH','MERCED']:
        spacecmd = r'\resizebox{\textwidth}{!}{'
    else:
        spacecmd = r'{'

    dbval = collect(self,plotorder,start_date,end_date,'swe_z',run_name,dfindt,'end')
    swe_byelev = dbval.round(decimals=dpts)
    swe_byelev.rename(index={'total':'mean'},inplace=True)
    swe_byelev.index.name = 'Elevation'
    variables['SWE_BYELEV'] = (
                                r'  \textbf{SWE [%s], %s}\\ \vspace{0.1cm} \\'
                                %(self.depthlbl,self.report_date.date().strftime("%Y-%-m-%-d")) +
                                spacecmd + swe_byelev[plotorder].to_latex(na_rep='-', column_format=colstr) +
                                r'} \\ \footnotesize{\textbf{Table 2:} Mean depth of SWE by elevation band.}'
                                )

    dbval = collect(self,plotorder,start_date,end_date,'swi_vol',run_name,dfindt,'sum')
    swi_byelev = dbval.round(decimals=dpts)
    variables['ACCUM_BYELEV'] = (
                                r'  \textbf{SWI [%s] by elevation, %s to %s}\\ \vspace{0.1cm} \\'
                                %(self.vollbl,self.report_start.date().strftime("%Y-%-m-%-d"),
                                  self.report_date.date().strftime("%Y-%-m-%-d")) + spacecmd +
                                  swi_byelev[plotorder].to_latex(na_rep='-', column_format=colstr) +
                                  r'} \\ \footnotesize{\textbf{Table 6:} Volume of SWI during the report period.} ' +
                                  r'\\ \clearpage'
                                )

    dbval = collect(self,plotorder,start_date,end_date,'swe_vol',run_name,dfindt,'end')
    swe_byelev = dbval.round(decimals=dpts)
    swe_byelev.rename(index={'total':'mean'},inplace=True)
    swe_byelev.index.name = 'Elevation'
    variables['SWEVOL_BYELEV'] = (
                                r'  \textbf{SWE [%s], %s}\\ \vspace{0.1cm} \\'
                                %(self.vollbl,self.report_date.date().strftime("%Y-%-m-%-d")) +
                                spacecmd + swe_byelev[plotorder].to_latex(na_rep='-', column_format=colstr) +
                                r'} \\ \footnotesize{\textbf{Table 5:} Volume of SWE by elevation band.}'
                                )

    dbval = collect(self,plotorder,start_date,start_date,'swe_z',run_name,dfindt,'end')
    start_swe = dbval.round(decimals=dpts)
    dbval = collect(self,plotorder,start_date,end_date,'swe_z',run_name,dfindt,'end')
    end_swe = dbval.round(decimals=dpts)
    dswe_byelev = end_swe - start_swe
    dswe_byelev.rename(index={'total':'mean'},inplace=True)
    dswe_byelev.index.name = 'Elevation'
    variables['DSWE_BYELEV'] = (
                                r'  \textbf{Change in SWE [%s], %s to %s}\\ \vspace{0.1cm} \\'
                                %(self.depthlbl,self.report_start.date().strftime("%Y-%-m-%-d"),
                                  self.report_date.date().strftime("%Y-%-m-%-d")) + spacecmd +
                                  dswe_byelev[self.plotorder].to_latex(na_rep='-', column_format=colstr) +
                                  r'} \\ \footnotesize{\textbf{Table 3:} Change in depth of SWE by elevation band.} ' +
                                  r'\\ \clearpage'
                                )

    dbval = collect(self,plotorder,start_date,end_date,'swe_vol',run_name,dfind,'end')
    swe_byelev = dbval.round(decimals=dpts)
    value = swe_byelev.iloc[:-1].sum()
    sweper_byelev = (swe_byelev/value*100).round(decimals = dpts)
    sweper_byelev.index.name = 'Elevation'

    variables['SWEPER_BYELEV'] = (
                                r'  \textbf{SWE volume, percent of basin total by elevation band, %s}\\ \vspace{0.1cm} \\'
                                %(self.report_date.date().strftime("%Y-%-m-%-d"))+
                                spacecmd + sweper_byelev[plotorder].round(1).to_latex(na_rep='-', column_format=colstr) +
                                r'}  \\ \footnotesize{\textbf{Table 4:} Percent of SWE volume by elevation.} '
                                )
    variables['SWEPER_BYELEV'] = variables['SWEPER_BYELEV'].replace('inf','-')

    # If the basin total percents will be 0/nan, omit that table
    # if ((sub == self.plotorder[0]) and
    #     ((np.nansum(sum) < 0.001) or (np.nansum(sum) == np.nan))):
    #     sum_flag = False

    variables['TOT_LBL'] = self.plotorder[0]

    for n in range(1,len(self.plotorder)):
        s = 'SUB' + str(n) + '_LBL'
        variables[s] = self.plotorder[n]

    # Convert floats to strings
    for name in variables:
        if isinstance(variables[name], float):
            if self.dplcs == 0:
                tmp = str(int(variables[name]))
            else:
                tmp = str(round(variables[name],self.dplcs))
            variables[name] = tmp

    # Summary sections and fig template have variable strings
    # (e.g. CHANGES_FIG) that need to be replaced
    section_dict = {'SUMMARY':self.summary_file,
                    'CHANGES_FIG_TPL':self.figs_tpl_path + 'changes_fig_tpl.txt',
                    'SWI_FIG_TPL':self.figs_tpl_path + 'swi_fig_tpl.txt',
                    'MEAN_FIG_TPL':self.figs_tpl_path + 'mean_fig_tpl.txt',
                    'TOTALS_FIG_TPL':self.figs_tpl_path + 'totals_fig_tpl.txt',
                    'MULTITOTSWE_FIG_TPL':self.figs_tpl_path + 'multitotswe_fig_tpl.txt',
                    'VALID_FIG_TPL':self.figs_tpl_path + 'valid_fig_tpl.txt',
                    'FLTCHANGES_FIG_TPL':self.figs_tpl_path + 'flt_fig_tpl.txt',
                    'PDEP_FIG_TPL':self.figs_tpl_path + 'pdep_fig_tpl.txt',
                    'COLD_FIG_TPL':self.figs_tpl_path + 'cold_fig_tpl.txt',
                    'SWE_FIG_TPL':self.figs_tpl_path + 'swe_fig_tpl.txt',
                    'SUBBASINS_FIG_TPL':self.figs_tpl_path + 'subbasins_fig_tpl.txt',
                    'INFLOW_FIG_TPL':self.figs_tpl_path + 'inflow_fig_tpl.txt',
                    'SWEFORECAST_FIG_TPL':self.figs_tpl_path + 'forecastswe_fig_tpl.txt',
                    'SWIFORECAST_FIG_TPL':self.figs_tpl_path + 'forecastswi_fig_tpl.txt',
                    'CHANGESFORECAST_FIG_TPL':self.figs_tpl_path + 'forecastchanges_fig_tpl.txt',
                    'TOTALSFORECAST_FIG_TPL':self.figs_tpl_path + 'forecasttotals_fig_tpl.txt',
                    'PDEPFORECAST_FIG_TPL':self.figs_tpl_path + 'forecastpdep_fig_tpl.txt'
                    }

    # Define and load summary tables depending on number of subbasins
    section_dict['PRECIP_SUMMARY_TPL'] = (self.figs_tpl_path
                                          + 'precip_summary_%ssub.txt'%str(len(self.plotorder)) )
    section_dict['SWE_SUMMARY_TPL'] = (self.figs_tpl_path
                                          + 'swe_summary_%ssub.txt'%str(len(self.plotorder)) )

    # Remove if no flight options
    if self.flt_flag is False:
        del section_dict['FLTCHANGES_FIG_TPL']

    # Remove if no flight options
    if self.forecast_flag is False:
        del section_dict['SWEFORECAST_FIG_TPL']
        del section_dict['SWIFORECAST_FIG_TPL']
        del section_dict['CHANGESFORECAST_FIG_TPL']
        del section_dict['TOTALSFORECAST_FIG_TPL']
        del section_dict['PDEPFORECAST_FIG_TPL']

        variables['SWIFORECAST_FIG'] = ''
        variables['SWIFORECAST_FIG_TPL'] = ''
        variables['PDEPFORECAST_FIG'] = ''
        variables['PDEPFORECAST_FIG_TPL'] = ''
        variables['SWEFORECAST_FIG'] = ''
        variables['SWEFORECAST_FIG_TPL'] = ''
        variables['CHANGESFORECAST_FIG'] = ''
        variables['CHANGESFORECAST_FIG_TPL'] = ''
        variables['TOTALSFORECAST_FIG'] = ''
        variables['TOTALSFORECAST_FIG_TPL'] = ''
        variables['FORE_TITLE'] = ''

    else:
        variables['FORE_TITLE'] = 'with WRF forecast, {} to {}'.format(self.for_start_date.date().strftime("%b %-d"),
                    self.for_end_date.date().strftime("%b %-d"))

    # Remove if not available
    if self.plot_flag is True:
        del section_dict['PDEP_FIG_TPL']

    for rep in section_dict.keys():

        # for variable numbers of fligth figures
        if rep == 'FLTCHANGES_FIG_TPL':
            fid = open(section_dict[rep],'r')
            var = fid.read()
            fid.close()

            for name in sorted(variables):
                if name == 'DFLT_FIG':
                    for i,fltname in enumerate(self.flight_diff_fig_names):
                        flt_date = self.flight_outputs['dates'][i].date().strftime("%Y%m%d")
                        flt_num = self.flight_delta_vol_df[flt_date].round(1).to_latex(na_rep='-', column_format=colstr)
                        table = (r'{ \vspace{0.5cm} \textbf{Change in SWE [%s] by elevation, '%(self.vollbl) +
                                 r'from %s update} \\ \vspace{0.1cm} \\'%(flt_date) +
                                 r' %s %s }'%(spacecmd,flt_num))

                        if i == 0:
                            tmp = var.replace(name,fltname)
                        else:
                            tmp += var.replace(name,fltname)

                        tmp += table

                    var = tmp.replace(name,variables[name])

                else:
                    var = var.replace(name,variables[name])

            variables[rep] = var

        else:
            fid = open(section_dict[rep],'r')
            var = fid.read()
            fid.close()

            for name in sorted(variables):
                var = var.replace(name,variables[name])
            variables[rep] = var

    # Remove if not BRB
    if self.compare_runs_flag is False:
        variables['MULTITOTSWE_FIG'] = ''
        variables['MULTITOTSWE_FIG_TPL'] = ''

    # If figs are listed in exclude, replace with empty string in latex file
    if self.exclude_figs != None:
        for name in self.exclude_figs:
            variables[name + '_FIG'] = ''
            variables[name + '_TPL'] = ''
            variables[name + '_FIG_TPL'] = ''
            # variables[name] = ' '

    # Make the report
    env = make_env(loader = FileSystemLoader(self.templ_path))
    tpl = env.get_template(self.tex_file)
    # print(tpl.render(variables))
    pdf = build_pdf(tpl.render(variables))

    # Save in reports and with figs
    rpath_1 = os.path.join(self.rep_path, '' + self.report_name)
    rpath_2 = os.path.join(self.figs_path, '' + self.report_name)
    self._logger.info('Saving {}\nSaving {}'.format(rpath_1,rpath_2))

    if not os.path.isdir(os.path.join(self.rep_path,'')):
        os.makedirs(os.path.join(self.rep_path,''))

    pdf.save_to(rpath_1)
    pdf.save_to(rpath_2)
