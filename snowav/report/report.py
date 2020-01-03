
from jinja2 import FileSystemLoader
from latex.jinja2 import make_env
from latex import build_pdf
import pandas as pd
import numpy as np
from datetime import datetime
import copy
import os
import snowav
from snowav.database.database import collect

def report(self):
    '''
    Create the pdf report.

    See CoreConfig.ini and README.md for more information.

    The latex formats for figures, tables, and summary paragraphs are created
    by templates in snowav/report/. Data is pulled from the snowav database.

    To add a figure to the report:
    - create the plot (of course)
    - add to variables dict in this file, with the same naming convention
        variables['PDEP_FIG'] = 'precip_depth%s.png'%(self.directory)
    - add figure template to snowav_report.tex
        \VAR{PDEP_FIG_TPL}
    - add to section_dict
    - add template to snowav/report/figs/
    - add to snowav/config/CoreConfig.py [report] exclude_figs
    - if the figure may not exist (such as flt_diff, or those that require
      process() to run), address that with some form of exception before
      creating the pdf

    '''

    bid = self.basins[self.plotorder[0]]['basin_id']
    basins = self.basins
    wy_start = datetime(self.wy-1,10,1)
    start_date = self.start_date
    end_date = self.end_date
    run_name = self.run_name
    edges = self.edges
    plotorder = self.plotorder
    dpts = int(self.dplcs)
    ddpts = int(self.rep_dplcs)
    cnx = self.connector

    # variables to pass to latex file
    variables = {}

    dbval = collect(cnx,plotorder[0],basins,wy_start,end_date,'swi_vol',
                    run_name,'total','sum')
    variables['TOTAL_SWI'] = dbval.sum().values.round(decimals=dpts)[0]

    dbval = collect(cnx,plotorder[0],basins,start_date,end_date,'swi_vol',
                    run_name,'total','sum')
    variables['PER_SWI'] = dbval.sum().values.round(decimals=dpts)[0]

    dbval = collect(cnx,plotorder[0],basins,start_date,end_date,'swe_vol',
                    run_name,'total','end')
    variables['TOTAL_SWE'] = dbval.sum().values.round(decimals=dpts)[0]

    dbval = collect(cnx,plotorder[0],basins,start_date,end_date,'swe_avail',
                    run_name,'total','end')
    variables['TOTAL_SAV'] = dbval.sum().values.round(decimals=dpts)[0]

    dbval = collect(cnx,plotorder[0],basins,start_date,end_date,'swe_z',
                    run_name,'total','end')
    variables['TOTAL_PM'] = dbval.sum().values.round(decimals=ddpts)[0]

    if self.precip_flag:
        dbval = collect(cnx,plotorder[0],basins,wy_start,end_date,'precip_z',
                        run_name,'total','sum')
        variables['TOTALPRE_PM'] = dbval.sum().values.round(decimals=ddpts)[0]
    else:
        variables['TOTALPRE_PM'] = 0

    s = collect(cnx,plotorder[0],basins,start_date,start_date,'swe_vol',
                    run_name,'total','end')
    e = collect(cnx,plotorder[0],basins,start_date,end_date,'swe_vol',
                    run_name,'total','end')
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

    if self.precip_flag:
        dbval = collect(cnx,plotorder[0],basins,wy_start,end_date,'rain_z',
                        run_name,'total','sum')
        total_rain = dbval.sum().values.round(decimals=ddpts)[0]
    else:
        total_rain = 0.0

    if variables['TOTALPRE_PM'] != 0:
        variables['TOTAL_RAT'] = str(int((total_rain/variables['TOTALPRE_PM'])*100))
    else:
        variables['TOTAL_RAT'] = '-'
        variables['TOTALPRE_PM'] = '-'

    report_time = datetime.now().strftime("%Y-%-m-%-d %H:%M")

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

        dbval = collect(cnx,sub,basins,wy_start,end_date,'swi_vol',
                        run_name,'total','sum')
        variables[SWIIND] = dbval.sum().values.round(decimals=dpts)[0]

        dbval = collect(cnx,sub,basins,start_date,end_date,'swi_vol',
                        run_name,'total','sum')
        variables[PERSWIIND] = dbval.sum().values.round(decimals=dpts)[0]

        dbval = collect(cnx,sub,basins,start_date,end_date,'swe_vol',
                        run_name,'total','end')
        variables[SWEIND] = dbval.sum().values.round(decimals=dpts)[0]

        dbval = collect(cnx,sub,basins,start_date,end_date,'swe_avail',
                        run_name,'total','end')
        variables[AVSWEIND] = dbval.sum().values.round(decimals=dpts)[0]

        dbval = collect(cnx,sub,basins,start_date,start_date,'swe_vol',
                        run_name,'total','end')
        start_swe = dbval.sum().values.round(decimals=dpts)[0]
        dbval = collect(cnx,sub,basins,start_date,end_date,'swe_vol',
                        run_name,'total','end')
        end_swe = dbval.sum().values.round(decimals=dpts)[0]
        variables[SWEDEL] = end_swe - start_swe

        dbval = collect(cnx,sub,basins,start_date,end_date,'swe_z',
                        run_name,'total','end')
        variables[PM] = dbval.sum().values.round(decimals=ddpts)[0]

        dbval = collect(cnx,sub,basins,wy_start,end_date,'precip_z',
                        run_name,'total','sum')
        variables[PREPM] = dbval.sum().values.round(decimals=ddpts)[0]

        dbval = collect(cnx,sub,basins,wy_start,end_date,'rain_z',
                        run_name,'total','sum')
        variables[RAIN] = dbval.sum().values.round(decimals=ddpts)[0]

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
        if self.flt_flag and self.flight_figs:
            variables['REPORT_TITLE'] = self.rep_title + r' \\ ASO Updates'
            fst = ' Model snow depths were updated with ASO snow depths on '
            for i,d in enumerate(self.flight_diff_dates):
                if d <= self.end_date:
                    dn = d.date().strftime("%m/%d")
                    if len(self.flight_diff_dates) == 1:
                        fst = fst + dn + '.'

                    if ((len(self.flight_diff_dates) > 1) and
                        (i < len(self.flight_diff_dates) - 1)):
                        fst = fst + dn + ', '

                    if ((len(self.flight_diff_dates) > 1) and
                        (i == len(self.flight_diff_dates) - 1)):
                        fst = fst + 'and ' + dn + '.'

                else:
                    fst = fst.split(dn)[0] + 'and ' + dn + '.'
                    break

            variables['FLTSENT'] = fst

        else:
            variables['REPORT_TITLE'] = self.rep_title
            # variables['FLTSENT'] = (' Model snow depths will be updated with '
            #                         'ASO snow depths as they become available.')
            variables['FLTSENT'] = ('')

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
    variables['SWI_FIG'] = 'swi_{}.png'.format(self.directory)
    variables['CHANGES_FIG'] = 'swe_change_{}.png'.format(self.directory)
    variables['TOTALS_FIG'] = 'basin_total_{}.png'.format(self.directory)
    variables['MULTITOTSWE_FIG'] = 'compare_swe_vol_{}.png'.format(self.directory)
    variables['DENSITY_FIG'] = 'density_{}.png'.format(self.directory)
    variables['PDEP_FIG'] = 'precip_depth_{}.png'.format(self.directory)
    variables['VALID_FIG'] = 'validation_{}.png'.format(self.directory)
    variables['COLD_FIG'] = 'cold_content_{}.png'.format(self.directory)
    variables['SWE_FIG'] = 'swe_volume_{}.png'.format(self.directory)
    variables['VERSION'] = snowav.__version__

    if (self.update_file is not None) and self.flt_flag and self.flight_figs:
        for name in self.flight_diff_fig_names:
            variables['DFLT_FIG'] = name

    if self.forecast_flag:
        variables['FORE_START_DATE'] = self.for_start_date.date().strftime("%B %-d")
        variables['FORE_DATE'] = self.for_end_date.date().strftime("%B %-d")
        variables['SWEFORECAST_FIG'] = 'swe_volume_{}.png'.format(self.directory + '_forecast')
        variables['SWIFORECAST_FIG'] = 'swi_{}.png'.format(self.directory + '_forecast')
        variables['CHANGESFORECAST_FIG'] = 'swe_change_{}.png'.format(self.directory + '_forecast')
        variables['TOTALSFORECAST_FIG'] = 'basin_total_{}.png'.format(self.directory + '_forecast')
        variables['PDEPFORECAST_FIG'] = 'precip_depth_{}.png'.format(self.directory + '_forecast')

    if self.report_diagnostics:
        variables['DIAGNOSTICS_FIG'] = 'diagnostics_{}'.format(self.directory)
        variables['INPUTS_FIG'] = 'inputs_period_{}'.format(self.directory)

    if self.subs_fig is not None:
        variables['SUBBASINS_FIG'] = '{}'.format(self.subs_fig)

    # Logos
    variables['ARSLOGO'] = os.path.join(self.figs_tpl_path, 'ARS.jpg')
    variables['ASOLOGO'] = os.path.join(self.figs_tpl_path, 'ASO.jpg')
    variables['USDALOGO'] = os.path.join(self.figs_tpl_path, 'USDA.png')
    variables['JPLLOGO'] = os.path.join(self.figs_tpl_path, 'JPL.jpg')
    variables['CDWRLOGO'] = os.path.join(self.figs_tpl_path, 'CDWR.png')
    variables['USBRLOGO'] = os.path.join(self.figs_tpl_path, 'USBR.jpg')
    variables['NRCSLOGO'] = os.path.join(self.figs_tpl_path, 'NRCS.jpg')
    variables['KRWALOGO'] = os.path.join(self.figs_tpl_path, 'KRWA.jpg')
    variables['FRIANTLOGO'] = os.path.join(self.figs_tpl_path, 'FRIANT.jpg')
    variables['AWSMLOGO'] = os.path.join(self.figs_tpl_path, 'logo.png')

    dfind = [str(i) for i in self.edges]
    dfindt = [str(i) for i in self.edges] + ['total']
    colstr = 'l' + 'r'*len(self.plotorder)

    if len(self.plotorder) > 5:
        spacecmd = r'\resizebox{\textwidth}{!}{'
    else:
        spacecmd = r'{'

    ntables = len(self.tables)
    mtables = 2
    ptables = 0

    if 'swe_depth' in self.tables:

        dbval = collect(cnx,plotorder,basins,start_date,end_date,'swe_z',run_name,dfindt,'end')
        dbval = dbval.rename(columns=self.labels)
        swe_byelev = dbval.round(decimals=dpts)
        swe_byelev.rename(index={'total':'mean'},inplace=True)
        swe_byelev.index.name = 'Elevation'
        variables['SWE_BYELEV'] = (r' \normalsize \textbf{SWE [%s], %s}\\ \vspace{0.1cm} \\'
                                    %(self.depthlbl,self.report_date.date().strftime("%Y-%-m-%-d")) +
                                    spacecmd + swe_byelev.to_latex(na_rep='-', column_format=colstr) +
                                    r'} \\ \footnotesize{\textbf{Table %s:} SWE depth.}'%(str(mtables)))
        mtables += 1
        ptables += 1

    else:
        variables['SWE_BYELEV'] = ''

    if 'swe_vol' in self.tables:
        ptables += 1

        if ptables == 2:
            clrpage = r'\clearpage'
        else:
            clrpage = ''
        dbval = collect(cnx,plotorder,basins,start_date,end_date,'swe_vol',run_name,dfindt,'end')
        dbval = dbval.rename(columns=self.labels)
        swe_byelev = dbval.round(decimals=dpts)
        swe_byelev.index.name = 'Elevation'
        variables['SWEVOL_BYELEV'] = (r' \normalsize \textbf{SWE [%s], %s}\\ \vspace{0.1cm} \\'
                                    %(self.vollbl,self.report_date.date().strftime("%Y-%-m-%-d")) +
                                    spacecmd + swe_byelev.to_latex(na_rep='-', column_format=colstr) +
                                    r'} \\ \footnotesize{\textbf{Table %s:} SWE volume.}%s'%(str(mtables),clrpage))
        mtables += 1

    else:
        variables['SWEVOL_BYELEV'] = ''

    if 'swe_change' in self.tables:
        ptables += 1

        if ptables == 2:
            clrpage = r'\clearpage'
        else:
            clrpage = ''
        dbval = collect(cnx,plotorder,basins,start_date,start_date,'swe_z',run_name,dfindt,'end')
        dbval = dbval.rename(columns=self.labels)
        start_swe = dbval.round(decimals=dpts)
        dbval = collect(cnx,plotorder,basins,start_date,end_date,'swe_z',run_name,dfindt,'end')
        dbval = dbval.rename(columns=self.labels)
        end_swe = dbval.round(decimals=dpts)
        dswe_byelev = end_swe - start_swe
        dswe_byelev.rename(index={'total':'mean'},inplace=True)
        dswe_byelev.index.name = 'Elevation'
        variables['DSWE_BYELEV'] = (r'  \normalsize \textbf{Change in SWE [%s], %s to %s}\\ \vspace{0.1cm} \\'
                                    %(self.depthlbl,self.report_start.date().strftime("%Y-%-m-%-d"),
                                      self.report_date.date().strftime("%Y-%-m-%-d")) + spacecmd +
                                      dswe_byelev.to_latex(na_rep='-', column_format=colstr) +
                                      r'} \\ \footnotesize{\textbf{Table %s:} Change in SWE.} %s'%(str(mtables),clrpage))
        mtables += 1

    else:
        variables['DSWE_BYELEV'] = ''

    if 'swe_percent' in self.tables:
        ptables += 1
        if (ptables % 2) == 0 and ptables != 0:
            clrpage = r'\clearpage'
        else:
            clrpage = ''

        dbval = collect(cnx,plotorder,basins,start_date,end_date,'swe_vol',run_name,dfind,'end')
        dbval = dbval.rename(columns=self.labels)
        swe_byelev = dbval.round(decimals=dpts)
        value = swe_byelev.iloc[:-1].sum()
        sweper_byelev = (swe_byelev/value*100).round(decimals = dpts)
        sweper_byelev.index.name = 'Elevation'

        variables['SWEPER_BYELEV'] = (r'  \normalsize \textbf{SWE volume, percent of basin total, %s}\\ \vspace{0.1cm} \\'
                                    %(self.report_date.date().strftime("%Y-%-m-%-d"))+
                                    spacecmd + sweper_byelev.round(1).to_latex(na_rep='-', column_format=colstr) +
                                    r'}  \\ \footnotesize{\textbf{Table %s:} Percent of total SWE volume.}%s'%(str(mtables),clrpage))
        variables['SWEPER_BYELEV'] = variables['SWEPER_BYELEV'].replace('inf','-')

        mtables += 1

    else:
        variables['SWEPER_BYELEV'] = ''

    if 'swi_vol' in self.tables:
        dbval = collect(cnx,plotorder,basins,start_date,end_date,'swi_vol',run_name,dfindt,'sum')
        dbval = dbval.rename(columns=self.labels)
        swi_byelev = dbval.round(decimals=dpts)
        variables['ACCUM_BYELEV'] = (r' \normalsize \textbf{SWI [%s] by elevation, %s to %s}\\ \vspace{0.1cm} \\'
                                    %(self.vollbl,self.report_start.date().strftime("%Y-%-m-%-d"),
                                      self.report_date.date().strftime("%Y-%-m-%-d")) + spacecmd +
                                      swi_byelev.to_latex(na_rep='-', column_format=colstr) +
                                      r'} \\ \footnotesize{\textbf{Table %s:} SWI volume. }'%(str(mtables)))
        mtables += 1

    else:
        variables['ACCUM_BYELEV'] = ''

    variables['TOT_LBL'] = self.plotorder[0]

    # for n in range(1,len(self.plotorder)):
    for n in range(1,len(self.labels)):
        s = 'SUB' + str(n) + '_LBL'
        variables[s] = self.labels[self.plotorder[n]]

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
                    'CHANGES_FIG_TPL':os.path.join(self.figs_tpl_path, 'changes_fig_tpl.txt'),
                    'SWI_FIG_TPL':os.path.join(self.figs_tpl_path, 'swi_fig_tpl.txt'),
                    'TOTALS_FIG_TPL':os.path.join(self.figs_tpl_path, 'totals_fig_tpl.txt'),
                    'MULTITOTSWE_FIG_TPL':os.path.join(self.figs_tpl_path, 'multitotswe_fig_tpl.txt'),
                    'VALID_FIG_TPL':os.path.join(self.figs_tpl_path, 'valid_fig_tpl.txt'),
                    'FLTCHANGES_FIG_TPL':os.path.join(self.figs_tpl_path, 'flt_fig_tpl.txt'),
                    'PDEP_FIG_TPL':os.path.join(self.figs_tpl_path, 'pdep_fig_tpl.txt'),
                    'COLD_FIG_TPL':os.path.join(self.figs_tpl_path, 'cold_fig_tpl.txt'),
                    'SWE_FIG_TPL':os.path.join(self.figs_tpl_path, 'swe_fig_tpl.txt'),
                    'SUBBASINS_FIG_TPL':os.path.join(self.figs_tpl_path, 'subbasins_fig_tpl.txt'),
                    'SWEFORECAST_FIG_TPL':os.path.join(self.figs_tpl_path, 'forecastswe_fig_tpl.txt'),
                    'SWIFORECAST_FIG_TPL':os.path.join(self.figs_tpl_path, 'forecastswi_fig_tpl.txt'),
                    'CHANGESFORECAST_FIG_TPL':os.path.join(self.figs_tpl_path, 'forecastchanges_fig_tpl.txt'),
                    'TOTALSFORECAST_FIG_TPL':os.path.join(self.figs_tpl_path, 'forecasttotals_fig_tpl.txt'),
                    'PDEPFORECAST_FIG_TPL':os.path.join(self.figs_tpl_path, 'forecastpdep_fig_tpl.txt'),
                    'DIAGNOSTICS_FIG_TPL':os.path.join(self.figs_tpl_path, 'diagnostics_fig_tpl.txt')
                    }

    # Define and load summary tables depending on number of subbasins
    section_dict['SWE_SUMMARY_TPL'] = os.path.join(self.figs_tpl_path,
        'swe_summary_{}sub.txt'.format(str(len(self.plotorder))))

    # Remove if no flight
    if not self.flt_flag or not self.flight_figs:
        del section_dict['FLTCHANGES_FIG_TPL']

    if not self.report_diagnostics:
        del section_dict['DIAGNOSTICS_FIG_TPL']
        variables['DIAGNOSTICS_FIG'] = ''

    # Remove if no forecast
    if not self.forecast_flag:
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

    if not self.subs_fig:
        variables['SUBBASINS_FIG_TPL'] = ''

    if not self.rep_compare_runs_flag:
        variables['MULTITOTSWE_FIG_TPL'] = ''

    if not self.rep_swi_flag:
        variables['SWI_FIG_TPL'] = ''

    if not self.rep_image_change_flag:
        variables['CHANGES_FIG_TPL'] = ''

    if not self.rep_cold_content_flag:
        variables['COLD_FIG_TPL'] = ''

    if not self.rep_swe_volume_flag:
        variables['SWE_FIG_TPL'] = ''

    if not self.rep_basin_total_flag:
        variables['TOTALS_FIG_TPL'] = ''

    if not self.rep_stn_validate_flag:
        variables['VALID_FIG_TPL'] = ''

    if not self.rep_compare_runs_flag:
        variables['MULTITOTSWE_FIG_TPL'] = ''

    if not self.rep_precip_depth_flag:
        variables['PDEP_FIG_TPL'] = ''

    # Make the report
    env = make_env(loader = FileSystemLoader(self.templ_path))
    tpl = env.get_template(self.tex_file)

    # print VAR-replaced latex file for debugging if desired
    if self.print_latex:
        print(tpl.render(variables))

    pdf = build_pdf(tpl.render(variables))

    # Save in reports and with figs
    rpath = os.path.join(self.figs_path, '' + self.report_name)
    pdf.save_to(rpath)
    self._logger.info(' Saved {}'.format(rpath))

    if self.rep_path is not None:
        rpath = os.path.join(self.rep_path, '' + self.report_name)
        pdf.save_to(rpath)
        self._logger.info(' Saved {}'.format(rpath))
