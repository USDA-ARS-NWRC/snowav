from datetime import datetime
from jinja2 import FileSystemLoader
from latex.jinja2 import make_env
from latex import build_pdf
import os

import snowav
from snowav.database.database import collect


def report(cfg):
    """ Create the pdf report.

    The latex formats for figures, tables, and summary paragraphs are created
    by templates in snowav/report/. Data is pulled from the snowav database.

    To add a figure to the report:
    - create the plot (of course)
    - add to variables dict in this file, with the same naming convention
        variables['PDEP_FIG'] = 'precip_depth%s.png'%(cfg.directory)
    - add figure template to snowav_report.tex
        \VAR{PDEP_FIG_TPL}
    - add to section_dict
    - add template to snowav/report/figs/
    - add to snowav/config/CoreConfig.py [report] exclude_figs
    - if the figure may not exist (such as flt_diff, or those that require
      process() to run), address that with some form of exception before
      creating the pdf

    Args
    ------
    cfg {class}: config class
    """

    # bid = cfg.basins[cfg.plotorder[0]]['basin_id']
    basins = cfg.basins
    wy_start = datetime(cfg.wy-1, 10, 1)
    start_date = cfg.start_date
    end_date = cfg.end_date
    run_name = cfg.run_name
    # edges = cfg.edges
    plotorder = cfg.plotorder
    dpts = int(cfg.dplcs)
    ddpts = int(cfg.rep_dplcs)
    cnx = cfg.connector

    # variables to pass to latex file
    variables = {}

    dbval = collect(cnx, plotorder[0], basins, wy_start, end_date, 'swi_vol',
                    run_name, 'total', 'sum')
    variables['TOTAL_SWI'] = dbval.sum().values.round(decimals=dpts)[0]

    dbval = collect(cnx, plotorder[0], basins, start_date, end_date, 'swi_vol',
                    run_name, 'total', 'sum')
    variables['PER_SWI'] = dbval.sum().values.round(decimals=dpts)[0]

    dbval = collect(cnx, plotorder[0], basins, start_date, end_date, 'swe_vol',
                    run_name, 'total', 'end')
    variables['TOTAL_SWE'] = dbval.sum().values.round(decimals=dpts)[0]

    dbval = collect(cnx, plotorder[0], basins, start_date, end_date, 'swe_avail',
                    run_name, 'total', 'end')
    variables['TOTAL_SAV'] = dbval.sum().values.round(decimals=dpts)[0]

    dbval = collect(cnx, plotorder[0], basins, start_date, end_date, 'swe_z',
                    run_name, 'total', 'end')
    variables['TOTAL_PM'] = dbval.sum().values.round(decimals=ddpts)[0]

    dbval = collect(cnx, plotorder[0], basins, wy_start, end_date, 'precip_z',
                    run_name, 'total', 'sum')
    variables['TOTALPRE_PM'] = dbval.sum().values.round(decimals=ddpts)[0]

    s = collect(cnx, plotorder[0], basins, start_date, start_date, 'swe_vol',
                run_name, 'total', 'end')
    e = collect(cnx, plotorder[0], basins, start_date, end_date, 'swe_vol',
                run_name, 'total', 'end')
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

    dbval = collect(cnx, plotorder[0], basins, wy_start, end_date, 'rain_z',
                    run_name, 'total', 'sum')
    total_rain = dbval.sum().values.round(decimals=ddpts)[0]

    if variables['TOTALPRE_PM'] != 0:
        variables['TOTAL_RAT'] = str(int((total_rain / variables['TOTALPRE_PM']) * 100))
    else:
        variables['TOTAL_RAT'] = '-'
        variables['TOTALPRE_PM'] = '-'

    report_time = datetime.now().strftime("%Y-%-m-%-d %H:%M")

    numsubs = range(1, len(cfg.plotorder))

    for n, sub in zip(numsubs, cfg.plotorder[1:]):
        swiind = 'SUB' + str(n) + '_SWI'
        perswiind = 'SUB' + str(n) + '_PERSWI'
        sweind = 'SUB' + str(n) + '_SWE'
        avsweind = 'SUB' + str(n) + '_SAV'
        swedel = 'SUB' + str(n) + '_SDEL'
        pm = 'SUB' + str(n) + '_PM'
        prepm = 'SUB' + str(n) + 'PRE_PM'
        rain = 'SUB' + str(n) + 'RAI'
        ratio = 'SUB' + str(n) + '_RAT'
        pvol = 'SUB' + str(n) + '_PVOL'

        dbval = collect(cnx, sub, basins, wy_start, end_date, 'swi_vol',
                        run_name, 'total', 'sum')
        variables[swiind] = dbval.sum().values.round(decimals=dpts)[0]

        dbval = collect(cnx, sub, basins, start_date, end_date, 'swi_vol',
                        run_name, 'total', 'sum')
        variables[perswiind] = dbval.sum().values.round(decimals=dpts)[0]

        dbval = collect(cnx, sub, basins, start_date, end_date, 'swe_vol',
                        run_name, 'total', 'end')
        variables[sweind] = dbval.sum().values.round(decimals=dpts)[0]

        dbval = collect(cnx, sub, basins, start_date, end_date, 'swe_avail',
                        run_name, 'total', 'end')
        variables[avsweind] = dbval.sum().values.round(decimals=dpts)[0]

        dbval = collect(cnx, sub, basins, start_date, start_date, 'swe_vol',
                        run_name, 'total', 'end')
        start_swe = dbval.sum().values.round(decimals=dpts)[0]
        dbval = collect(cnx, sub, basins, start_date, end_date, 'swe_vol',
                        run_name, 'total', 'end')
        end_swe = dbval.sum().values.round(decimals=dpts)[0]
        variables[swedel] = end_swe - start_swe

        dbval = collect(cnx, sub, basins, start_date, end_date, 'swe_z',
                        run_name, 'total', 'end')
        variables[pm] = dbval.sum().values.round(decimals=ddpts)[0]

        dbval = collect(cnx, sub, basins, wy_start, end_date, 'precip_z',
                        run_name, 'total', 'sum')
        variables[prepm] = dbval.sum().values.round(decimals=ddpts)[0]

        dbval = collect(cnx, sub, basins, wy_start, end_date, 'rain_z',
                        run_name, 'total', 'sum')
        variables[rain] = dbval.sum().values.round(decimals=ddpts)[0]

        if (end_swe > 0) and (variables['TOTAL_SWE'] > 0):
            variables[pvol] = end_swe / variables['TOTAL_SWE'] * 100
        else:
            variables[pvol] = '-'

        if variables[prepm] != 0.0:
            variables[ratio] = str(int((variables[rain] / variables[prepm]) * 100))
        else:
            variables[ratio] = '0'

    # Upper case variables are used in the LaTex file,
    # lower case versions are assigned here

    # untested - if report title contains comma?
    if isinstance(cfg.rep_title, list):
        title = cfg.rep_title[0]
        for s in cfg.rep_title[1::]:
            title = title + ', ' + s
        variables['REPORT_TITLE'] = title

    else:
        if cfg.flt_flag and cfg.flight_figs:
            fst = ' Model snow depths were updated with ASO snow depths on '
            tst = ''
            for i, d in enumerate(cfg.flight_diff_dates):
                if d <= cfg.end_date:
                    dn = d.date().strftime("%m/%d")
                    if len(cfg.flight_diff_dates) == 1:
                        fst = fst + dn + '.'
                        tst += dn

                    if ((len(cfg.flight_diff_dates) > 1) and
                            (i < len(cfg.flight_diff_dates) - 1)):
                        fst = fst + dn + ', '
                        tst += dn + ', '

                    if ((len(cfg.flight_diff_dates) > 1) and
                            (i == len(cfg.flight_diff_dates) - 1)):
                        fst = fst + 'and ' + dn + '.'
                        tst += dn

                else:
                    fst = fst.split(dn)[0] + 'and ' + dn + '.'
                    break

            variables['REPORT_TITLE'] = (cfg.rep_title +
                                         r' \\ ASO Updates {}'.format(tst))
            variables['FLTSENT'] = fst

        else:
            variables['REPORT_TITLE'] = cfg.rep_title
            variables['FLTSENT'] = ''

    variables['REPORT_TIME'] = report_time
    variables['WATERYEAR'] = str(cfg.wy)
    variables['UNITS'] = cfg.vollbl
    variables['VOLLBL'] = cfg.vollbl
    variables['DEPLBL'] = cfg.depthlbl
    variables['START_DATE'] = cfg.report_start.date().strftime("%B %-d")
    variables['END_DATE'] = cfg.report_date.date().strftime("%B %-d")
    variables['SWE_IN'] = variables['TOTAL_PM']
    variables['SWI_IN'] = variables['TOTAL_SWI']
    variables['FIG_PATH'] = cfg.figs_path
    variables['SWI_FIG'] = 'swi_{}.png'.format(cfg.directory)
    variables['CHANGES_FIG'] = 'swe_change_{}.png'.format(cfg.directory)
    variables['TOTALS_FIG'] = 'basin_total_{}.png'.format(cfg.directory)
    variables['MULTITOTSWE_FIG'] = 'compare_swe_vol_{}.png'.format(cfg.directory)
    variables['DENSITY_FIG'] = 'density_{}.png'.format(cfg.directory)
    variables['PDEP_FIG'] = 'precip_depth_{}.png'.format(cfg.directory)
    variables['VALID_FIG'] = cfg.stn_validate_fig_name
    variables['COLD_FIG'] = 'cold_content_{}.png'.format(cfg.directory)
    variables['SWE_FIG'] = 'swe_volume_{}.png'.format(cfg.directory)
    variables['VERSION'] = snowav.__version__

    if (cfg.update_file is not None) and cfg.flt_flag and cfg.flight_figs:
        for name in cfg.flight_diff_fig_names:
            variables['DFLT_FIG'] = name

    if cfg.forecast_flag:
        variables['FORE_START_DATE'] = cfg.for_start_date.date().strftime("%B %-d")
        variables['FORE_DATE'] = cfg.for_end_date.date().strftime("%B %-d")
        variables['SWEFORECAST_FIG'] = 'swe_volume_{}.png'.format(cfg.directory + '_forecast')
        variables['SWIFORECAST_FIG'] = 'swi_{}.png'.format(cfg.directory + '_forecast')
        variables['CHANGESFORECAST_FIG'] = 'swe_change_{}.png'.format(cfg.directory + '_forecast')
        variables['TOTALSFORECAST_FIG'] = 'basin_total_{}.png'.format(cfg.directory + '_forecast')
        variables['PDEPFORECAST_FIG'] = 'precip_depth_{}.png'.format(cfg.directory + '_forecast')

    if cfg.report_diagnostics:
        variables['DIAGNOSTICS_FIG'] = 'diagnostics_{}'.format(cfg.directory)
        variables['INPUTS_FIG'] = 'inputs_period_{}'.format(cfg.directory)

    if cfg.subs_fig is not None:
        variables['SUBBASINS_FIG'] = '{}'.format(cfg.subs_fig)

    # Logos
    variables['ARSLOGO'] = os.path.join(cfg.figs_tpl_path, 'ARS.jpg')
    variables['ASOLOGO'] = os.path.join(cfg.figs_tpl_path, 'ASO.jpg')
    variables['USDALOGO'] = os.path.join(cfg.figs_tpl_path, 'USDA.png')
    variables['JPLLOGO'] = os.path.join(cfg.figs_tpl_path, 'JPL.jpg')
    variables['CDWRLOGO'] = os.path.join(cfg.figs_tpl_path, 'CDWR.png')
    variables['USBRLOGO'] = os.path.join(cfg.figs_tpl_path, 'USBR.jpg')
    variables['NRCSLOGO'] = os.path.join(cfg.figs_tpl_path, 'NRCS.jpg')
    variables['KRWALOGO'] = os.path.join(cfg.figs_tpl_path, 'KRWA.jpg')
    variables['FRIANTLOGO'] = os.path.join(cfg.figs_tpl_path, 'FRIANT.jpg')
    variables['AWSMLOGO'] = os.path.join(cfg.figs_tpl_path, 'logo.png')

    dfind = [str(i) for i in cfg.edges]
    dfindt = [str(i) for i in cfg.edges] + ['total']
    colstr = 'l' + 'r' * len(cfg.plotorder)

    if len(cfg.plotorder) > 5:
        spacecmd = r'\resizebox{\textwidth}{!}{'
    else:
        spacecmd = r'{'

    # ntables = len(cfg.tables)
    mtables = 2
    ptables = 0

    if 'swe_depth' in cfg.tables:

        dbval = collect(cnx, plotorder, basins, start_date, end_date, 'swe_z', run_name, dfindt, 'end')
        dbval = dbval.rename(columns=cfg.labels)
        swe_byelev = dbval.round(decimals=dpts)
        swe_byelev.rename(index={'total': 'mean'}, inplace=True)
        swe_byelev.index.name = 'Elevation'
        variables['SWE_BYELEV'] = (r' \normalsize \textbf{SWE [%s], %s}\\ \vspace{0.1cm} \\'
                                   % (cfg.depthlbl, cfg.report_date.date().strftime("%Y-%-m-%-d")) +
                                   spacecmd + swe_byelev.to_latex(na_rep='-', column_format=colstr) +
                                   r'} \\ \footnotesize{\textbf{Table %s:} SWE depth.}' % (str(mtables)))
        mtables += 1
        ptables += 1

    else:
        variables['SWE_BYELEV'] = ''

    if 'swe_vol' in cfg.tables:
        ptables += 1

        if ptables == 2:
            clrpage = r'\clearpage'
        else:
            clrpage = ''
        dbval = collect(cnx, plotorder, basins, start_date, end_date, 'swe_vol', run_name, dfindt, 'end')
        dbval = dbval.rename(columns=cfg.labels)
        swe_byelev = dbval.round(decimals=dpts)
        swe_byelev.index.name = 'Elevation'
        variables['SWEVOL_BYELEV'] = (r' \normalsize \textbf{SWE [%s], %s}\\ \vspace{0.1cm} \\'
                                      % (cfg.vollbl, cfg.report_date.date().strftime("%Y-%-m-%-d")) +
                                      spacecmd + swe_byelev.to_latex(na_rep='-', column_format=colstr) +
                                      r'} \\ \footnotesize{\textbf{Table %s:} SWE volume.}%s' % (str(mtables), clrpage))
        mtables += 1

    else:
        variables['SWEVOL_BYELEV'] = ''

    if 'swe_change' in cfg.tables:
        ptables += 1

        if ptables == 2:
            clrpage = r'\clearpage'
        else:
            clrpage = ''
        dbval = collect(cnx, plotorder, basins, start_date, start_date, 'swe_z', run_name, dfindt, 'end')
        dbval = dbval.rename(columns=cfg.labels)
        start_swe = dbval.round(decimals=dpts)
        dbval = collect(cnx, plotorder, basins, start_date, end_date, 'swe_z', run_name, dfindt, 'end')
        dbval = dbval.rename(columns=cfg.labels)
        end_swe = dbval.round(decimals=dpts)
        dswe_byelev = end_swe - start_swe
        dswe_byelev.rename(index={'total': 'mean'}, inplace=True)
        dswe_byelev.index.name = 'Elevation'
        variables['DSWE_BYELEV'] = (r'  \normalsize \textbf{Change in SWE [%s], %s to %s}\\ \vspace{0.1cm} \\'
                                    % (cfg.depthlbl, cfg.report_start.date().strftime("%Y-%-m-%-d"),
                                       cfg.report_date.date().strftime("%Y-%-m-%-d")) + spacecmd +
                                    dswe_byelev.to_latex(na_rep='-', column_format=colstr) +
                                    r'} \\ \footnotesize{\textbf{Table %s:} Change in SWE.} %s' % (
                                    str(mtables), clrpage))
        mtables += 1

    else:
        variables['DSWE_BYELEV'] = ''

    if 'swe_percent' in cfg.tables:
        ptables += 1
        if (ptables % 2) == 0 and ptables != 0:
            clrpage = r'\clearpage'
        else:
            clrpage = ''

        dbval = collect(cnx, plotorder, basins, start_date, end_date, 'swe_vol', run_name, dfind, 'end')
        dbval = dbval.rename(columns=cfg.labels)
        swe_byelev = dbval.round(decimals=dpts)
        value = swe_byelev.iloc[:-1].sum()
        sweper_byelev = (swe_byelev / value * 100).round(decimals=dpts)
        sweper_byelev.index.name = 'Elevation'

        variables['SWEPER_BYELEV'] = (
                    r'  \normalsize \textbf{SWE volume, percent of basin total, %s}\\ \vspace{0.1cm} \\'
                    % (cfg.report_date.date().strftime("%Y-%-m-%-d")) +
                    spacecmd + sweper_byelev.round(1).to_latex(na_rep='-', column_format=colstr) +
                    r'}  \\ \footnotesize{\textbf{Table %s:} Percent of total SWE volume.}%s' % (str(mtables), clrpage))
        variables['SWEPER_BYELEV'] = variables['SWEPER_BYELEV'].replace('inf', '-')

        mtables += 1

    else:
        variables['SWEPER_BYELEV'] = ''

    if 'swi_vol' in cfg.tables:
        dbval = collect(cnx, plotorder, basins, start_date, end_date, 'swi_vol', run_name, dfindt, 'sum')
        dbval = dbval.rename(columns=cfg.labels)
        swi_byelev = dbval.round(decimals=dpts)
        variables['ACCUM_BYELEV'] = (r' \normalsize \textbf{SWI [%s] by elevation, %s to %s}\\ \vspace{0.1cm} \\'
                                     % (cfg.vollbl, cfg.report_start.date().strftime("%Y-%-m-%-d"),
                                        cfg.report_date.date().strftime("%Y-%-m-%-d")) + spacecmd +
                                     swi_byelev.to_latex(na_rep='-', column_format=colstr) +
                                     r'} \\ \footnotesize{\textbf{Table %s:} SWI volume. }' % (str(mtables)))
        mtables += 1

    else:
        variables['ACCUM_BYELEV'] = ''

    variables['TOT_LBL'] = cfg.plotorder[0]

    # for n in range(1,len(cfg.plotorder)):
    for n in range(1, len(cfg.labels)):
        s = 'SUB' + str(n) + '_LBL'
        variables[s] = cfg.labels[cfg.plotorder[n]]

    # Convert floats to strings
    for name in variables:
        if isinstance(variables[name], float):
            if cfg.dplcs == 0:
                tmp = str(int(variables[name]))
            else:
                tmp = str(round(variables[name], cfg.dplcs))
            variables[name] = tmp

    # Summary sections and fig template have variable strings
    # (e.g. CHANGES_FIG) that need to be replaced
    section_dict = {'SUMMARY': cfg.summary_file,
                    'CHANGES_FIG_TPL': os.path.join(cfg.figs_tpl_path, 'changes_fig_tpl.txt'),
                    'SWI_FIG_TPL': os.path.join(cfg.figs_tpl_path, 'swi_fig_tpl.txt'),
                    'TOTALS_FIG_TPL': os.path.join(cfg.figs_tpl_path, 'totals_fig_tpl.txt'),
                    'MULTITOTSWE_FIG_TPL': os.path.join(cfg.figs_tpl_path, 'multitotswe_fig_tpl.txt'),
                    'VALID_FIG_TPL': os.path.join(cfg.figs_tpl_path, 'valid_fig_tpl.txt'),
                    'FLTCHANGES_FIG_TPL': os.path.join(cfg.figs_tpl_path, 'flt_fig_tpl.txt'),
                    'PDEP_FIG_TPL': os.path.join(cfg.figs_tpl_path, 'pdep_fig_tpl.txt'),
                    'COLD_FIG_TPL': os.path.join(cfg.figs_tpl_path, 'cold_fig_tpl.txt'),
                    'SWE_FIG_TPL': os.path.join(cfg.figs_tpl_path, 'swe_fig_tpl.txt'),
                    'SUBBASINS_FIG_TPL': os.path.join(cfg.figs_tpl_path, 'subbasins_fig_tpl.txt'),
                    'SWEFORECAST_FIG_TPL': os.path.join(cfg.figs_tpl_path, 'forecastswe_fig_tpl.txt'),
                    'SWIFORECAST_FIG_TPL': os.path.join(cfg.figs_tpl_path, 'forecastswi_fig_tpl.txt'),
                    'CHANGESFORECAST_FIG_TPL': os.path.join(cfg.figs_tpl_path, 'forecastchanges_fig_tpl.txt'),
                    'TOTALSFORECAST_FIG_TPL': os.path.join(cfg.figs_tpl_path, 'forecasttotals_fig_tpl.txt'),
                    'PDEPFORECAST_FIG_TPL': os.path.join(cfg.figs_tpl_path, 'forecastpdep_fig_tpl.txt'),
                    'DIAGNOSTICS_FIG_TPL': os.path.join(cfg.figs_tpl_path, 'diagnostics_fig_tpl.txt')
                    }

    # Define and load summary tables depending on number of subbasins
    section_dict['SWE_SUMMARY_TPL'] = os.path.join(cfg.figs_tpl_path,
                                                   'swe_summary_{}sub.txt'.format(str(len(cfg.plotorder))))

    # Remove if no flight
    if not cfg.flt_flag or not cfg.flight_figs:
        del section_dict['FLTCHANGES_FIG_TPL']

    if not cfg.report_diagnostics:
        del section_dict['DIAGNOSTICS_FIG_TPL']
        variables['DIAGNOSTICS_FIG'] = ''

    # Remove if no forecast
    if not cfg.forecast_flag:
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
        variables['FORE_TITLE'] = 'with WRF forecast, {} to {}'.format(cfg.for_start_date.date().strftime("%b %-d"),
                                                                       cfg.for_end_date.date().strftime("%b %-d"))

    for rep in section_dict.keys():

        # for variable numbers of fligth figures
        if rep == 'FLTCHANGES_FIG_TPL':
            fid = open(section_dict[rep], 'r')
            var = fid.read()
            fid.close()

            for name in sorted(variables):
                if name == 'DFLT_FIG':
                    for i, fltname in enumerate(cfg.flight_diff_fig_names):
                        flt_date = cfg.flight_outputs['dates'][i].date().strftime("%Y%m%d")
                        flt_num = cfg.flight_delta_vol_df[flt_date].round(1).to_latex(na_rep='-', column_format=colstr)
                        table = (r'{ \vspace{0.5cm} \textbf{Change in SWE [%s] by elevation, ' % (cfg.vollbl) +
                                 r'from %s update} \\ \vspace{0.1cm} \\' % (flt_date) +
                                 r' %s %s }' % (spacecmd, flt_num))

                        if i == 0:
                            tmp = var.replace(name, fltname)
                        else:
                            tmp += var.replace(name, fltname)

                        tmp += table

                    var = tmp.replace(name, variables[name])

                else:
                    var = var.replace(name, variables[name])

            variables[rep] = var

        else:
            fid = open(section_dict[rep], 'r')
            var = fid.read()
            fid.close()

            for name in sorted(variables):
                var = var.replace(name, variables[name])
            variables[rep] = var

    if not cfg.subs_fig:
        variables['SUBBASINS_FIG_TPL'] = ''

    if not cfg.rep_compare_runs_flag:
        variables['MULTITOTSWE_FIG_TPL'] = ''

    if not cfg.rep_swi_flag:
        variables['SWI_FIG_TPL'] = ''

    if not cfg.rep_image_change_flag:
        variables['CHANGES_FIG_TPL'] = ''

    if not cfg.rep_cold_content_flag:
        variables['COLD_FIG_TPL'] = ''

    if not cfg.rep_swe_volume_flag:
        variables['SWE_FIG_TPL'] = ''

    if not cfg.rep_basin_total_flag:
        variables['TOTALS_FIG_TPL'] = ''

    if not cfg.rep_stn_validate_flag:
        variables['VALID_FIG_TPL'] = ''

    if not cfg.rep_compare_runs_flag:
        variables['MULTITOTSWE_FIG_TPL'] = ''

    if not cfg.rep_precip_depth_flag:
        variables['PDEP_FIG_TPL'] = ''

    # Make the report
    env = make_env(loader=FileSystemLoader(cfg.templ_path))
    tpl = env.get_template(cfg.tex_file)

    # print VAR-replaced latex file for debugging if desired
    if cfg.print_latex:
        print(tpl.render(variables))

    pdf = build_pdf(tpl.render(variables))

    # Save in reports and with figs
    rpath = os.path.join(cfg.figs_path, '' + cfg.report_name)
    pdf.save_to(rpath)
    cfg._logger.info(' Saved {}'.format(rpath))

    if cfg.rep_path is not None:
        rpath = os.path.join(cfg.rep_path, '' + cfg.report_name)
        pdf.save_to(rpath)
        cfg._logger.info(' Saved {}'.format(rpath))
