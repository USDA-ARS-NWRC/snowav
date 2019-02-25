
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


def report(obj):
    '''
    This function makes the pdf report.

    The latex formats for figures, tables, and summary paragraphs are created
    by templates in snowav/report/. Data is now pulled from the SNOWAV
    database.

    To add a figure to the report:
    - create the plot (of course)
    - add to variables dict in this file, with the same naming convention
        variables['PDEP_FIG'] = 'precip_depth%s.png'%(obj.name_append)
    - add figure template to snowav_report.tex
        \VAR{PDEP_FIG_TPL}
    - add to section_dict in this file
    - add template to snowav/report/figs/
    - add to snowav/config/CoreConfig.py [report] exclude_figs
    - if the figure may not exist (such as flt_diff, or those that require
      process() to run), address that with some form of exception before
      creating the pdf

    '''

    bid = Basins.basins[obj.plotorder[0]]['basin_id']
    wy_start = datetime(obj.wy-1,10,1)
    start_date = obj.start_date
    end_date = obj.end_date

    r = database.database.query(obj, wy_start, obj.end_date, obj.run_name,
                                bid = obj.plotorder)

    # Initialize variables to pass to latex file
    variables = {}
    variables['TOTAL_SWI'] = r[ (r['basin_id']==bid)
                                & (r['date_time']>=wy_start)
                                & (r['date_time']<=end_date)
                                & (r['elevation']=='total')
                                & (r['variable']=='swi_vol')]['value'].sum().round(decimals = int(obj.dplcs))
    variables['PER_SWI'] = r[ (r['basin_id']==bid)
                                & (r['date_time']>=start_date)
                                & (r['date_time']<=end_date)
                                & (r['elevation']=='total')
                                & (r['variable']=='swi_vol')]['value'].sum().round(decimals = int(obj.dplcs))
    variables['TOTAL_SWE'] = r[ (r['basin_id']==bid)
                                & (r['date_time']==end_date)
                                & (r['elevation']=='total')
                                & (r['variable']=='swe_vol')]['value'].values[0].round(decimals = int(obj.dplcs))
    variables['TOTAL_SAV'] = r[  (r['basin_id']==bid)
                                & (r['date_time']==end_date)
                                & (r['elevation']=='total')
                                & (r['variable']=='swe_avail')]['value'].values[0].round(decimals = int(obj.dplcs))
    variables['TOTAL_PVOL'] = '100'

    start_swe = r[ (r['basin_id']==bid)
                    & (r['date_time']==start_date)
                    & (r['elevation']=='total')
                    & (r['variable']=='swe_vol')]['value'].values[0].round(decimals = int(obj.dplcs))
    end_swe = r[ (r['basin_id']==bid)
                 & (r['date_time']==end_date)
                 & (r['elevation']=='total')
                 & (r['variable']=='swe_vol')]['value'].values[0].round(decimals = int(obj.dplcs))

    variables['TOTAL_SDEL'] = end_swe - start_swe

    if float(end_swe) - float(start_swe) > 0:
        variables['SIGN'] = r'$+$'

    if float(end_swe) - float(start_swe) == 0.0:
        variables['SIGN'] = ''

    if float(end_swe) - float(start_swe) < 0.0:
        variables['SIGN'] = r'-'

    variables['TOTAL_PM'] = r[ (r['basin_id']==bid)
                                & (r['date_time']==end_date)
                                & (r['elevation']=='total')
                                & (r['variable']=='swe_z')]['value'].values[0].round(decimals = int(obj.dplcs))
    variables['TOTALPRE_PM'] = r[ (r['basin_id']==bid)
                                & (r['date_time']>=wy_start)
                                & (r['date_time']<=end_date)
                                & (r['elevation']=='total')
                                & (r['variable']=='precip_z')]['value'].sum().round(decimals = int(obj.dplcs))
    total_rai = r[ (r['basin_id']==bid)
                    & (r['date_time']>=wy_start)
                    & (r['date_time']<=end_date)
                    & (r['elevation']=='total')
                    & (r['variable']=='rain_z')]['value'].sum().round(decimals = int(obj.dplcs))

    if variables['TOTALPRE_PM'] != 0.0:
        variables['TOTAL_RAT'] = str(int((total_rai/variables['TOTALPRE_PM'])*100))
    else:
        variables['TOTAL_RAT'] = '0'

    report_time = datetime.now().strftime("%Y-%-m-%-d %H:%M")
    if hasattr(obj,'orig_date'):
        report_time = obj.orig_date + ', revised ' + report_time

    numsubs = range(1,len(obj.plotorder))

    for n,sub in zip(numsubs,obj.plotorder[1:]):
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

        swival = r[ (r['basin_id']==Basins.basins[sub]['basin_id'])
                     & (r['date_time']>=wy_start)
                     & (r['date_time']<=end_date)
                     & (r['elevation']=='total')
                     & (r['variable']=='swi_vol')]['value'].sum()
        perswival = r[ (r['basin_id']==Basins.basins[sub]['basin_id'])
                     & (r['date_time']>=start_date)
                     & (r['date_time']<=end_date)
                     & (r['elevation']=='total')
                     & (r['variable']=='swi_vol')]['value'].sum()
        sweval = r[ (r['basin_id']==Basins.basins[sub]['basin_id'])
                     & (r['date_time']==end_date)
                     & (r['elevation']=='total')
                     & (r['variable']=='swe_vol')]['value'].values[0]
        avsweval = r[ (r['basin_id']==Basins.basins[sub]['basin_id'])
                     & (r['date_time']==end_date)
                     & (r['elevation']=='total')
                     & (r['variable']=='swe_avail')]['value'].values[0]
        start_swe = r[ (r['basin_id']==Basins.basins[sub]['basin_id'])
                        & (r['date_time']==start_date)
                        & (r['elevation']=='total')
                        & (r['variable']=='swe_vol')]['value'].values[0]
        end_swe = r[ (r['basin_id']==Basins.basins[sub]['basin_id'])
                     & (r['date_time']==end_date)
                     & (r['elevation']=='total')
                     & (r['variable']=='swe_vol')]['value'].values[0]
        swedelval = end_swe - start_swe
        pmval = r[ (r['basin_id']==Basins.basins[sub]['basin_id'])
                    & (r['date_time']==end_date)
                    & (r['elevation']=='total')
                    & (r['variable']=='swe_z')]['value'].values[0]
        prepmval = r[ (r['basin_id']==Basins.basins[sub]['basin_id'])
                    & (r['date_time']>=wy_start)
                    & (r['date_time']<=end_date)
                    & (r['elevation']=='total')
                    & (r['variable']=='precip_z')]['value'].sum()
        rainval = r[ (r['basin_id']==Basins.basins[sub]['basin_id'])
                    & (r['date_time']>=wy_start)
                    & (r['date_time']<=end_date)
                    & (r['elevation']=='total')
                    & (r['variable']=='rain_z')]['value'].sum()

        if (end_swe > 0) and (variables['TOTAL_SWE'] > 0):
            pvol = end_swe/variables['TOTAL_SWE']*100
        else:
            pvol = '-'

        if prepmval != 0.0:
            ratval = str(int((rainval/prepmval)*100))
        else:
            ratval = '0'

        variables[SWIIND] = swival
        variables[PERSWIIND] = perswival
        variables[SWEIND] = sweval
        variables[AVSWEIND] = avsweval
        variables[SWEDEL] = swedelval
        variables[PM] = pmval
        variables[PREPM] = prepmval
        variables[RAIN] = rainval
        variables[RATIO] = ratval
        variables[PVOL] = pvol

    # Upper case variables are used in the LaTex file,
    # lower case versions are assigned here

    # untested - if report title contains comma?
    if isinstance(obj.rep_title, list):
        title = obj.rep_title[0]
        for s in obj.rep_title[1::]:
            title = title + ', ' + s
        variables['REPORT_TITLE'] = title
    else:
        variables['REPORT_TITLE'] = obj.rep_title

    variables['REPORT_TIME'] = report_time
    variables['WATERYEAR'] = str(obj.wy)
    variables['UNITS'] = obj.vollbl
    variables['VOLLBL'] = obj.vollbl
    variables['DEPLBL'] = obj.depthlbl
    variables['START_DATE'] = obj.report_start.date().strftime("%B %-d")
    variables['END_DATE'] = obj.report_date.date().strftime("%B %-d")
    variables['FORE_DATE'] = ''
    variables['SWE_IN'] = variables['TOTAL_PM']
    variables['SWI_IN'] = variables['TOTAL_SWI']
    variables['FIG_PATH'] = obj.figs_path

    # Figures
    variables['SWI_FIG'] = 'swi_%s.png'%(obj.name_append)
    variables['CHANGES_FIG'] = 'swe_change_%s.png'%(obj.name_append)
    variables['DFLT_FIG'] = 'dflt_swe_change_%s.png'%(obj.name_append)
    variables['CHANGES_DEP_FIG'] = 'swe_change_depth_%s.png'%(obj.name_append)
    variables['TOTALS_FIG'] = 'basin_total_%s.png'%(obj.name_append)
    variables['MULTITOTSWE_FIG'] = 'compare_swe_vol_%s.png'%(obj.name_append)
    # variables['MULTITOTSWI_FIG'] = 'compare_swi_vol_%s.png'%(obj.name_append)
    variables['HYP_FIG'] = 'hypsometry_%s.png'%(obj.name_append)
    variables['MEAN_FIG'] = 'mean_swe_depth_%s.png'%(obj.name_append)
    variables['DETAIL_FIG'] = 'mean_detail_%s.png'%(obj.name_append)
    variables['DENSITY_FIG'] = 'density_%s.png'%(obj.name_append)
    variables['DENSITY_SUB_FIG'] = 'density_sub_%s.png'%(obj.name_append)
    variables['DENSITY_SWE_FIG'] = 'density_swe_%s.png'%(obj.name_append)
    variables['PDEP_FIG'] = 'precip_depth_%s.png'%(obj.name_append)
    variables['VALID_FIG'] = 'validation_%s.png'%(obj.name_append)
    variables['COLD_FIG'] = 'cold_content_{}.png'.format(obj.name_append)
    variables['SWE_FIG'] = 'swe_volume_{}.png'.format(obj.name_append)
    variables['INFLOW_FIG'] = 'inflow_{}.png'.format(obj.name_append)

    if obj.subs_fig is not None:
        variables['SUBBASINS_FIG'] = '{}'.format(obj.subs_fig)

    else:
        if obj.exclude_figs != None:
            obj.exclude_figs = obj.exclude_figs + ['SUBBASINS']

        else:
            obj.exclude_figs = ['SUBBASINS']

    # variables['MULTITOT_FIG'] = 'basin_total_multiyr_%s.png'%(obj.name_append)
    # variables['RESULTS_FIG'] = 'results_%s.png'%(obj.name_append)
    # variables['ELEV_FIG'] = 'swe_elev_%s.png'%(obj.name_append)

    # Logos
    variables['ARSLOGO'] = obj.figs_tpl_path + 'ARS.jpg'
    variables['ASOLOGO'] = obj.figs_tpl_path + 'ASO.jpg'
    variables['USDALOGO'] = obj.figs_tpl_path + 'USDA.png'
    variables['JPLLOGO'] = obj.figs_tpl_path + 'JPL.jpg'
    variables['CDWRLOGO'] = obj.figs_tpl_path + 'CDWR.png'
    variables['USBRLOGO'] = obj.figs_tpl_path + 'USBR.jpg'
    variables['NRCSLOGO'] = obj.figs_tpl_path + 'NRCS.jpg'
    variables['KRWALOGO'] = obj.figs_tpl_path + 'KRWA.jpg'
    variables['FRIANTLOGO'] = obj.figs_tpl_path + 'FRIANT.jpg'
    variables['AWSMLOGO'] = obj.figs_tpl_path + 'logo.png'

    # dfind = [str(i) for i in obj.edges] + ['total']
    # dfind = obj.edges
    dfind = [str(i) for i in obj.edges]
    # bmeanstr = 'basin mean'

    swe_byelev = pd.DataFrame(np.nan, index = dfind, columns = obj.plotorder)
    swe_total = pd.DataFrame(np.nan, index = ['basin mean'], columns = obj.plotorder)
    depth_byelev = pd.DataFrame(np.nan, index = dfind, columns = obj.plotorder)
    sweper_byelev = pd.DataFrame(np.nan, index = dfind, columns = obj.plotorder)
    swevol_byelev = pd.DataFrame(np.nan, index = [str(i) for i in obj.edges] + ['total'], columns = obj.plotorder)
    # swevol_total = pd.DataFrame(np.nan, index = ['total'], columns = obj.plotorder)
    sswe_byelev = pd.DataFrame(index = dfind, columns = obj.plotorder)
    sswe_total = pd.DataFrame(np.nan, index = ['basin mean'], columns = obj.plotorder)
    # dswe_byelev = pd.DataFrame(np.nan, index = dfind, columns = obj.plotorder)
    dswe_byelev = pd.DataFrame(np.nan, index = [str(i) for i in obj.edges] + ['total'], columns = obj.plotorder)

    dswe_byelev.index.name = 'Elevation'
    swe_byelev.index.name = 'Elevation'
    swevol_byelev.index.name = 'Elevation'
    sweper_byelev.index.name = 'Elevation'
    depth_byelev.index.name = 'Elevation'

    sum_flag = True

    if obj.basin in ['KINGS', 'SJ', 'KAWEAH','MERCED']:
        spacecmd = r'\resizebox{\textwidth}{!}{'
    else:
        spacecmd = r'{'

    for sub in obj.plotorder:

        swe_byelev.loc[dfind,sub] = r[ (r['date_time']==end_date)
                                    & (r['variable']=='swe_z')
                                    & (r['basin_id'] == Basins.basins[sub]['basin_id'])]['value'].values[:-1].round(decimals = int(obj.dplcs))
        swe_total.loc['basin mean',sub] = r[ (r['basin_id']==Basins.basins[sub]['basin_id'])
                                    & (r['date_time']==end_date)
                                    & (r['elevation']=='total')
                                    & (r['variable']=='swe_z')]['value'].values[0].round(decimals = int(obj.dplcs))

        swevol_byelev.loc[[str(i) for i in obj.edges] + ['total'],sub] = r[ (r['date_time']==end_date)
                                    & (r['variable']=='swe_vol')
                                    & (r['basin_id'] == Basins.basins[sub]['basin_id'])]['value'].values.round(decimals = int(obj.dplcs))


        sum = r[ (r['date_time']==end_date)
                                    & (r['variable']=='swe_vol')
                                    & (r['basin_id'] == Basins.basins[sub]['basin_id'])]['value'].values[:-1].round(decimals = int(obj.dplcs))

        sweper_byelev.loc[dfind,sub] = (r[ (r['date_time']==end_date)
                                    & (r['variable']=='swe_vol')
                                    & (r['basin_id'] == Basins.basins[sub]['basin_id'])]['value'].values[:-1]/np.nansum(sum)*100).round(decimals = int(obj.dplcs))

        sswe_byelev.loc[dfind,sub] = r[ (r['date_time']==start_date)
                                    & (r['variable']=='swe_z')
                                    & (r['basin_id'] == Basins.basins[sub]['basin_id'])]['value'].values[:-1].round(decimals = int(obj.dplcs))
        sswe_total.loc['basin mean',sub] = r[ (r['basin_id']==Basins.basins[sub]['basin_id'])
                                    & (r['date_time']==start_date)
                                    & (r['elevation']=='total')
                                    & (r['variable']=='swe_z')]['value'].values[0].round(decimals = int(obj.dplcs))

        depth_byelev.loc[dfind,sub] = r[ (r['date_time']==start_date)
                                    & (r['variable']=='depth')
                                    & (r['basin_id'] == Basins.basins[sub]['basin_id'])]['value'].values[:-1].round(decimals = int(obj.dplcs))

        # If the basin total percents will be 0/nan, omit that table
        if ((sub == obj.plotorder[0]) and
            ((np.nansum(sum) < 0.001) or (np.nansum(sum) == np.nan))):
            sum_flag = False

    swe_byelev = swe_byelev.append(swe_total)
    # swevol_byelev = swevol_byelev.append(swevol_total)
    sswe_byelev = sswe_byelev.append(sswe_total)
    dswe_byelev = swe_byelev - sswe_byelev

    colstr = 'l' + 'r'*len(obj.plotorder)

    variables['SWE_BYELEV'] = (
                                r'  \textbf{SWE [%s], %s}\\ \vspace{0.1cm} \\'
                                %(obj.depthlbl,obj.report_date.date().strftime("%Y-%-m-%-d")) +
                                spacecmd + swe_byelev[obj.plotorder].to_latex(na_rep='-', column_format=colstr) +
                                r'} \\ \footnotesize{\textbf{Table 2:} Mean depth of SWE by elevation band.}'
                                )

    variables['SWEVOL_BYELEV'] = (
                                r'  \textbf{SWE volume by elevation [%s], %s}\\ \vspace{0.1cm} \\'
                                %(obj.vollbl,obj.report_date.date().strftime("%Y-%-m-%-d")) +
                                spacecmd + swevol_byelev[obj.plotorder].to_latex(na_rep='-', column_format=colstr) +
                                r'} \\ \footnotesize{\textbf{Table 5:} Volume of SWE by elevation band.}'
                                )
    # variables['SWE_BYELEV'] = variables['SWE_BYELEV'].replace(r'\toprule','')

    variables['DSWE_BYELEV'] = (
                                r'  \textbf{Change in SWE [%s], %s to %s}\\ \vspace{0.1cm} \\'
                                %(obj.depthlbl,obj.report_start.date().strftime("%Y-%-m-%-d"),
                                  obj.report_date.date().strftime("%Y-%-m-%-d")) + spacecmd +
                                  dswe_byelev[obj.plotorder].to_latex(na_rep='-', column_format=colstr) +
                                  r'} \\ \footnotesize{\textbf{Table 3:} Mean change in depth of SWE by elevation band.} ' +
                                  r'\\ \clearpage'
                                )

    # variables['DSWE_BYELEV'] = variables['DSWE_BYELEV'].replace(r'\toprule','')

    if sum_flag is not True:
        variables['SWEPER_BYELEV'] = ' '

    else:
        variables['SWEPER_BYELEV'] = (
                                    r'  \textbf{SWE volume, percent of basin total by elevation band, %s}\\ \vspace{0.1cm} \\'
                                    %(obj.report_date.date().strftime("%Y-%-m-%-d"))+
                                    spacecmd + sweper_byelev[obj.plotorder].round(1).to_latex(na_rep='-', column_format=colstr) +
                                    r'}  \\ \footnotesize{\textbf{Table 4:} Percent of SWE volume by elevation band (totals may not add to 100 due to rounding).} '
                                    )
        variables['SWEPER_BYELEV'] = variables['SWEPER_BYELEV'].replace('inf','-')
        # variables['SWEPER_BYELEV'] = variables['SWEPER_BYELEV'].replace(r'\toprule','')

    # variables['DEPTH_BYELEV'] = (
    #                             r'  \textbf{Snow depth [%s], %s}\\ \vspace{0.1cm} \\'
    #                             %(obj.depthlbl,
    #                               obj.report_date.date().strftime("%Y-%-m-%-d")) + spacecmd +
    #                               depth_byelev[obj.plotorder].to_latex(na_rep='-', column_format=colstr) +
    #                               r'} \\ \footnotesize{\textbf{Table 5:} Mean snow depth by elevation band.} ' +
    #                               r'\\ \clearpage'
    #                             )

    variables['TOT_LBL'] = obj.plotorder[0]

    for n in range(1,len(obj.plotorder)):
        s = 'SUB' + str(n) + '_LBL'
        variables[s] = obj.plotorder[n]

    # Convert floats to strings
    for name in variables:
        if isinstance(variables[name], float):
            if obj.dplcs == 0:
                tmp = str(int(variables[name]))
            else:
                tmp = str(round(variables[name],obj.dplcs))
            variables[name] = tmp

    # Summary sections and fig template have variable strings
    # (e.g. CHANGES_FIG) that need to be replaced
    section_dict = {'SUMMARY':obj.summary_file,
                    'CHANGES_FIG_TPL':obj.figs_tpl_path + 'changes_fig_tpl.txt',
                    'SWI_FIG_TPL':obj.figs_tpl_path + 'swi_fig_tpl.txt',
                    'MEAN_FIG_TPL':obj.figs_tpl_path + 'mean_fig_tpl.txt',
                    'TOTALS_FIG_TPL':obj.figs_tpl_path + 'totals_fig_tpl.txt',
                    'MULTITOTSWE_FIG_TPL':obj.figs_tpl_path + 'multitotswe_fig_tpl.txt',
                    # 'MULTITOTSWI_FIG_TPL':obj.figs_tpl_path + 'multitotswi_fig_tpl.txt',
                    'VALID_FIG_TPL':obj.figs_tpl_path + 'valid_fig_tpl.txt',
                    'FLTCHANGES_FIG_TPL':obj.figs_tpl_path + 'flt_fig_tpl.txt',
                    'PDEP_FIG_TPL':obj.figs_tpl_path + 'pdep_fig_tpl.txt',
                    'COLD_FIG_TPL':obj.figs_tpl_path + 'cold_fig_tpl.txt',
                    'SWE_FIG_TPL':obj.figs_tpl_path + 'swe_fig_tpl.txt',
                    'SUBBASINS_FIG_TPL':obj.figs_tpl_path + 'subbasins_fig_tpl.txt',
                    'INFLOW_FIG_TPL':obj.figs_tpl_path + 'inflow_fig_tpl.txt'
                    # 'RESULTS_FIG_TPL':obj.figs_tpl_path + 'results_fig_tpl.txt',
                    # 'ELEV_FIG_TPL':obj.figs_tpl_path + 'elev_fig_tpl.txt',
                    }

    # Define and load summary tables depending on number of subbasins
    section_dict['PRECIP_SUMMARY_TPL'] = (obj.figs_tpl_path
                                          + 'precip_summary_%ssub.txt'%str(len(obj.plotorder)) )
    section_dict['SWE_SUMMARY_TPL'] = (obj.figs_tpl_path
                                          + 'swe_summary_%ssub.txt'%str(len(obj.plotorder)) )

    # Remove if no flight options
    if obj.flt_flag is False:
        del section_dict['FLTCHANGES_FIG_TPL']

    # Remove if not available
    if obj.plot_flag is True:
        del section_dict['PDEP_FIG_TPL']

    for rep in section_dict.keys():
        fid = open(section_dict[rep],'r')
        var = fid.read()
        fid.close()

        for name in sorted(variables):
            var = var.replace(name,variables[name])
        variables[rep] = var

    # Remove if not BRB
    if (obj.basin not in ['BRB', 'TUOL']) or (obj.basin_total_flag is False):
        variables['MULTITOTSWE_FIG'] = ' '
        variables['MULTITOTSWE_FIG_TPL'] = ' '

    # If figs are listed in exclude, replace with empty string in latex file
    if obj.exclude_figs != None:
        for name in obj.exclude_figs:
            variables[name + '_FIG'] = ' '
            variables[name + '_TPL'] = ' '
            variables[name + '_FIG_TPL'] = ' '
            # variables[name] = ' '

    # Make the report
    env = make_env(loader = FileSystemLoader(obj.templ_path))
    tpl = env.get_template(obj.tex_file)
    # print(tpl.render(variables))
    pdf = build_pdf(tpl.render(variables))
    # To see what's in latex  >>> print(tpl.render(variables))

    # Save in reports and with figs
    rpath_1 = os.path.join(obj.rep_path, '' + obj.report_name)
    rpath_2 = os.path.join(obj.figs_path, '' + obj.report_name)
    obj._logger.info('Saving report to {} and \n{}'.format(rpath_1,rpath_2))

    if not os.path.isdir(os.path.join(obj.rep_path,'')):
        os.makedirs(os.path.join(obj.rep_path,''))

    pdf.save_to(rpath_1)
    pdf.save_to(rpath_2)
