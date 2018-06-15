from jinja2 import FileSystemLoader
from latex.jinja2 import make_env
from latex import build_pdf
import pandas as pd
import numpy as np
from datetime import datetime
import copy


def report(obj):
    '''

    '''

    # Initialize all the variables to pass to latex file
    variables = {}
    variables['TOTAL_SWI'] = obj.accum_byelev[obj.plotorder[0]].sum()
    variables['TOTAL_SWE'] = obj.state_byelev[obj.plotorder[0]].sum()
    variables['TOTAL_SAV'] = obj.melt[obj.plotorder[0]].sum()
    variables['TOTAL_SDEL'] = obj.delta_state_byelev[obj.plotorder[0]].sum()
    variables['TOTAL_PM'] = (np.nansum(
                            np.multiply(obj.state
                            * obj.masks[obj.plotorder[0]]['mask'],1))
                            / obj.masks[obj.plotorder[0]]['mask'].sum())
    variables['TOTALPRE_PM'] = (np.nansum(
                                np.multiply(obj.precip
                                * obj.masks[obj.plotorder[0]]['mask'],1))
                                / obj.masks[obj.plotorder[0]]['mask'].sum() )
    total_rai = obj.rain_bg_byelev[obj.plotorder[0]].sum()
    variables['TOTAL_RAT'] = str(int((total_rai
                                      /obj.accum_byelev[obj.plotorder[0]].sum())*100))

    report_time = datetime.now().strftime("%Y-%-m-%-d %H:%M")
    if hasattr(obj,'orig_date'):
        report_time = obj.orig_date + ', revised ' + report_time

    # HACK, somewhere the rounding is different when we start from Oct 1...
    # if obj.dateFrom == datetime(2017,10,1,23,0):
    #     total_swe_del = copy.copy(total_swe)

    numsubs = range(1,len(obj.plotorder))

    for n,sub in zip(numsubs,obj.plotorder[1:]):
        SWIIND = 'SUB' + str(n) + '_SWI'
        SWEIND = 'SUB' + str(n) + '_SWE'
        AVSWEIND = 'SUB' + str(n) + '_SAV'
        SWEDEL = 'SUB' + str(n) + '_SDEL'
        PM = 'SUB' + str(n) + '_PM'
        PREPM = 'SUB' + str(n) + 'PRE_PM'
        RAIN = 'SUB' + str(n) + 'RAI'
        RATIO = 'SUB' + str(n) + '_RAT'

        swival = obj.accum_byelev[sub].sum()
        sweval = obj.state_byelev[sub].sum()
        avsweval = obj.melt[sub].sum()
        swedelval = obj.delta_state_byelev[sub].sum()
        pmval = (np.nansum(np.multiply(obj.state
                                      * obj.masks[sub]['mask'],1))
                                      / obj.masks[sub]['mask'].sum())
        prepmval = (np.nansum(np.multiply(obj.precip*obj.masks[sub]['mask'],1))
                        /obj.masks[sub]['mask'].sum())
        rainval = obj.rain_bg_byelev[sub].sum()
        ratval = str(int((rainval/swival)*100))

        variables[SWIIND] = swival
        variables[SWEIND] = sweval
        variables[AVSWEIND] = avsweval
        variables[SWEDEL] = swedelval
        variables[PM] = pmval
        variables[PREPM] = prepmval
        variables[RAIN] = rainval
        variables[RATIO]  = ratval

    # Assign some more variables
    start_date = obj.dateFrom.date().strftime("%B %-d")
    end_date = obj.dateTo.date().strftime("%B %-d")
    fore_date = ' '

    # Upper case variables are used in the LaTex file,
    # lower case versions are assigned here
    variables['REPORT_TITLE'] = obj.rep_title
    variables['REPORT_TIME'] = report_time
    variables['WATERYEAR'] = str(obj.wy)
    variables['UNITS'] = obj.vollbl
    variables['VOLLBL'] = obj.vollbl
    variables['DEPLBL'] = obj.depthlbl
    variables['START_DATE'] = start_date
    variables['END_DATE'] = end_date
    variables['FORE_DATE'] = fore_date
    variables['SWE_IN'] = obj.state_byelev[obj.plotorder[0]].sum()
    variables['SWI_IN'] = obj.accum_byelev[obj.plotorder[0]].sum()
    variables['FIG_PATH'] = obj.figs_path
    variables['SWI_FIG'] = 'swi%s.png'%(obj.name_append)
    variables['RESULTS_FIG'] = 'results%s.png'%(obj.name_append)
    variables['CHANGES_FIG'] = 'swe_change%s.png'%(obj.name_append)
    variables['FLT_CHANGES_FIG'] = 'flt_swe_change%s.png'%(obj.name_append)
    variables['CHANGES_DEP_FIG'] = 'swe_change_depth%s.png'%(obj.name_append)
    variables['ELEV_FIG'] = 'swe_elev%s.png'%(obj.name_append)
    variables['TOTALS_FIG'] = 'basin_total%s.png'%(obj.name_append)
    variables['TOTALSMY_FIG'] = 'basin_total_multiyr%s.png'%(obj.name_append)
    variables['HYP_FIG'] = 'hypsometry%s.png'%(obj.name_append)
    variables['MEAN_FIG'] = 'mean_swe_depth%s.png'%(obj.name_append)
    variables['DETAIL_FIG'] = 'mean_detail%s.png'%(obj.name_append)
    variables['DENSITY_FIG'] = 'density%s.png'%(obj.name_append)
    variables['DENSITY_SUB_FIG'] = 'density_sub%s.png'%(obj.name_append)
    variables['DENSITY_SWE_FIG'] = 'density_swe%s.png'%(obj.name_append)
    variables['VALID_FIG'] = 'validation%s.png'%(obj.name_append)
    variables['ARSLOGO'] = obj.figs_tpl_path + 'ars_logo.png'
    variables['AWSMLOGO'] = obj.figs_tpl_path + 'logo.png'

    # Put the by-elevation tables together
    # Eventually put this in framework? But _o is beating it to hell right now
    obj.state_mswe_byelev.index.name = 'Elevation'
    obj.state_mswe_byelev = obj.state_mswe_byelev.astype('float64')
    obj.state_mswe_byelev = obj.state_mswe_byelev.fillna(0).round(1)
    variables['SWE_BYELEV'] = (
                                r'\textbf{Mean SWE [%s], %s}\\ \vspace{0.1cm} \\'
                                %(obj.depthlbl,obj.dateTo.date().strftime("%Y-%-m-%-d"))
                                + obj.state_mswe_byelev[obj.plotorder].to_latex()
                                )

    obj.delta_swe_byelev.index.name = 'Elevation'
    obj.delta_swe_byelev = obj.delta_swe_byelev.astype('float64')
    obj.delta_swe_byelev = obj.delta_swe_byelev.fillna(0).round(1)
    variables['DSWE_BYELEV'] = (
                                r'\textbf{Change in SWE [%s], %s to %s}\\ \vspace{0.1cm} \\'
                                %(obj.depthlbl,obj.dateFrom.date().strftime("%Y-%-m-%-d"),
                                  obj.dateTo.date().strftime("%Y-%-m-%-d"))
                                + obj.delta_swe_byelev[obj.plotorder].to_latex()
                                )

    variables['TOT_LBL'] = obj.plotorder[0]
    if len(obj.plotorder) >= 2:
        variables['SUB1_LBL'] = obj.plotorder[1]
    if len(obj.plotorder) >= 3:
        variables['SUB2_LBL'] = obj.plotorder[2]
    if len(obj.plotorder) >= 4:
        variables['SUB3_LBL'] = obj.plotorder[3]
    if len(obj.plotorder) >= 5:
        variables['SUB4_LBL'] = obj.plotorder[4]

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
                    'RESULTS_FIG_TPL':obj.figs_tpl_path + 'results_fig_tpl.txt',
                    'ELEV_FIG_TPL':obj.figs_tpl_path + 'elev_fig_tpl.txt',
                    'MEAN_FIG_TPL':obj.figs_tpl_path + 'mean_fig_tpl.txt',
                    'TOTALS_FIG_TPL':obj.figs_tpl_path + 'totals_fig_tpl.txt',
                    'VALID_FIG_TPL':obj.figs_tpl_path + 'valid_fig_tpl.txt',
                    'FLT_CHANGES_FIG_TPL':obj.figs_tpl_path + 'flt_changes_fig_tpl.txt'
                    }

    # Define and load summary tables depending on number of subbasins
    section_dict['PRECIP_SUMMARY_TPL'] = (obj.figs_tpl_path
                                          + 'precip_summary_%ssub.txt'%str(len(obj.plotorder)) )
    section_dict['SWE_SUMMARY_TPL'] = (obj.figs_tpl_path
                                          + 'swe_summary_%ssub.txt'%str(len(obj.plotorder)) )


    # Remove if no flight options
    if not hasattr(obj,'flt_flag'):
        del section_dict['FLT_CHANGES_FIG_TPL']

    for rep in section_dict.keys():
        fid = open(section_dict[rep],'r')
        var = fid.read()
        fid.close()

        for name in variables:
            var = var.replace(name,variables[name])
        variables[rep] = var

    # If figs are listed in exclude, replace with empty string in latex file
    if obj.exclude_figs != None:
        for name in obj.exclude_figs:
            variables[name + '_FIG_TPL'] = ' '

    # Make the report
    env = make_env(loader = FileSystemLoader(obj.templ_path))
    tpl = env.get_template(obj.tex_file)
    pdf = build_pdf(tpl.render(variables))
    # To see what's in latex  >>> print(tpl.render(variables))

    # Save in reports and with figs
    obj._logger.info('Saving report to %s%s and \n%s%s'%(obj.rep_path,obj.report_name,obj.rep_path,obj.report_name))
    pdf.save_to('%s%s'%(obj.rep_path,obj.report_name))
    pdf.save_to('%s%s'%(obj.figs_path,obj.report_name))
