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
       
    variables['TOTAL_SWI'] = obj.accum_byelev[obj.total_lbl].sum()
    variables['TOTAL_SWE'] = obj.state_byelev[obj.total_lbl].sum()
    variables['TOTAL_SWE_AV'] = obj.melt[obj.total_lbl].sum()
    variables['TOTAL_SWEDEL'] = obj.delta_state_byelev[obj.total_lbl].sum()
    variables['TOTAL_PM'] = (np.nansum(
                            np.multiply(obj.state 
                            * obj.masks[obj.total_lbl]['mask'],1))
                            / obj.masks[obj.total_lbl]['mask'].sum())
    variables['TOTALPRE_PM'] = (np.nansum(
                                np.multiply(obj.precip 
                                * obj.masks[obj.total_lbl]['mask'],1))
                                / obj.masks[obj.total_lbl]['mask'].sum() )
    total_rai = obj.rain_bg_byelev[obj.total_lbl].sum()
    variables['TOTAL_RAT'] = str(int((total_rai 
                                      /obj.accum_byelev[obj.total_lbl].sum())*100))
  
    report_time = datetime.now().strftime("%Y-%-m-%-d %H:%M")   
    if hasattr(obj,'orig_date'):
        report_time = obj.orig_date + ', revised ' + report_time
   
    # HACK, somewhere the rounding is different when we start from Oct 1...
    # if obj.dateFrom == datetime(2017,10,1,23,0):
    #     total_swe_del = copy.copy(total_swe)
    
    numsubs = range(1,len(obj.plotorder))

    for n,sub in zip(numsubs,obj.suborder):
        SWIIND = 'SUB' + str(n) + '_SWI'
        SWEIND = 'SUB' + str(n) + '_SWE'
        AVSWEIND = 'SUB' + str(n) + '_SWE_AV'
        SWEDEL = 'SUB' + str(n) + '_SWEDEL'
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
    
    report_title = obj.rep_title
    fore_date = ' '
    swe_in = obj.state_byelev[obj.total_lbl].sum()
    swi_in = obj.accum_byelev[obj.total_lbl].sum()
    vollbl = obj.vollbl
    deplbl = obj.depthlbl        
    
    fig_path = obj.figs_path
    swi_fig = 'swi%s.png'%(obj.name_append)
    results_fig = 'results%s.png'%(obj.name_append)
    changes_fig = 'swe_change%s.png'%(obj.name_append)
    elev_fig = 'swe_elev%s.png'%(obj.name_append)
    totals_fig = 'basin_total%s.png'%(obj.name_append)
    totalsmy_fig = 'basin_total_multiyr%s.png'%(obj.name_append)
    valid_fig = 'validation%s.png'%(obj.name_append)
    hyp_fig  = 'hypsometry%s.png'%(obj.name_append)
    mean_fig = 'mean_swe_depth%s.png'%(obj.name_append)
    detail_fig = 'mean_detail%s.png'%(obj.name_append)
    density_fig = 'density%s.png'%(obj.name_append)
    density_sub_fig = 'density_sub%s.png'%(obj.name_append)
    density_swe_fig = 'density_swe%s.png'%(obj.name_append)
    
    if obj.tex_file == 'tuol_report_flt.tex':
        changes_flt_fig = 'swe_change_flt%s.png'%(obj.name_append)
        pre_swe = obj.pre_swe
        pre_pm = obj.pre_pm
        diff_pm = variables['TOTAL_PM'] - pre_pm
        
    if obj.units == 'SI':
        unitlbl = '$km^3$'
    else:
        unitlbl = vollbl
        
    # Upper case variables are used in the LaTex file, 
    # lower case versions are assigned here
    variables['REPORT_TITLE'] = report_title
    variables['REPORT_TIME'] = report_time
    variables['WATERYEAR'] = str(obj.wy)
    variables['UNITS'] = unitlbl
    variables['VOLLBL'] = vollbl
    variables['DEPLBL'] = deplbl
    variables['START_DATE'] = start_date
    variables['END_DATE'] = end_date
    variables['FORE_DATE'] = fore_date
    variables['SWE_IN'] = swe_in
    variables['SWI_IN'] = swi_in
    variables['FIG_PATH'] = fig_path
    variables['SWI_FIG'] = swi_fig
    variables['RESULTS_FIG'] = results_fig
    variables['CHANGES_FIG'] = changes_fig
    variables['ELEV_FIG'] = elev_fig
    variables['TOTALS_FIG'] = totals_fig
    variables['TOTALSMY_FIG'] = totalsmy_fig
    variables['HYP_FIG'] = hyp_fig         
    variables['MEAN_FIG'] = mean_fig
    variables['DETAIL_FIG'] = detail_fig 
    variables['DENSITY_FIG'] = density_fig    
    variables['DENSITY_SUB_FIG'] = density_sub_fig 
    variables['DENSITY_SWE_FIG'] = density_swe_fig                                  
    variables['VALID_FIG'] = valid_fig
        
    if obj.basin == 'TUOL' and obj.tex_file == 'tuol_report_flt.tex':
        variables['CHANGES_FLT_FIG'] = changes_flt_fig
        variables['PRE_SWE'] = pre_swe
        variables['PRE_PM'] = pre_pm
        variables['DIFF_PM'] = diff_pm
        variables['FLT_SWEDEL'] = variables['TOTAL_SWE'] - pre_swe
    
    # Convert floats to strings    
    for name in variables:
        if isinstance(variables[name], float):
            if obj.dplcs == 0:
                tmp = str(int(variables[name]))
            else:
                tmp = str(round(variables[name],obj.dplcs))
            variables[name] = tmp
  
    # Load in summary.txt files, and replace variables with values
    if obj.basin == 'BRB':
        summary = '%sbrb_summary.txt'%(obj.env_path)
        results_summary = '%sbrb_results_summary.txt'%(obj.env_path)
        
    if obj.basin == 'TUOL':
        if obj.tex_file == 'tuol_report_flt.tex':
            summary = '%stuol_summary_flt.txt'%(obj.env_path)
        else:
            summary = '%stuol_summary.txt'%(obj.env_path)
        results_summary = '%stuol_results_summary.txt'%(obj.env_path)

    if obj.basin == 'SJ':
        
        if obj.tex_file == 'sj_report_hrrr.tex':
            summary = '%ssj_summary_hrrr.txt'%(obj.env_path)
        else:
            summary = '%ssj_summary.txt'%(obj.env_path)
        
        results_summary = '%ssj_results_summary.txt'%(obj.env_path)
        
    if obj.basin == 'LAKES':
        summary = '%slakes_summary.txt'%(obj.env_path)
        results_summary = '%slakes_results_summary.txt'%(obj.env_path)   
    if obj.basin == 'RCEW':
        summary = '%srcew_summary.txt'%(obj.env_path)
        results_summary = '%srcew_results_summary.txt'%(obj.env_path)               
    
    # Read in both section summaries and then replace variables
    fid = open(summary,'r')
    fid1 = open(results_summary,'r')
    summary = fid.read()
    results_summary = fid1.read()
    fid.close()
    fid1.close()  
    
    for name in variables:
        summary = summary.replace(name,variables[name])         
        results_summary = results_summary.replace(name,variables[name]) 
            
    # Add the section text variables to what we'll pass to the document
    variables['SUMMARY'] = summary 
    variables['RESULTS_SUMMARY'] = results_summary 
           
    # Make the report
    env = make_env(loader = FileSystemLoader(obj.templ_path))
    tpl = env.get_template(obj.tex_file)        
    pdf = build_pdf(tpl.render(variables))
    # To see what the latex version >>> print(tpl.render(variables))
    
    # Save in reports and with figs
    print('Saving report to %s%s'%(obj.rep_path,obj.report_name))
    pdf.save_to('%s%s'%(obj.rep_path,obj.report_name)) 
    pdf.save_to('%s%s'%(obj.figs_path,obj.report_name)) 

