from jinja2 import FileSystemLoader
from latex.jinja2 import make_env
from latex import build_pdf
import pandas as pd
import numpy as np
from datetime import datetime
import copy


def report(obj,*args):
    '''
    Need to add flexibility/checking with what band was specified and actually using SWE/SWI
    -forecast date
    -currently the string replacer only handles str and float
    -simple check for if figures even exist
    -make report environment paths relative
    -consider printing out summary info for sanity check...
    
    '''
    
    if len(args) != 0:
        # parts = obj.report_name
        obj.report_name = obj.report_name.split('.')[0] + args[0] + '.pdf'
    
    total_swi       = obj.accum_byelev[obj.total_lbl].sum()
    total_swe       = obj.state_byelev[obj.total_lbl].sum()
    totalav_swe     = obj.melt[obj.total_lbl].sum()
    total_swe_del   = obj.delta_state_byelev[obj.total_lbl].sum()
    report_time     = datetime.now().strftime("%Y-%-m-%-d %H:%M")
    
    # HACK, somewhere the rounding is different when we start from Oct 1...
    if obj.dateFrom == datetime(2017,10,1,23,0):
        total_swe_del = copy.copy(total_swe)
    
    # Total hack for different number of subbasins, needs improvement...
    if hasattr(obj, 'sub1_lbl') and hasattr(obj, 'sub2_lbl') and hasattr(obj, 'sub3_lbl'):
        
        # values
        sub1_swi        = obj.accum_byelev[obj.sub1_lbl].sum()
        sub2_swi        = obj.accum_byelev[obj.sub2_lbl].sum()
        sub3_swi        = obj.accum_byelev[obj.sub3_lbl].sum()
        sub1_swe        = obj.state_byelev[obj.sub1_lbl].sum()
        sub2_swe        = obj.state_byelev[obj.sub2_lbl].sum()
        sub3_swe        = obj.state_byelev[obj.sub3_lbl].sum()
        sub1av_swe      = obj.melt[obj.sub1_lbl].sum()
        sub2av_swe      = obj.melt[obj.sub2_lbl].sum()
        sub3av_swe      = obj.melt[obj.sub3_lbl].sum()
        sub1_swe_del    = obj.delta_state_byelev[obj.sub1_lbl].sum()
        sub2_swe_del    = obj.delta_state_byelev[obj.sub2_lbl].sum()
        sub3_swe_del    = obj.delta_state_byelev[obj.sub3_lbl].sum()
        
        
        # Calculate mean basin depths in mm
        # if want mm, np.multiply with (1/obj.depth_factor)
        total_pm        = np.nansum(np.multiply(obj.state*obj.masks[obj.total_lbl]['mask'],1))/obj.masks[obj.total_lbl]['mask'].sum()
        sub1_pm         = np.nansum(np.multiply(obj.state*obj.masks[obj.sub1_lbl]['mask'],1))/obj.masks[obj.sub1_lbl]['mask'].sum()
        sub2_pm         = np.nansum(np.multiply(obj.state*obj.masks[obj.sub2_lbl]['mask'],1))/obj.masks[obj.sub2_lbl]['mask'].sum()
        sub3_pm         = np.nansum(np.multiply(obj.state*obj.masks[obj.sub3_lbl]['mask'],1))/obj.masks[obj.sub3_lbl]['mask'].sum()
        
        # More SWI stuff
        total_mel       = obj.snowmelt_byelev[obj.total_lbl].sum()
        sub1_mel        = obj.snowmelt_byelev[obj.sub1_lbl].sum()
        sub2_mel        = obj.snowmelt_byelev[obj.sub2_lbl].sum()
        sub3_mel        = obj.snowmelt_byelev[obj.sub3_lbl].sum()
        
        # total_rai       = total_swi - total_mel
        # sub1_rai        = sub1_swi - sub1_mel
        # sub2_rai        = sub2_swi - sub2_mel
        # sub3_rai        = sub3_swi - sub3_mel
        
        # total_rai       = obj.rain_bg_byelev[obj.total_lbl].sum()
        # sub1_rai        = obj.rain_bg_byelev[obj.sub1_lbl].sum()
        # sub2_rai        = obj.rain_bg_byelev[obj.sub2_lbl].sum()
        # sub3_rai        = obj.rain_bg_byelev[obj.sub3_lbl].sum()
        
        total_rai       = obj.rain_bg_byelev[obj.total_lbl].sum()
        sub1_rai        = obj.rain_bg_byelev[obj.sub1_lbl].sum()
        sub2_rai        = obj.rain_bg_byelev[obj.sub2_lbl].sum()
        sub3_rai        = obj.rain_bg_byelev[obj.sub3_lbl].sum()        
        
        # Calculate mean rain in mm
        totalpre_pm    = np.nansum(np.multiply(obj.precip*obj.masks[obj.total_lbl]['mask'],1))/obj.masks[obj.total_lbl]['mask'].sum()
        sub1pre_pm     = np.nansum(np.multiply(obj.precip*obj.masks[obj.sub1_lbl]['mask'],1))/obj.masks[obj.sub1_lbl]['mask'].sum()
        sub2pre_pm     = np.nansum(np.multiply(obj.precip*obj.masks[obj.sub2_lbl]['mask'],1))/obj.masks[obj.sub2_lbl]['mask'].sum()
        sub3pre_pm     = np.nansum(np.multiply(obj.precip*obj.masks[obj.sub3_lbl]['mask'],1))/obj.masks[obj.sub3_lbl]['mask'].sum()
            
        
        # DEBUG for rain calculation
        for s in obj.plotorder:
            sub = obj.accum_byelev[s].sum() - obj.snowmelt_byelev[s].sum()
            cal = obj.rain_bg_byelev[s].sum()
            print(sub,cal)
        
        total_rat       = str(int((total_rai/total_swi)*100))
        sub1_rat        = str(int((sub1_rai/sub1_swi)*100))
        sub2_rat        = str(int((sub2_rai/sub2_swi)*100))
        sub3_rat        = str(int((sub3_rai/sub3_swi)*100))
    
    else: 

        sub3_swi        = 0.0
        sub3_swe        = 0.0
        sub3av_swe      = 0.0     
        sub3_swe_del    = 0.0 
    

    # Assign some more variables
    start_date      = obj.dateFrom.date().strftime("%B %-d")
    end_date        = obj.dateTo.date().strftime("%B %-d")  
    
    report_title    = obj.rep_title
    fore_date       = ' '
    swe_in          = total_swe
    swi_in          = total_swi 
    vollbl          = obj.vollbl
    deplbl          = obj.depthlbl        
    
    fig_path        = obj.figs_path
    swi_fig         = 'swi%s.png'%(obj.name_append)
    results_fig     = 'results%s.png'%(obj.name_append)
    changes_fig     = 'swe_change%s.png'%(obj.name_append)
    elev_fig        = 'swe_elev%s.png'%(obj.name_append)
    totals_fig      = 'basin_total%s.png'%(obj.name_append)
    totalsmy_fig    = 'basin_total_multiyr%s.png'%(obj.name_append)
    valid_fig       = 'validation%s.png'%(obj.name_append)
    hyp_fig         = 'hypsometry%s.png'%(obj.name_append)
    
    if obj.tex_file == 'tuol_report_flt.tex':
        changes_flt_fig = 'swe_change_flt%s.png'%(obj.name_append)
        pre_swe         = obj.pre_swe
        pre_pm          = obj.pre_pm
        diff_pm         = total_pm - pre_pm
        
    if obj.units == 'SI':
        unitlbl = '$km^3$'
    else:
        unitlbl = vollbl
    # Upper case variables are used in the LaTex file, lower case versions are assigned here
    variables = {
                    'REPORT_TITLE':report_title,
                    'REPORT_TIME':report_time,
                    'WATERYEAR':str(obj.wy),
                    'UNITS':unitlbl,
                    'VOLLBL':vollbl,
                    'DEPLBL':deplbl,
                    'START_DATE':start_date,
                    'END_DATE':end_date,
                    'FORE_DATE':fore_date,
                    'SWE_IN':swe_in,
                    'SWI_IN':swi_in,
                    'FIG_PATH':fig_path,
                    'SWI_FIG':swi_fig,
                    'RESULTS_FIG':results_fig,
                    'CHANGES_FIG':changes_fig,
                    'ELEV_FIG':elev_fig,
                    'TOTALS_FIG':totals_fig,
                    'TOTALSMY_FIG':totalsmy_fig,
                    'HYP_FIG':hyp_fig,
                                    
                    'TOTAL_SWI':total_swi,'SUB1_SWI':sub1_swi,'SUB2_SWI':sub2_swi,'SUB3_SWI':sub3_swi,
                    'TOTAL_MEL':total_mel,'SUB1_MEL':sub1_mel,'SUB2_MEL':sub2_mel,'SUB3_MEL':sub3_mel,
                    'TOTAL_RAI':total_rai,'SUB1_RAI':sub1_rai,'SUB2_RAI':sub2_rai,'SUB3_RAI':sub3_rai,
                    'TOTAL_RAT':total_rat,'SUB1_RAT':sub1_rat,'SUB2_RAT':sub2_rat,'SUB3_RAT':sub3_rat,
                    
                    'TOTAL_SWE':total_swe,'SUB1_SWE':sub1_swe,'SUB2_SWE':sub2_swe,'SUB3_SWE':sub3_swe,
                    'TOTAL_SWE_AV':totalav_swe,'SUB1_SWE_AV':sub1av_swe,'SUB2_SWE_AV':sub2av_swe,'SUB3_SWE_AV':sub3av_swe,
                    'TOTAL_PM':total_pm,'SUB1_PM':sub1_pm,'SUB2_PM':sub2_pm,'SUB3_PM':sub3_pm,
                    'TOTALPRE_PM':totalpre_pm,'SUB1PRE_PM':sub1pre_pm,'SUB2PRE_PM':sub2pre_pm,'SUB3PRE_PM':sub3pre_pm,
                    'TOTAL_SWEDEL':total_swe_del,'SUB1_SWEDEL':sub1_swe_del,'SUB2_SWEDEL':sub2_swe_del,'SUB3_SWEDEL':sub3_swe_del
                    
                    }
    
    if obj.basin == 'BRB':
        variables['VALID_FIG'] = valid_fig
        
    if obj.basin == 'TUOL' and obj.tex_file == 'tuol_report_flt.tex':
        variables['CHANGES_FLT_FIG'] = changes_flt_fig
        variables['PRE_SWE']        = pre_swe
        variables['PRE_PM']         = pre_pm
        variables['DIFF_PM']        = diff_pm
        variables['FLT_SWEDEL']     = total_swe - pre_swe
    
    # Convert floats to strings    
    for name in variables:
        try:
            if isinstance(variables[name], float):
                if obj.dplcs == 0:
                    tmp     = str(int(variables[name]))
                else:
                    tmp     = str(round(variables[name],obj.dplcs))
                variables[name] = tmp
        except:
            print('Failed converting variables to strings for report...')
    
    # Load in summary.txt files, and replace variables with values
    if obj.basin == 'BRB':
        summary         = '%sbrb_summary.txt'%(obj.env_path)
        results_summary = '%sbrb_results_summary.txt'%(obj.env_path)
        
    if obj.basin == 'TUOL':
        if obj.tex_file == 'tuol_report_flt.tex':
            summary         = '%stuol_summary_flt.txt'%(obj.env_path)
        else:
            summary         = '%stuol_summary.txt'%(obj.env_path)
        results_summary = '%stuol_results_summary.txt'%(obj.env_path)

    if obj.basin == 'SJ':
        summary         = '%ssj_summary.txt'%(obj.env_path)
        results_summary = '%ssj_results_summary.txt'%(obj.env_path)
    
    # Read in both section summaries and then replace variables
    fid                 = open(summary,'r')
    fid1                = open(results_summary,'r')
    summary             = fid.read()
    results_summary     = fid1.read()
    fid.close()
    fid1.close()  
    
    for name in variables:
        summary         = summary.replace(name,variables[name])         
        results_summary = results_summary.replace(name,variables[name]) 
            
    # Add the section text variables to what we'll pass to the document
    variables['SUMMARY']            = summary 
    variables['RESULTS_SUMMARY']    = results_summary 
           
    # Make the report
    env         = make_env(loader = FileSystemLoader(obj.templ_path))
    tpl         = env.get_template(obj.tex_file)        
    pdf         = build_pdf(tpl.render(variables))
    # To see what the latex version >>> print(tpl.render(variables))
    
    print('Saving report to %s%s'%(obj.rep_path,obj.report_name))
    pdf.save_to('%s%s'%(obj.rep_path,obj.report_name)) 

