from jinja2 import FileSystemLoader
from latex.jinja2 import make_env
from latex import build_pdf
import pandas as pd
import numpy as np


def report(obj):
    '''
    Need to add flexibility/checking with what band was specified and actually using SWE/SWI
    -forecast date
    -currently the string replacer only handles str and float
    -simple check for if figures even exist
    -make report environment paths relative
    -consider printing out summary info for sanity check...
    
    '''
    
    # SWI values
    total_swi       = obj.accum_byelev[obj.total_lbl].sum()
    sub1_swi        = obj.accum_byelev[obj.sub1_lbl].sum()
    sub2_swi        = obj.accum_byelev[obj.sub2_lbl].sum()
    sub3_swi        = obj.accum_byelev[obj.sub3_lbl].sum()
    
    # SWE values
    total_swe       = obj.state_byelev[obj.total_lbl].sum()
    sub1_swe        = obj.state_byelev[obj.sub1_lbl].sum()
    sub2_swe        = obj.state_byelev[obj.sub2_lbl].sum()
    sub3_swe        = obj.state_byelev[obj.sub3_lbl].sum()
    
    # SWE available for melt
    totalav_swe     = obj.melt[obj.total_lbl].sum()
    sub1av_swe      = obj.melt[obj.sub1_lbl].sum()
    sub2av_swe      = obj.melt[obj.sub2_lbl].sum()
    sub3av_swe      = obj.melt[obj.sub3_lbl].sum()
    
    # Change in SWE
    total_swe_del   = obj.delta_state_byelev[obj.total_lbl].sum()
    sub1_swe_del    = obj.delta_state_byelev[obj.sub1_lbl].sum()
    sub2_swe_del    = obj.delta_state_byelev[obj.sub2_lbl].sum()
    sub3_swe_del    = obj.delta_state_byelev[obj.sub3_lbl].sum()

    # Assign some more variables
    start_date      = obj.dateFrom.date().strftime("%B %-d")
    end_date        = obj.dateTo.date().strftime("%B %-d")  
    
    report_title    = obj.rep_title
    fore_date       = ' '
    swe_in          = total_swe
    swi_in          = total_swi       
    
    fig_path        = obj.figs_path
    swi_fig         = 'swi%s.png'%(obj.name_append)
    results_fig     = 'results%s.png'%(obj.name_append)
    changes_fig     = 'swe_change%s.png'%(obj.name_append)
    # elev_fig    = 'swe_elev%s.png'%(obj.name_append)
    
    # Check that figures actually exist
    # for name in []
    #     os.path.isfile() 
    
    # Upper case variables are used in the LaTex file, lower case versions are assigned here
    variables = {
                    'REPORT_TITLE':report_title,
                    'WATERYEAR':str(obj.wy),
                    'UNITS':obj.reportunits,
                    'START_DATE':start_date,
                    'END_DATE':end_date,
                    'FORE_DATE':fore_date,
                    'SWE_IN':swe_in,
                    'SWI_IN':swi_in,
                    'FIG_PATH':fig_path,
                    'SWI_FIG':swi_fig,
                    'RESULTS_FIG':results_fig,
                    'CHANGES_FIG':changes_fig,
                    #'ELEV_FIG':elev_fig,
                                    
                    'TOTAL_SWI':total_swi,'SUB1_SWI':sub1_swi,'SUB2_SWI':sub2_swi,'SUB3_SWI':sub3_swi,
                    'TOTAL_SWE':total_swe,'SUB1_SWE':sub1_swe,'SUB2_SWE':sub2_swe,'SUB3_SWE':sub3_swe,
                    'TOTAL_SWE_AV':totalav_swe,'SUB1_SWE_AV':sub1av_swe,'SUB2_SWE_AV':sub2av_swe,'SUB3_SWE_AV':sub3av_swe,
                    'TOTAL_SWEDEL':total_swe_del,'SUB1_SWEDEL':sub1_swe_del,'SUB2_SWEDEL':sub2_swe_del,'SUB3_SWEDEL':sub3_swe_del
                    
                    }
    
    # Convert floats to strings    
    for name in variables:
        try:
            if isinstance(variables[name], float):
                tmp     = str(int(variables[name]))
                variables[name] = tmp
        except:
            print('Failed converting variables to strings for report...')
    
    # Load in summary.txt files, and replace variables with values
    if obj.basin == 'BRB':
        summary         = '%sbrb_summary.txt'%(obj.env_path)
        results_summary = '%sbrb_results_summary.txt'%(obj.env_path)
        
    if obj.basin == 'TUOL':
        summary         = '%stuol_summary.txt'%(obj.env_path)
        results_summary = '%stuol_results_summary.txt'%(obj.env_path)
    
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
    # print(tpl.render(variables))
    
    print('Saving report to %s%s'%(obj.rep_path,obj.report_name))
    pdf.save_to('%s%s'%(obj.rep_path,obj.report_name)) 

