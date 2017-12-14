
def report(object):
    '''
    Need to add flexibility/checking with what band was specified and actually using SWE/SWI
    -forecast date
    
    '''
    
    # SWI values
    total_swi       = object.accum_byelev[object.total_lbl].sum()
    sub1_swi        = object.accum_byelev[object.sub1_lbl].sum()
    sub2_swi        = object.accum_byelev[object.sub2_lbl].sum()
    sub3_swi        = object.accum_byelev[object.sub3_lbl].sum()
    
    # SWE values
    total_swe       = object.state_byelev[object.total_lbl].sum()
    sub1_swe        = object.state_byelev[object.sub1_lbl].sum()
    sub2_swe        = object.state_byelev[object.sub2_lbl].sum()
    sub3_swe        = object.state_byelev[object.sub3_lbl].sum() 
    
    # SWE available for melt
    totalav_swe     = object.melt[object.total_lbl].sum()
    sub1av_swe      = object.melt[object.sub1_lbl].sum()
    sub2av_swe      = object.melt[object.sub2_lbl].sum()
    sub3av_swe      = object.melt[object.sub3_lbl].sum()
    
    # Change in SWE
    total_swe_del   = object.delta_state_byelev[object.total_lbl].sum()
    sub1_swe_del    = object.delta_state_byelev[object.sub1_lbl].sum()
    sub2_swe_del    = object.delta_state_byelev[object.sub2_lbl].sum()
    sub3_swe_del    = object.delta_state_byelev[object.sub3_lbl].sum()


    # Assign variables and put in dict
    start_date  = object.dateFrom.date().strftime("%B %-d")
    end_date    = object.dateTo.date().strftime("%B %-d")  
    
    report_title = object.rep_title
    fore_date   = ''
    swe_in      = total_swe
    swi_in      = total_swi       
    
    fig_path    = object.figs_path
    swi_fig     = 'swi%s.png'%(object.name_append)
    results_fig = 'results%s.png'%(object.name_append)
    changes_fig = 'swe_change%s.png'%(object.name_append)
    # elev_fig    = 'swe_elev%s.png'%(object.name_append)
    
    # Upper case variables are used in the LaTex file, lower case versions are assigned here
    variables = {
                    'REPORT_TITLE':report_title,
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
                    'TOTAL_SWE_DEL':total_swe_del,'SUB1_SWE_DEL':sub1_swe_del,'SUB2_SWE_DEL':sub2_swe_del,'SUB3_SWE_DEL':sub3_swe_del
                    
                    }
    
    #########################################################
    #         SUMMARY
    #########################################################
    
    # Replace variables (like START_DATE) that exist in the summary paragraph
    object.env_path = '/home/markrobertson/mrworkspace/code/SNOWAV/report/brb_template/section_text/'
    
    fid         = open('%ssummary.txt'%(object.env_path),'r')
    summary     = fid.read()
    fid.close()
    
    for name in variables:
        print(name,variables[name])
        summary = summary.replace(name,variables[name])
                 
    fid         = open('%sresults_summary.txt'%(self.env_path),'r')
    results_summary     = fid.read()
    fid.close()
    
    for name in variables:
        results_summary = sr.str_replacer(results_summary,name,variables[name])

           
    # Add the section text variables to what we'll pass to the document
    variables['SUMMARY']            = summary 
    variables['RESULTS_SUMMARY']    = results_summary 
           
    # Make the report
    # env         = make_env(loader=FileSystemLoader('/home/markrobertson/mrworkspace/reports/Tuolumne/tuol_report_template/'))
    # tpl         = env.get_template('tuol_report.tex')
    env         = make_env(loader=FileSystemLoader(self.templ_path))
    tpl         = env.get_template(self.tex_file)        
    pdf         = build_pdf(tpl.render(variables))
    # If you want to see what is being sent to the LaTex environment...
    # print(tpl.render(variables))
    pdf.save_to(self.rep_path) # /mnt/volumes/wkspace/reports/wy17/brb/BRB_Snowpack_Summary_test.pdf

