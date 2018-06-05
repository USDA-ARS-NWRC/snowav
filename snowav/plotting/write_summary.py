
# import pandas as pd

def write_summary(snow,df):

    try:
        if df == 'state':
            snow.state_summary.to_csv('%s%s_summary%s.csv'%(snow.figs_path,df,snow.name_append))
        if df == 'accum':
            snow.accum_summary.to_csv('%s%s_summary%s.csv'%(snow.figs_path,df,snow.name_append))
        if df == 'precip':
            snow.precip_summary.to_csv('%s%s_summary%s.csv'%(snow.figs_path,df,snow.name_append))
        if df not in ['accum','state','precip']:
            print('Unable to write %s, maybe check yo config file'%(df))    
    except:    
        print('Could not write %s dataframe...'%(df))