
import pandas as pd

def write_summary(snow,df):

    dataframes  = ['accum_summary','state_summary']

    if df not in dataframes:
        print('df needs to be one of: %s'%(dataframes))

    print('Writing summary file to %s%s%s.csv'%(snow.figs_path,df,snow.name_append))

    if df == 'state_summary':
        snow.state_summary.to_csv('%s%s%s.csv'%(snow.figs_path,df,snow.name_append))
    if df == 'accum_summary':
        snow.accum_summary.to_csv('%s%s%s.csv'%(snow.figs_path,df,snow.name_append))
