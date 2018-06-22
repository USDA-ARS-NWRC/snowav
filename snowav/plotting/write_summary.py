import os

def write_summary(snow,df):

    fp = os.path.join(snow.figs_path, df+'_summary'+snow.name_append+'.csv')
    if df == 'state':
        snow.state_summary.to_csv(fp)
    if df == 'accum':
        snow.accum_summary.to_csv(fp)
    if df == 'precip':
        snow.precip_summary.to_csv(fp)
    if df not in ['accum','state','precip']:
        snow._logger.info('Unable to write %s, maybe check yo config file'%(df))
