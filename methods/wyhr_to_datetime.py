
import datetime

def wyhr_to_datetime(wy,wyhr):
    # wydt = wyhr_to_datetime(wy,wyhr)
    wyst        = datetime.datetime(wy-1,10,1,0,0,0)
    hrs         = datetime.timedelta(0,wyhr*60*60)
    wydatetime  = wyst + hrs
    
    return wydatetime