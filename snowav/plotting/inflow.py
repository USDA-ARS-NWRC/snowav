
import pandas as pd
from matplotlib import pyplot as plt
from snowav import database
from snowav.database.tables import Basins
import numpy as np
from datetime import datetime, timedelta
import seaborn as sns
from climata.usgs import DailyValueIO
from sys import exit

def inflow(self):
    '''
    Data from:
    http://cdec.water.ca.gov/dynamicapp/QueryDaily?s=MIL&end=2018-10-18&span=30days

    To csv:
    /wkspace/results/sj/wy19/friant_dam_wy2019_flow.csv


    From Chad:
    pd.read_excel('Downloads/SJRRP Operations 2018.xlsx',sheet_name= 'Data Entry')
    pd.read_excel('Downloads/SJRRP Operations 2018.xlsx',sheet_name= 'SJRRP Flow Summary')

    '''

    # load inflow data
    inflow_path = self.flow_file
    inflow = pd.read_csv(inflow_path)
    inflow.set_index('DATE / TIME (PST)', inplace=True)
    inflow.index = inflow.index.to_datetime()

    cfs_to_TAF = 1.98/1000

    # get swi data from database
    swi = pd.DataFrame(columns = self.plotorder)

    end_date = inflow.index.max().to_pydatetime() + timedelta(days = 1)

    for bid in self.plotorder:
        r2 = database.database.query(self,
                                    datetime(self.wy-1,10,1),
                                    end_date,
                                    self.run_name,
                                    # 'sj_wy2019_ops',
                                    bid,
                                    'swi_vol')

        v2 = r2[(r2['elevation'] == 'total')]

        for iter,d in enumerate(v2['date_time'].values):
            swi.loc[d,bid] = v2['value'].values[iter]

    # print(inflow.index.max().to_pydatetime(),v2)

    swi.sort_index(inplace=True)
    swi.index = swi.index.to_datetime()
    percent = np.cumsum(inflow['INFLOW CFS'].values*cfs_to_TAF)/np.cumsum(swi['San Joaquin River Basin'].values)*100

    if self.basin == 'SJ':
        percent[0:60] = np.nan
        res = 'Millerton'
    else:
        res = ''

    sns.set_style('white')
    sns.set_context("notebook")

    plt.close(0)
    f = plt.figure(num=0, figsize=(6,4), dpi=self.dpi)
    a = plt.gca()
    b = a.twinx()

    a.plot(swi.index,np.nancumsum(swi['San Joaquin River Basin'].values),'r',label = 'SWI')
    a.plot(inflow.index,inflow['INFLOW CFS'].cumsum().values*cfs_to_TAF,'k',label='{} inflow'.format(res))
    # a.plot(dates,flow, 'g')

    p, = b.plot(inflow.index,
           percent,
           linestyle = ':',
           color = 'b',
           label = 'runoff\nefficiency')

    # p.axes.get_xaxis().set_ticks([])
    # p.axes.get_yaxis().set_ticks([])
    # b.xaxis.set_visible(True)

    a.legend()
    a.set_ylabel('volume [TAF]')

    for tick in a.get_xticklabels():
        tick.set_rotation(30)

    b.set_ylabel('runoff efficiency [%]')
    b.yaxis.label.set_color('b')
    b.tick_params(axis='y', colors=p.get_color())
    b.set_ylim((0,50))

    plt.tight_layout()
    plt.savefig('{}inflow_{}.png'.format(self.figs_path,self.name_append))
