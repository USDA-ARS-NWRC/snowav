import os
import datetime
from sqlalchemy import Column, ForeignKey, Integer, String
from sqlalchemy.ext.declarative import declarative_base
from sqlalchemy import create_engine
from sqlalchemy.orm import sessionmaker
from sqlalchemy.orm import backref
from snowav.database.tables import RunMetadata, Watershed, Basin, Results, VariableUnits, Watersheds, Basins, Base
from sqlalchemy import and_
from snowav import database
from snowav.database.tables import Basins
import pandas as pd
import snowav
import numpy as np
import matplotlib
from matplotlib import pyplot as plt
import seaborn as sns
import copy
import matplotlib.dates as mdates
import datetime


def compare_runs(self):

    '''
    Figure not complete for all basins, see bid...

    '''

    start_date = datetime.datetime(self.wy - 1,10,1)
    end_date = datetime.datetime(self.wy,8,1)

    bid = self.plotorder[0]

    sns.set_style('darkgrid')
    sns.set_context("notebook")

    for i,var in enumerate(self.plot_variables):

        plt.close(i)
        plt.figure(num = i, figsize = self.figsize, dpi = self.dpi)
        ax = plt.gca()

        formatter = mdates.DateFormatter('%b')
        ax.xaxis.set_major_formatter(formatter)

        if 'z' in var:
            lbl = 'mm'
        if ('vol' in var) or 'avail' in var:
            lbl = self.vollbl
        if var == 'density':
            lbl = 'kg/m^3'

        for iters,run in enumerate(self.plot_runs):
            qry = self.session.query(Results).join(RunMetadata).filter(and_(#(Results.date_time >= start_date),
                                                      #(Results.date_time <= end_date),
                                                      (Results.variable == var),
                                                      (RunMetadata.run_name == run),
                                                      (Results.basin_id == Basins.basins[bid]['basin_id'])))

            df = pd.read_sql(qry.statement, qry.session.connection())

            if df.empty:
                print('No results from database in compare_runs() query, '
                      'database may not have results for:\n'
                      '{} to {}\n'
                      '{}\n'
                      '{}\n'
                      '{}'.format(start_date,end_date,var,run,bid))

            df.index = pd.to_datetime(df.date_time)
            df['month'] = df.index.month
            df['day'] = df.index.day
            df['year'] = self.wy
            i = df['month'] >= 10
            df.year[i] = self.wy - 1
            df.index = pd.to_datetime(df[['year','month','day']])

            if var == 'swi_vol':
                ax.plot(df[df['elevation'] == 'total']['value'].cumsum(),
                        label = self.plot_labels[iters])
            else:
                ax.plot(df[df['elevation'] == 'total']['value'],
                        label = self.plot_labels[iters])

            if var == 'swe_vol':
                title = 'Basin Total SWE'
            if var == 'swi_vol':
                title = 'Basin Total SWI'

        if self.flight_dates is not None:
            for i,d in enumerate(self.flight_dates):
                if i == 0:
                    lb = 'flight update'
                else:
                    lb = '__nolabel__'
                ax.axvline(x=d,linestyle = ':',linewidth = 0.75, color = 'k',label = lb)

        for tick in ax.get_xticklabels():
            tick.set_rotation(30)

        formatter = mdates.DateFormatter('%b')
        ax.xaxis.set_major_formatter(formatter)
        ax.legend()
        ax.set_ylabel('[%s]'%(lbl))
        # ax.set_ylim((-5,2600))
        ax.set_xlim((datetime.datetime(self.wy-1,10,1),datetime.datetime(self.wy,9,30)))

        ax.set_title(title)
        plt.tight_layout()

        self._logger.info('saving figure to {}compare_{}_{}.png'.format(self.figs_path,var,self.name_append))
        plt.savefig('{}compare_{}_{}.png'.format(self.figs_path,var,self.name_append))

    # This will need to change if framework.py is updated with different workflow
    # self.session.close()
