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


def compare_runs(self):

    '''
    Figure not complete for all basins, see bid...

    '''

    start_date = datetime.datetime(self.wy - 1,10,1)
    end_date = datetime.datetime(self.wy,6,10)

    # hack
    bid = self.plotorder[0]

    sns.set_style('darkgrid')
    sns.set_context("notebook")

    for i,var in enumerate(self.plot_variables):

        plt.close(i)
        plt.figure(num = i, figsize = (6,4), dpi = 200)
        ax = plt.gca()

        if 'z' in var:
            lbl = 'mm'
        if ('vol' in var) or 'avail' in var:
            lbl = self.vollbl
        if var == 'density':
            lbl = 'kg/m^3'

        for run in self.plot_runs:
            qry = self.session.query(Results).join(RunMetadata).filter(and_((Results.date_time >= start_date),
                                                      (Results.date_time <= end_date),
                                                      (Results.variable == var),
                                                      (RunMetadata.run_name == run),
                                                      (Results.basin_id == BASINS.basins[bid]['basin_id'])))


            df = pd.read_sql(qry.statement, qry.session.connection())

            if df.empty:
                print('No results from database in compare_runs() query, '
                      'database may not have results for:\n'
                      '{} to {}\n'
                      '{}\n'
                      '{}\n'
                      '{}'.format(start_date,end_date,var,run,bid))

            df.index = df['date_time']
            df['month'] = df.index.to_series().dt.strftime('%b')

            if var == 'swi_vol':
                ax.plot(df[df['elevation'] == 'total']['date_time'],
                        df[df['elevation'] == 'total']['value'].cumsum(),label = run)
            else:
                ax.plot(df[df['elevation'] == 'total']['date_time'],
                        df[df['elevation'] == 'total']['value'],label = run)

            for tick in ax.get_xticklabels():
                tick.set_rotation(30)

            formatter = mdates.DateFormatter('%b')
            ax.xaxis.set_major_formatter(formatter)
            ax.legend()
            ax.set_ylabel('[%s]'%(lbl))
            ax.set_title(var)
            plt.tight_layout()

            plt.savefig('{}{}_{}{}.png'.format(self.figs_path,run,var,self.name_append))

    plt.show()
    # session.close()
