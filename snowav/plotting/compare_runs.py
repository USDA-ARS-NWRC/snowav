import os
import datetime
from sqlalchemy import Column, ForeignKey, Integer, String
from sqlalchemy.ext.declarative import declarative_base
from sqlalchemy import create_engine
from sqlalchemy.orm import sessionmaker
from sqlalchemy.orm import backref
from snowav.database.tables import Basin_Metadata, Base, Results, Run_Metadata, BASINS
from sqlalchemy import and_
from snowav import database
from snowav.database.tables import BASINS
import pandas as pd
import snowav
import numpy as np
from matplotlib import pyplot as plt
import seaborn as sns
import copy
import matplotlib.dates as mdates


def compare_runs(self):

    '''

    '''

    engine = create_engine(self.database)
    connection = engine.connect()
    DBSession = sessionmaker(bind=engine)
    session = DBSession()

    start_date = datetime.datetime(self.wy - 1,10,1)
    end_date = datetime.datetime(self.wy,6,10)

    # hack
    if self.basin == 'BRB':
        bid = 'Boise River Basin'

    sns.set_style('darkgrid')
    sns.set_context("notebook")

    plt.close(0)
    plt.figure(num = 0)
    ax = plt.gca()

    var = self.plot_variables[0]

    for run in self.plot_runs:
        qry = session.query(Results).filter(and_((Results.date_time >= start_date),
                                                  (Results.date_time <= end_date),
                                                  (Results.variable == var),
                                                  (Results.run_name == run),
                                                  (Results.basin_id == BASINS.basins[bid]['basin_id'])))

        df = pd.read_sql(qry.statement, qry.session.connection())

        df.index = df['date_time']
        df['month'] = df.index.to_series().dt.strftime('%b')

        ax.plot(df[df['elevation'] == 'total']['date_time'],
                df[df['elevation'] == 'total']['value'],label = run)

    session.close()

    for tick in ax.get_xticklabels():
        tick.set_rotation(30)

    formatter = mdates.DateFormatter('%b')
    ax.xaxis.set_major_formatter(formatter)
    ax.legend()
    ax.set_ylabel('[%s]'%(self.vollbl))

    # plt.savefig('%sbasin_total_multiyr%s.png'%(snow.figs_path,snow.name_append))
    plt.tight_layout()
    plt.show()
