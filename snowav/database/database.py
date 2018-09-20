
import os
import datetime
from sqlalchemy import Column, ForeignKey, Integer, String
from sqlalchemy.ext.declarative import declarative_base
from sqlalchemy import create_engine
from sqlalchemy.orm import sessionmaker
from sqlalchemy.orm import backref
from snowav.database.tables import BasinMetadata, Base, Results, RunMetadata, BASINS
from sqlalchemy import and_
import pandas as pd
import smrf
import awsm
import snowav
from snowav import database
import numpy as np


def insert(self, table, values):
    '''
    Inserts standard run results to the database. Uses database session created
    in read_config()

    Args
        table: database table (currently RunMetadata or Results)
        values: dictionary of values

    '''

    dbtable = getattr(snowav.database.tables, table)
    my_dbtable = dbtable()

    for k,v in values.items():
        setattr(my_dbtable, k, v)

    self.session.add(my_dbtable)
    self.session.commit()


def query(self, start_date, end_date, run_name, bid = None, value = None):
    '''
    Query and retrieve results from the database, using database session created
    in read_config().

    Args
        start_date: (datetime)
        end_date: (datetime)
        run_name: identifier for run, specified in config file
        bid: basin id in string format ('Boise River Basin')
        value: value to query (i.e. 'swi_z')

    Returns
        dataframe of query results

    '''

    # Subset query with value and bid if desired
    if (value != None) and (bid != None):
        qry = self.session.query(Results).filter(and_((Results.date_time >= start_date),
                                              (Results.date_time <= end_date),
                                              (Results.run_name == run_name),
                                              (Results.variable == value),
                                              (Results.basin_id == BASINS.basins[bid]['basin_id'])))

    # Otherwise you're gonna get a bunch of stuff...
    else:
        qry = self.session.query(Results).filter(and_((Results.date_time >= start_date),
                                              (Results.date_time <= end_date),
                                              (Results.run_name == run_name)))

    # x = self.session.query(Results).get(1)
    # print(x.run_metadata.data_location)

    df = pd.read_sql(qry.statement, qry.session.connection())
    return df

def delete(self, start_date, end_date, bid, run_name):
    '''
    Delete results from the databaseusing database session created
    in read_config(). Currently this deletes all values with basin_id == bid
    and run_name = run_name within the date range.

    Args
        start_date: (datetime)
        end_date: (datetime)
        bid: basin id in string format ('Boise River Basin')
        run_name: identifier for run, specified in config file

    '''

    try:
        qry = self.session.query(Results).filter(and_((Results.date_time >= start_date),
                                          (Results.date_time <= end_date),
                                          (Results.run_name == run_name),
                                          (Results.basin_id == BASINS.basins[bid]['basin_id'])))

        # Only ever figured out how to do this one record at a time
        for record in qry:
            self.session.delete(record)

        self.session.commit()

    except:
        print('Failed during attempted deletion in database.delete()')


def run_metadata(self):
    '''
    Create run metadata entry for each snowav run, using database session
    created in read_config().

    '''

    qry = self.session.query(RunMetadata)
    df = pd.read_sql(qry.statement, qry.session.connection())

    # Increase each runid by 1
    if not df.run_id.values.size:
        self.runid = 1
    else:
        self.runid = int(df['run_id'].max() + 1)

    values = {
              'run_id':self.runid,
              'basin_id':BASINS.basins[self.plotorder[0]]['basin_id'],
              'run_name':self.run_name,
              'basin_name':self.plotorder[0],
              'description':'',
              'smrf_version':'smrf'+ smrf.__version__,
              'awsm_version':'aswm'+ awsm.__version__,
              'snowav_version':'snowav'+ snowav.__version__,
              'data_type':'',
              'data_location':','.join(self.run_dirs),
              'file_type':'',
              'config_file':self.config_file,
              'proc_time':datetime.datetime.now(),
              'var_units': (self.depthlbl + ', ' + self.vollbl + ', ' + 'MJ'
                            + ', ' + 'kg/m^3'),
              'elev_units':self.elevlbl
                }

    snowav.database.database.insert(self,'RunMetadata',values)


def check_fields(self, start_date, end_date, bid, run_name, value):
    '''
    This functions queries the database and returns True if any value exists
    in the given date range, basin_id = id and run_name = run_name,
    and False if not.

    Args
        start_date: (datetime)
        end_date: (datetime)
        bid: basin id in string format ('Boise River Basin')
        run_name: identifier for run, specified in config file
        value: value to query (i.e. 'swi_z')

    Returns:
        flag: boolean (True if results exist)

    '''

    qry = self.session.query(Results).filter(and_((Results.date_time >= start_date),
                                          (Results.date_time <= end_date),
                                          (Results.run_name == run_name),
                                          (Results.basin_id == BASINS.basins[bid]['basin_id']),
                                          (Results.variable == value))).first()

    if qry is not None:
        flag = True
    else:
        flag = False

    return flag

def write_csv(self):
    '''
    This writes out variables, specified in the config file, from the database
    to the figs_path directory

    '''

    # Initialize
    out = pd.DataFrame(columns = self.plotorder)

    # By variable
    for var in self.write_csv:

        # For all subbasins
        for bid in self.plotorder:
            r = database.database.query(self, datetime.datetime(self.wy-1,10,1),
                                        self.end_date, self.run_name, bid, var)

            # For now, just write out totals
            v = r[(r['elevation'] == 'total')]

            for iter,d in enumerate(v['date_time'].values):
                out.loc[d,bid] = v['value'].values[iter]

            # SWI is cumulative
            if var == 'swi_vol':
                out = out.cumsum()

        if self.vollbl == 'KAF':
            out.to_csv(os.path.join(self.figs_path,var + '_' + self.vollbl + '.csv'))

        if self.vollbl == 'SI':
            out.to_csv(os.path.join(self.figs_path,var + '_' + 'm3' + '.csv'))
