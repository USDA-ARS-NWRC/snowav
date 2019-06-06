
import os
import datetime
from sqlalchemy import Column, ForeignKey, Integer, String
from sqlalchemy.ext.declarative import declarative_base
from sqlalchemy import create_engine
from sqlalchemy.orm import sessionmaker
from sqlalchemy.orm import backref
from snowav.database.tables import Base, RunMetadata, Watershed, Basin, Results, VariableUnits, Watersheds, Basins
from sqlalchemy import and_
import pandas as pd
import smrf
import awsm
import snowav
from snowav import database
import numpy as np
import copy
import urllib.parse
import mysql.connector


def insert(self, table, values):
    '''
    Inserts standard run results to the database. Uses database session created
    in read_config()

    Args
        table: database table (RunMetadata, Results, or VariableUnits)
        values: dictionary of values

    '''

    dbtable = getattr(snowav.database.tables, table)
    my_dbtable = dbtable()

    for k,v in values.items():
        setattr(my_dbtable, k, v)

    self.session.add(my_dbtable)
    self.session.commit()

def collect(self, plotorder, start_date, end_date, value, run_name, edges, method):
    '''
    Collect snowav database results into our standard dataframe
    with elevation rows and subbasin columns, with end date,
    difference between dates, or sum between dates.

           Extended Tuolumne  Tuolumne  Cherry Creek  Eleanor
    3000      0.008              0.008         0.000    0.000

    '''

    if type(plotorder) != list:
        plotorder = [plotorder]

    if edges == 'total':
        edges = ['total']

    if method == 'daily':
        df = pd.DataFrame(index = ['date_time'], columns=plotorder)

    else:
        df = pd.DataFrame(np.nan, index=edges, columns=plotorder)

    for bid in plotorder:
        results = query(self, start_date, end_date, run_name, bid, value)

        if results.empty:
            raise IOError('Empty results for database query in '
                          'database.collect(), check config file '
                          'start_date, end_date, and that the run_name={} '
                          'has results for {} to {}'. format(run_name,
                          start_date, end_date))

        if method == 'daily':
            e = results[(results['elevation'] == 'total')
                        & (results['date_time'] >= start_date)
                        & (results['date_time'] <= end_date)]
            if e.empty:
                raise IOError('Empty results for database query in '
                              'database.collect(), for run_name={}, '
                              'elev={}, {} to {}'. format(run_name, elev,
                              start_date, end_date))
            else:
                e = e.set_index('date_time')
                e.sort_index(inplace=True)
                if bid == plotorder[0]:
                    df = e[['value']].copy()
                    df = df.rename(columns={'value':bid})

                else:
                    df.loc[e.index,bid] = e['value'].values

            df.sort_index(inplace=True)

        if method == 'end':
            for elev in edges:
                e = results[(results['elevation'] == str(elev))
                            & (results['date_time'] == end_date)]
                if e.empty:
                    raise IOError('Empty results for database query in '
                                  'database.collect(), for run_name={}, '
                                  'elev={}, {} to {}'. format(run_name, elev,
                                  start_date, end_date))
                else:
                    df.loc[elev,bid] = e['value'].values

        if method == 'difference':
            for elev in edges:
                s = results[(results['elevation'] == str(elev))
                            & (results['date_time'] == start_date)]
                e = results[(results['elevation'] == str(elev))
                            & (results['date_time'] == end_date)]
                if e.empty or s.empty:
                    raise IOError('Empty results for database query in '
                                  'database.collect(), for run_name={}, '
                                  'elev={}, {} to {}'. format(run_name, elev,
                                  start_date, end_date))
                else:
                    df.loc[elev,bid] = np.nansum(e['value'].values-s['value'].values)

        if method == 'sum':
            for elev in edges:
                e = results[(results['elevation'] == str(elev))
                            & (results['date_time'] >= start_date)
                            & (results['date_time'] <= end_date)]
                if e.empty:
                    raise IOError('Empty results for database query in '
                                  'database.collect(), for run_name={}, '
                                  'elev={}, {} to {}'. format(run_name, elev,
                                  start_date, end_date))
                else:
                    df.loc[elev,bid] = e['value'].sum()

    return df

def query(self, start_date, end_date, run_name, bid=None, value=None, rid=None):
    '''
    Query and retrieve results from the database, using database session created
    in read_config().

    Args
        start_date: (datetime)
        end_date: (datetime)
        run_name: identifier for run, specified in config file
        bid: basin id in string format ('Boise River Basin')
        value: value to query (i.e. 'swi_z')
        rid: run_id

    Returns
        df: dataframe of query results

    '''

    # Subset query with value and bid if desired
    if (value != None) and (bid != None):
        qry = self.session.query(Results).join(RunMetadata).filter(and_(
                        (Results.date_time >= start_date),
                        (Results.date_time <= end_date),
                        (RunMetadata.run_name == run_name),
                        (Results.variable == value),
                        (Results.basin_id == Basins.basins[bid]['basin_id'])))

    # Report queries here
    elif (value == None) and (bid != None) and rid != None:
        bids = [Basins.basins[n]['basin_id'] for n in bid]
        qry = self.session.query(Results).join(RunMetadata).filter(and_(
                                            (Results.date_time >= start_date),
                                            (Results.date_time <= end_date),
                                            (Results.run_id == rid),
                                            (RunMetadata.run_name == run_name),
                                            (Results.basin_id.in_(bids))))

    else:
        qry = self.session.query(Results).join(RunMetadata).filter(and_(
                                            (Results.date_time >= start_date),
                                            (Results.date_time <= end_date),
                                            (RunMetadata.run_name == run_name)))

    # x = self.session.query(Results).get(1)
    # print(x.runmetadata.run_name)
    df = pd.read_sql(qry.statement, qry.session.connection())

    return df

def delete(self, start_date, end_date, bid, run_name):
    '''
    Delete results from the database using database session created
    in read_config(). This deletes values with basin_id == bid
    and run_name = run_name within the date range. If there are no more
    remaining records for that run_id (the entire run was deleted), then the
    RunMetadata and VariableUnits data are also removed.

    Args
        start_date: (datetime)
        end_date: (datetime)
        bid: basin id in string format ('Boise River Basin')
        run_name: identifier for run, specified in config file

    Table deletion order matters, by relationship.

    '''
    try:
        # print('Deleting records from bid={}, run_name={}, from {} '
        #       'to {}'.format(bid, run_name,start_date.date(),end_date.date()))
        self._logger.info(' Deleting records from bid={}, run_name={}, from {} '
              'to {}'.format(bid, run_name,start_date.date(),end_date.date()))

        # Since we know run_name, but not run_id, get the run_id
        qry = self.session.query(Results).join(RunMetadata).filter(and_(
                        (Results.date_time >= start_date),
                        (Results.date_time <= end_date),
                        (RunMetadata.run_name == run_name),
                        (Results.basin_id == Basins.basins[bid]['basin_id'])))

        df = pd.read_sql(qry.statement, qry.session.connection())

        # Delete by date range, run_id, and basin
        for r in df.run_id.unique():
            self.session.query(Results).filter(and_(
                (Results.date_time >= start_date),
                (Results.date_time <= end_date),
                (Results.run_id == int(r)),
                (Results.basin_id == Basins.basins[bid]['basin_id']))).delete()

            self.session.flush()

            # Query again, if no results from those run_id exist still, also remove
            # the metadata and variableUnits
            qry2 = self.session.query(Results).filter(Results.run_id == int(r))
            df2 = pd.read_sql(qry2.statement, qry2.session.connection())

            if df2.empty:
                # print('Deleting RunMetadata run_name={}, run_id={}, from {} '
                #       'to {}'.format(run_name,str(r),start_date.date(),
                #       end_date.date()))
                self._logger.info(' Deleting RunMetadata run_name={}, run_id={},'
                                  ' from {} to {}'.format(run_name,
                                                          str(r),
                                                          start_date.date(),
                                                          end_date.date()))

                q1 = self.session.query(VariableUnits).\
                                filter(VariableUnits.run_id == int(r)).first()
                self.session.delete(q1)

                q2 = self.session.query(RunMetadata).\
                                filter(RunMetadata.run_id == int(r)).first()
                self.session.delete(q2)

                self.session.flush()

        self.session.commit()

    except:
        print('Failed attempting to delete records from run_name={}, bid = {}, '
              'from {} to {}'.format(run_name,bid,
                                     start_date.date(),
                                     end_date.date()))

def create_tables(self=None, url=None):
    '''
    This function creates the Watersheds and Basins tables in the database
    from classes defined in /snowav/database/tables.py

    '''

    # Make database connection for duration of snowav processing
    if url is not None:
        engine = create_engine(url)

    else:
        engine = create_engine(self.database)

    Base.metadata.create_all(engine)
    DBSession = sessionmaker(bind=engine)
    session = DBSession()

    # Initialize watersheds
    for ws in Watersheds.watersheds:
        wsid = Watersheds.watersheds[ws]['watershed_id']
        wsn = Watersheds.watersheds[ws]['watershed_name']
        wval = Watershed(watershed_id = wsid, watershed_name = wsn)
        session.add(wval)

        # Initialize basins within the watershed
        for bid in Basins.basins:
            if Basins.basins[bid]['watershed_id'] == wsid:
                bval = Basin(watershed_id = wsid,
                             basin_id = Basins.basins[bid]['basin_id'],
                             basin_name = Basins.basins[bid]['basin_name'])
                session.add(bval)

    session.commit()

    if self is not None:
        self.session = session

    print('Completed initializing basins and watersheds')

def run_metadata(self, forecast=None):
    '''
    Create run metadata entry for each snowav run, using database session
    created in read_config().

    '''

    qry = self.session.query(RunMetadata)
    df = pd.read_sql(qry.statement, qry.session.connection())

    # Increase each runid by 1
    if not df.run_id.values.size:
        runid = 1
    else:
        runid = int(df['run_id'].max() + 1)

    if forecast is not None:
        run_name = self.for_run_name
        self.for_run_id = runid

    else:
        run_name = self.run_name
        self.run_id = runid

    # From RunMetadata table, snowav/database/tables.py
    rdir_strlim = 1200
    trdir = ','.join(self.run_dirs)
    if len(trdir) > rdir_strlim:
        trdir = self.run_dirs[0]

    values = {
              'run_id':runid,
              'run_name':run_name,
              'watershed_id':Watersheds.watersheds[self.plotorder[0]]['watershed_id'],
              'pixel':self.pixel,
              'description':'',
              'smrf_version':'smrf'+ smrf.__version__,
              'awsm_version':'aswm'+ awsm.__version__,
              'snowav_version':'snowav'+ snowav.__version__,
              'data_type':'',
              'data_location':trdir,
              'file_type':'',
              'config_file':self.config_file,
              'proc_time':datetime.datetime.now()
                }

    snowav.database.database.insert(self, 'RunMetadata', values)

    # self.vars is defined in read_config(), and is also used in process()
    self.vid = {}
    for v in self.vars.keys():
        if (('vol' in v) or ('avail' in v)):
            u = self.units
        if (('z' in v) or (v == 'depth')):
            u = self.depthlbl
        if v == 'density':
            u = 'kg m^-3'
        if v == 'coldcont':
            u = 'MJ'

        variables = {
                     'run_id':runid,
                     'variable':v,
                     'unit':u,
                     'name':self.vars[v]
                    }

        snowav.database.database.insert(self,'VariableUnits',variables)

        # After inserting into VariableUnits, get the existing id
        x = self.session.query(VariableUnits).all()
        self.vid[v] = copy.deepcopy(x[-1].id)

    # Add post-processing values
    sum_vals = {'swi_vol_wy':'surface water input volume, water year total',
                'swi_z_wy':'surface water input depth, water year total',
                'evap_z_wy':'evaporation depth, water year total'}

    for v in sum_vals.keys():
        if ('vol' in v):
            u = self.units
        if ('z' in v) :
            u = self.depthlbl

        variables = {
                     'run_id':runid,
                     'variable':v,
                     'unit':u,
                     'name':sum_vals[v]
                    }

        snowav.database.database.insert(self,'VariableUnits',variables)

        # After inserting into VariableUnits, get the existing id
        x = self.session.query(VariableUnits).all()
        self.vid[v] = copy.deepcopy(x[-1].id)

def check_fields(self,start_date,end_date,bid,run_name,value,forecast=None):
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

    '''

    qry = self.session.query(Results).join(RunMetadata).filter(and_(
                        (Results.date_time >= start_date),
                        (Results.date_time <= end_date),
                        (RunMetadata.run_name == run_name),
                        (Results.basin_id == Basins.basins[bid]['basin_id']),
                        (Results.variable == value))).first()

    if forecast is None:
        if qry is not None:
            self.pflag = True
        else:
            self.pflag = False

    if forecast is not None:
        if qry is not None:
            self.for_pflag = True
        else:
            self.for_pflag = False


def connect(self):
    '''
    This establishes a connection with a database for results. If the specified
    sqlite database doesn't exist, it will be created.

    '''

    if self.sqlite is not None:
        fp = urllib.parse.urlparse(self.sqlite)

        if (not os.path.isfile(fp.path)):
            self.database = self.sqlite
            self.tmp_log.append('Creating {} for results'.format(self.database))
            create_tables(self)
            engine = create_engine(self.sqlite)

        else:
            engine = create_engine(self.sqlite)

        DBSession = sessionmaker(bind=engine)
        self.session = DBSession()
        self.tmp_log.append('Using {} for results...'.format(self.sqlite))

    if (self.mysql is not None) and (self.sqlite is None):

        cnx = mysql.connector.connect(user=self.db_user,
                                      password=self.db_password,
                                      host=self.db_host,
                                      port=self.db_port)
        cursor = cnx.cursor()

        # Check if database exists, create if necessary
        query = ("SHOW DATABASES")
        cursor.execute(query)
        dbs = cursor.fetchall()

        try:
            dbs = [i[0].decode("utf-8") for i in dbs]
        except:
            dbs = [i[0] for i in dbs]

        db_engine = 'mysql+mysqlconnector://{}:{}@{}/{}'.format(self.db_user,
                                                                self.db_password,
                                                                self.db_host,
                                                                self.mysql)

        # If the database doesn't exist, create it, otherwise connect
        if (self.mysql not in dbs):
            self.tmp_log.append('Specified mysql database {} '.format(self.mysql) +
                  ' does not exist, it is being created...')
            query = ("CREATE DATABASE {};".format(self.mysql))
            cursor.execute(query)
            create_tables(self, url = db_engine)

        else:
            query = ("SHOW DATABASES")
            cursor.execute('USE {}'.format(self.mysql))
            cursor.execute('SHOW TABLES')
            tbls = cursor.fetchall()

            # This hasn't been tested with new docker mysql instance, could be
            # a point of failure...
            if tbls is None:
                create_tables(self, url = db_engine)

            try:
                engine = create_engine(db_engine)
                Base.metadata.create_all(engine)
                DBSession = sessionmaker(bind=engine)
                self.session = DBSession()

            except:
                self.tmp_log.append('Failed trying to make database connection '
                      'to {}'.format(self.mysql))

        cursor.close()
        cnx.close()

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
            out.to_csv(os.path.join(self.figs_path,var+'_'+self.vollbl+'.csv'))

        if self.vollbl == 'SI':
            out.to_csv(os.path.join(self.figs_path,var+'_'+'m3'+'.csv'))
