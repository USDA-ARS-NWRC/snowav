
import os
from datetime import datetime
from sqlalchemy import Column, ForeignKey, Integer, String, create_engine, and_
from sqlalchemy.ext.declarative import declarative_base
from sqlalchemy.orm import backref, sessionmaker
from snowav.database.tables import Base, RunMetadata, Watershed, Basin, Results, VariableUnits
import pandas as pd
import snowav
import numpy as np
import copy
import urllib.parse
import mysql.connector
from sys import exit
import warnings

# Fix these two by pulling smrf and awsm versions from netcdf
try:
    import smrf
    smrf_version = smrf.__version__
except:
    print('Could not import smrf, database smrf version will be "unknown"')
    smrf_version = 'unknown'

try:
    warnings.filterwarnings("ignore")
    import awsm
    warnings.filterwarnings("default")
    awsm_version = awsm.__version__

except:
    print('Could not import awsm, database awsm version will be "unknown"')
    awsm_version = 'unknown'

def make_session(connector):
    '''
    Make snowav database session.

    Args
    -------
    connector : str
        database connector, either path to sqlite database or mysqlconnector
        string

    Returns
    -------
    session : object
        snowav database session
    '''

    try:
        engine = create_engine(connector)
        Base.metadata.create_all(engine)

    except:
        raise Exception('Failed to make database connection with '
                        '{}'.format(connector))

    DBSession = sessionmaker(bind=engine)
    session = DBSession()

    return session

def insert(connector, table, values):
    '''
    Inserts standard run results to the database.

    Args
    --------
    session : object
        database session
    table : str
        database table (RunMetadata, Results, or VariableUnits)
    values : dict
        dictionary of values

    '''

    dbtable = getattr(snowav.database.tables, table)
    my_dbtable = dbtable()

    for k,v in values.items():
        setattr(my_dbtable, k, v)

    session = make_session(connector)
    session.add(my_dbtable)
    session.commit()
    session.close()

def collect(connector, plotorder, basins, start_date, end_date, value,
            run_name, edges, method):
    '''

    Collect snowav database results into our standard dataframe with elevation
    rows and subbasin columns.

    Args
    ---------
    session : object
        sqlalchemy snowav database session
    plotorder : str
        list of basins
    start_date : datetime
    end_date : datetime
    value : str
        database value identifier
    run_name : str
        database run name identifier
    edges : list
        elevation band, can include 'total'
    method : str
        options [end sum difference daily]

    Returns
    ----------
    df : DataFrame
        standard subbasin columns and elevation labels

               Extended Tuolumne  Tuolumne  Cherry Creek  Eleanor
        3000               0.000     0.000         0.000    0.000
        ...                  ...       ...           ...      ...

    '''

    value_options = ['swe_z','swe_vol','density','precip_z','precip_vol',
                     'rain_z','rain_vol','swe_avail','swe_unavail','coldcont',
                     'swi_z','swi_vol','depth','snow_line','mean_air_temp',
                     'evap_z', 'L_v_E','melt','lwc','temp_surface',
                     'temp_lower','temp_bulk','depth_lower_layer','h20_sat',
                     'R_n','H','G','M','delta_Q']

    method_options = ['sum','difference','daily','end']

    if value not in value_options:
        raise Exception('value={} call to database.collect() '.format(value) +
                        'must be one of {}'.format(value_options))

    if method not in method_options:
        raise Exception('method={} call to database collect() '.format(method) +
                        'must be one of {}'.format(method_options))

    if type(plotorder) != list:
        plotorder = [plotorder]

    if edges == 'total':
        edges = ['total']

    if method == 'daily':
        df = pd.DataFrame(index = ['date_time'], columns=plotorder)

    else:
        df = pd.DataFrame(np.nan, index=edges, columns=plotorder)

    for bid in plotorder:
        results = query(connector, start_date, end_date, run_name, basins, bid, value)

        if results.empty:
            raise Exception('Empty results for database query in '
                            'database.collect(), check config file '
                            'start_date, end_date, and that the run_name={} '
                            'has results for {} to {}'. format(run_name,
                            start_date, end_date))

        if method == 'daily':
            e = results[(results['elevation'] == 'total')
                        & (results['date_time'] >= start_date)
                        & (results['date_time'] <= end_date)]
            if e.empty:
                raise Exception('Empty results for database query in '
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
                    df[bid] = e['value']

            df.sort_index(inplace=True)

        if method == 'end':
            for elev in edges:
                e = results[(results['elevation'] == str(elev))
                            & (results['date_time'] == end_date)]
                if e.empty:
                    raise Exception('Empty results for database query in '
                                    'database.collect(), for run_name={}, '
                                    'elev={}, {} to {}'. format(run_name, elev,
                                    start_date, end_date))
                else:
                    if e['value'].values[0] is None:
                        df.loc[elev,bid] = np.nan
                    else:

                        # The database can occasionally get multiple values
                        # for the same record if it doesn't exit cleanly
                        if len(e['value'].values) == 1:
                            df.loc[elev,bid] = e['value'].values
                        else:
                            raise Exception('Multiple database entries for a '
                                'single field, consider running with '
                                '[database] overwrite: True')

        if method == 'difference':
            for elev in edges:
                s = results[(results['elevation'] == str(elev))
                            & (results['date_time'] == start_date)]
                e = results[(results['elevation'] == str(elev))
                            & (results['date_time'] == end_date)]
                if e.empty or s.empty:
                    raise Exception('Empty results for database query in '
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
                    raise Exception('Empty results for database query in '
                                    'database.collect(), for run_name={}, '
                                    'elev={}, {} to {}'. format(run_name, elev,
                                    start_date, end_date))
                else:
                    df.loc[elev,bid] = e['value'].sum(skipna=False)

    df = df[plotorder]

    return df

def query(connector, start_date, end_date, run_name, basins, bid = None,
          value = None, rid = None):
    '''
    Retrieve results from snowav database.

    Args
    --------
    start_date : datetime
    end_date : datetime
    run_name : str
        identifier for run, specified in config file
    bid : str
        basin id in string format ('Boise River Basin')
    value : str
        value to query (i.e. 'swi_z')
    rid : run_id

    Returns
    ---------
    df : DataFrame
        dataframe of query results

    '''

    basin_id = int(basins[bid]['basin_id'])
    session = make_session(connector)

    if (value != None) and (bid != None):
        qry = session.query(Results).join(RunMetadata).filter(and_(
                        (Results.date_time >= start_date),
                        (Results.date_time <= end_date),
                        (RunMetadata.run_name == run_name),
                        (Results.variable == value),
                        (Results.basin_id == basin_id)))

    elif (value == None) and (bid != None) and rid != None:
        qry = session.query(Results).join(RunMetadata).filter(and_(
                                            (Results.date_time >= start_date),
                                            (Results.date_time <= end_date),
                                            (Results.run_id == int(rid)),
                                            (RunMetadata.run_name == run_name),
                                            (Results.basin_id.in_(bids))))

    else:
        qry = session.query(Results).join(RunMetadata).filter(and_(
                                            (Results.date_time >= start_date),
                                            (Results.date_time <= end_date),
                                            (RunMetadata.run_name == run_name)))

    df = pd.read_sql(qry.statement, qry.session.connection())

    session.close()

    return df

def delete(connector, basins, start_date, end_date, bid, run_name):
    '''
    Delete results from the database. This deletes values with basin_id == bid
    and run_name = run_name within the date range. If there are no more
    remaining records for that run_id (the entire run was deleted), then the
    RunMetadata and VariableUnits data are also removed.

    Note: deletion order matters because of table relationships.

    Args
    --------
    start_date : datetime
    end_date : datetime
    bid : str
        basin name
    run_name : str
        identifier for run, specified in config file

    Returns
    --------
    logger : list

    '''

    logger = []
    basin_id = int(basins[bid]['basin_id'])
    session = make_session(connector)

    logger.append(' Deleting existing records for {}, {}, {} '.format(
                  bid, run_name, start_date.date()))

    # Get the run_id
    qry = session.query(Results).join(RunMetadata).filter(and_(
                                        (Results.date_time >= start_date),
                                        (Results.date_time <= end_date),
                                        (RunMetadata.run_name == run_name),
                                        (Results.basin_id == basin_id)))

    df = pd.read_sql(qry.statement, qry.session.connection())
    unique_runs = df.run_id.unique()

    if not df.empty:
        for r in unique_runs:
            session.query(Results).filter(and_((Results.date_time >= start_date),
                                    (Results.date_time <= end_date),
                                    (Results.run_id == int(r)),
                                    (Results.basin_id == basin_id))).delete()

            session.commit()

        # Query again, if no results from those run_id exist still, also remove
        # the metadata and variableUnits
        qrym = session.query(Results).filter(Results.run_id == int(r))
        dfm = pd.read_sql(qrym.statement, qrym.session.connection())

        if dfm.empty:
            logger.append(' Deleting RunMetadata run_name={}, run_id={}, from {} '
                          'to {}'.format(run_name, str(r), start_date.date(),
                          end_date.date()))

            session.query(VariableUnits).filter(VariableUnits.run_id == int(r)).delete()
            session.query(RunMetadata).filter(RunMetadata.run_id == int(r)).delete()
            session.commit()

    session.close()

    return logger

def create_tables(database, plotorder):
    '''
    Creates database Watersheds and Basins tables. This is currently only
    being called on database creation.

    Args
    -------
    database : str
        database identifier
    plotorder : list

    Returns
    --------
    logger : list


    '''

    logger = []

    engine = create_engine(database)
    Base.metadata.create_all(engine)
    DBSession = sessionmaker(bind=engine)
    session = DBSession()

    logger.append(' Adding watershed {} to database'.format(plotorder[0]))

    # Initialize watershed
    wval = Watershed(watershed_id = 1, watershed_name = plotorder[0])
    session.add(wval)
    session.commit()

    # Initialize basins within the watershed
    for i,name in enumerate(plotorder):
        logger.append(' Adding basin {} to database'.format(name))

        bval = Basin(watershed_id = 1,
                     basin_id = i + 1,
                     basin_name = name)

        session.add(bval)

    session.commit()
    session.close()

    logger.append(' Completed initializing basins and watersheds')

    return logger

def basins_dict():
    '''
    Make the basins lookup dictionary.

    '''

def run_metadata(self, run_name):
    '''
    Create database RunMetadata for each snowav run.

    Args
    ------
    run_name : str

    '''

    watershed_id = int(self.basins[self.plotorder[0]]['watershed_id'])

    session = make_session(self.connector)
    qry = session.query(RunMetadata)
    df = pd.read_sql(qry.statement, qry.session.connection())
    session.close()

    # Increase each runid by 1
    if not df.run_id.values.size:
        self.run_id = 1

    else:
        self.run_id = int(df['run_id'].max() + 1)

    # in most cases it's not practical to save all run_dirs
    trdir = ','.join(self.run_dirs)
    if len(trdir) > 1200:
        trdir = self.run_dirs[0]

    values = {'run_id':int(self.run_id),
              'run_name':run_name,
              'watershed_id':int(watershed_id),
              'pixel':int(self.pixel),
              'description':'',
              'smrf_version':'smrf'+ smrf_version,
              'awsm_version':'aswm'+ awsm_version,
              'snowav_version':'snowav'+ snowav.__version__,
              'data_type':'',
              'data_location':trdir,
              'file_type':'',
              'config_file':self.config_file,
              'proc_time':datetime.now() }

    insert(self.connector, 'RunMetadata', values)

    self.vid = {}
    for v in [*self.vars.keys()]:
        u = self.vars[v]['units']
        if (('vol' in v) or ('avail' in v)):
            u = self.units
        if (('z' in v) or (v == 'depth')):
            u = self.depthlbl
        if v == 'density':
            u = 'kg m^-3'
        if v == 'coldcont':
            u = 'MJ'

        variables = {'run_id':self.run_id,
                     'variable':v,
                     'unit':u,
                     'name':self.vars[v]['description']}

        insert(self.connector,'VariableUnits',variables)

        # After inserting into VariableUnits, get the existing id
        session = make_session(self.connector)
        qry = session.query(VariableUnits).filter(VariableUnits.run_id == self.run_id)
        session.close()
        df = pd.read_sql(qry.statement, qry.session.connection())
        self.vid[v] = df[df['variable'] == v]['id'].values[0]

    # snow_line
    variables = {'run_id':self.run_id,
                 'variable':'snow_line',
                 'unit':self.elevlbl,
                 'name':'snow_line'}

    insert(self.connector,'VariableUnits',variables)
    session = make_session(self.connector)
    qry = session.query(VariableUnits).filter(VariableUnits.run_id == self.run_id)
    session.close()
    df = pd.read_sql(qry.statement, qry.session.connection())
    self.vid['snow_line'] = df['id'].values[0]

    # Add post-processing values
    sum_vals = {'swi_vol_wy':'surface water input volume, water year total',
                'swi_z_wy':'surface water input depth, water year total',
                'evap_z_wy':'evaporation depth, water year total'}

    for v in sum_vals.keys():
        if ('vol' in v):
            u = self.units
        if ('z' in v) :
            u = self.depthlbl

        variables = {'run_id':self.run_id,
                     'variable':v,
                     'unit':u,
                     'name':sum_vals[v]}

        insert(self.connector,'VariableUnits',variables)

        session = make_session(self.connector)
        qry = session.query(VariableUnits).filter(VariableUnits.run_id == self.run_id)
        session.close()
        df = pd.read_sql(qry.statement, qry.session.connection())

        self.vid[v] = df[df['variable'] == v]['id'].values[0]


def connect(sqlite = None, sql = None, plotorder = None, user = None,
            password = None, host = None, port = None, convert = False,
            add = False):
    '''
    This establishes a connection with a database for results. If the specified
    sqlite database doesn't exist, it will be created.

    Note: consider moving creation of mysql database to Makefile

    Args
    --------
    sqlite : str
        path to sqlite database
    sql : str
        snowav mysql identifier
    plotorder : list
        list of basins
    user : str
        database user
    password : str
        database password
    host : int
        database host
    port : int
        database port

    Returns
    ---------
    basins : dict
        dictionary of snowav basins and ids for database connection
    logger : list

    '''

    logger = []

    if sqlite is not None:
        fp = urllib.parse.urlparse(sqlite)

        if (not os.path.isfile(fp.path)):
            database = sqlite
            logger.append(' Creating {} for results'.format(database))
            log = create_tables(database, plotorder)

            for out in log:
                logger.append(out)

            engine = create_engine(sqlite)

        else:
            engine = create_engine(sqlite)

        DBSession = sessionmaker(bind=engine)
        session = DBSession()
        logger.append(' Using {} for results'.format(sqlite))
        connector = sqlite

    if (sql is not None) and (sqlite is None):

        try:
            cnx = mysql.connector.connect(user=user,
                                          password=password,
                                          host=host,
                                          port=port)
        except:
            print('Failed attempting to make database connection with '
                  '{}:{}@{}/{}\nCheck config options in [database] '
                  'section'.format(user,password,host,port))
            exit()
        cursor = cnx.cursor()

        # Check if database exists, create if necessary
        query = ("SHOW DATABASES")
        cursor.execute(query)
        dbs = cursor.fetchall()

        try:
            dbs = [i[0].decode("utf-8") for i in dbs]
        except:
            dbs = [i[0] for i in dbs]

        db_engine = 'mysql+mysqlconnector://{}:{}@{}:{}/{}'.format(user,
                                                                password,
                                                                host,
                                                                port,
                                                                sql)

        # If the database doesn't exist, create it, otherwise connect
        if (sql not in dbs):
            logger.append(' Specified mysql database {} '.format(sql) +
                  ' does not exist, it is being created...')
            query = ("CREATE DATABASE {};".format(sql))
            cursor.execute(query)
            log = create_tables(db_engine, plotorder)

            for out in log:
                logger.append(out)

            logger.append(log)

        else:
            query = ("SHOW DATABASES")
            cursor.execute('USE {}'.format(sql))
            cursor.execute('SHOW TABLES')
            tbls = cursor.fetchall()

            # This hasn't been tested with new docker mysql instance, could be
            # a point of failure...
            if tbls is None:
                log = create_tables(db_engine, plotorder)

                for out in log:
                    logger.append(out)

                logger.append(log)

            try:
                logger.append(' Using database connection {}@{}/{} for '
                              'results'.format(user, host, sql))
                engine = create_engine(db_engine)
                Base.metadata.create_all(engine)
                DBSession = sessionmaker(bind=engine)
                session = DBSession()

            except:
                logger.append(' Failed trying to make database connection '
                      'to {}'.format(sql))

        cursor.close()
        cnx.close()

        connector = db_engine

    qry = session.query(Watershed)
    ws = pd.read_sql(qry.statement, qry.session.connection())

    qry = session.query(Basin)
    bs = pd.read_sql(qry.statement, qry.session.connection())

    if convert:
        plotorder[0] = convert_watershed_names(plotorder[0])

    if plotorder[0] in ws['watershed_name'].values:
        f = ws[ws['watershed_name'].values == plotorder[0]]
        wid = f['watershed_id'].values[0]

    else:
        logger.append(' {} not found as a watershed in '.format(plotorder[0]) +
                      'existing database, it is being added')

        if not ws.empty:
            wid = np.max(ws['watershed_id'].values) + 1
        else:
            wid = 1

        wval = Watershed(watershed_id = int(wid),
                         watershed_name = plotorder[0])
        if add:
            session.add(wval)

        else:
            logger.append(' Config option [database] add_basins: False, so '
                          '{} will not be added'.format(plotorder[0]))

            print('WARNING! Given current config options and topo.nc file, '
                  'trying to use {} as a watershed, but it does not currently '
                  'exist as a watershed in the database...\nCorrect [snowav] '
                  'masks or consider changing [database] add_basins: True '
                  'if you know what you are doing.'.format(plotorder[0]))

        session.commit()

    basins = {}
    watersheds = {plotorder[0]:{'watershed_id':wid, 'watershed_name':'',
                                'basins': '', 'shapefile':''}}

    # Initialize basins within the watershed
    for i,name in enumerate(plotorder):

        # This logic doesn't work if trying to 'move' a subbasin from one watershed to another,
        # i.e., 'Cherry Creek' going from part of the Extended Tuolumne to Tuolumne River
        # Basin in wy2020
        if name in bs['basin_name'].values:
            bid = int(bs[(bs['watershed_id'] == wid) & (bs['basin_name'] == name)]['basin_id'].values[0])

        else:
            logger.append(' {} not found as a basin in '.format(name) +
                          'existing database, it is being added')

            if not bs.empty:
                bid = np.max(bs['basin_id'].values) + 1
            else:
                bid = 1

            bval = Basin(watershed_id = int(wid),
                         basin_id = int(bid + i),
                         basin_name = name)

            if add:
                session.add(bval)

            else:
                logger.append(' Config option [database] add_basins: False, so '
                              '{} will not be added'.format(name))

                print('WARNING! Given current config options and topo.nc file, '
                      'trying to use {} as a basin, but it does not currently '
                      'exist as a basin in the database...\nConsider '
                      'changing [database] add_basins: True if you know what '
                      'you are doing.'.format(name))

        basins[name] = {'watershed_id':wid, 'basin_id':bid}

    session.commit()
    session.close()

    return basins, connector, logger

def convert_watershed_names(name):
    '''
    Convert watershed from topo.nc mask to existing snowav database version.

    '''
    convert_list = ['Boise River Basin','Lakes Basin','Merced River Basin',
                    'Kaweah River Basin','Kings River Basin']

    if name in convert_list:
        watershed = name.split(' ')[0]

    elif name == 'San Joaquin River Basin':
        watershed = 'San Joaquin'

    else:
        watershed = name

    return watershed
