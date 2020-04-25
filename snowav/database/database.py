from datetime import datetime
from itertools import count, filterfalse
import mysql.connector
import numpy as np
import os
import pandas as pd
import sqlalchemy as sa
from sqlalchemy import create_engine, and_
from sqlalchemy.orm import sessionmaker
from sys import exit
import urllib.parse
import warnings

from snowav.database import tables as ta
import snowav
from snowav.database.tables import Base, RunMetadata, Watershed, Basin, \
    Results, VariableUnits

# Fix these two by pulling smrf and awsm versions from netcdf
try:
    import smrf

    smrf_version = smrf.__version__
except:
    smrf_version = 'unknown'

try:
    warnings.filterwarnings("ignore")
    import awsm

    warnings.filterwarnings("default")
    awsm_version = awsm.__version__
except:
    awsm_version = 'unknown'


class Database(object):
    """ Database class for snowav.

    Args
    ------
    db_type {string}: database type
    db_name {string}: mysql database name, if using db_type='sql'
    user {string}: database user, if using db_type='sql'
    password {string}: database password, if using db_type='sql'
    host {string}: database host, if using db_type='sql'
    port {string}: database port, if using db_type='sql'
    database {string}: sqlite database file path, if using db_type='sqlite'
    """

    def __init__(self, db_type='sql', db_name='snowav', user=None,
                 password=None, host=None, port=None, database=None):

        db_type_options = ['sql', 'sqlite']
        if db_type not in db_type_options:
            raise ValueError('Database input "db_type" must be one of '
                             '{}'.format(db_type_options))

        if db_type == 'sql':
            input_list = [user, password, host, port]
            if sum([x is not None for x in input_list]) != len(input_list):
                raise Exception('Database input db_type "sql" requires '
                                '"user", "password", "host", and "port"')

        if db_type == 'sqlite' and database is None:
            raise Exception('must supply "database" with db_type=sqlite')

        self.db_type = db_type
        self.user = user
        self.password = password
        self.host = host
        self.port = port
        self.db_name = db_name
        self.database = database

    def assign_vars(self, kwargs):
        """ Assign Database attributes.

        Args
        ------
        kwargs {dict}: {item, value}
        """

        if not isinstance(kwargs, dict):
            raise TypeError("kwargs must be dict")

        for item, v in kwargs.items():
            setattr(self, item, v)

    def make_connection(self, logger=None):
        """ Make and check database connection.

        Assigns self.engine

        Args
        ------
        logger {class}: logger
        """

        if self.db_type == 'sql':
            connector = '{}{}:{}@{}:{}/{}'.format('mysql+mysqlconnector://',
                                                  self.user,
                                                  self.password,
                                                  self.host,
                                                  self.port,
                                                  self.db_name)
            if logger is not None:
                logger.info(" Using database connection "
                            "{}@{}:{}".format(self.user, self.host, self.port))

        if self.db_type == 'sqlite':
            connector = 'sqlite:///' + self.database
            if logger is not None:
                logger.info(" Using database connector: {}".format(connector))

        try:
            engine = create_engine(connector)
        except Exception as e:
            print(e)

        self.assign_vars({'engine': engine})

    def check_tables(self, logger=None):
        """ Make database tables.

        Args
        ------
        logger {class}: logger
        """

        if not hasattr(self, 'engine'):
            if logger is not None:
                logger.error(" Database connector required...")
            else:
                print('Database connector required')

        Base.metadata.create_all(self.engine)

    def insert(self, dbtable, kwargs, logger=None):
        """ Put data on database.

        Args
        ------
        dbtable {string}: string format of database table name (i.e., 'Pixels')
        kwargs {dict}: data to put on the database {field: value}
        logger {class}: logger
        """

        if type(kwargs) != dict:
            raise TypeError("kwargs must be a dict")

        # get database table
        tbl = sa.Table(dbtable, sa.MetaData(), autoload_with=self.engine)

        # table columns
        fields = list(tbl.columns.keys())

        for key in kwargs:
            if key not in fields:
                if logger is not None:
                    logger.warning(" kwargs field '{}' not in table columns "
                                   "{}".format(key, fields))
                else:
                    print("query() kwargs field '{}' not in table columns "
                          "{}".format(key, fields))

        with self.engine.connect() as dbcon:
            dbcon.execute(tbl.insert(), kwargs)

    def operators(self, operator):
        """ Translate sqlalchemy query operators for readability.

            eq for ==
            lt for <
            ge for >=
            in for in_
            like for like

        """

        if not isinstance(operator, str):
            raise TypeError("operators must be a str")

        omap = {'==': 'eq',
                '<': 'lt',
                '>': 'gt',
                '>=': 'ge',
                '<=': 'le',
                'in_': 'in',
                'like': 'like'}

        if operator not in list(omap.keys()):
            raise ValueError("operator must be in "
                             "'{}'".format(list(omap.keys())))

        return omap[operator]

    def query(self, params, columns=None, index=None, logger=None):
        """ Query database values.

        Columns and index are used to filter results after query. Columns
        is applied before index, so if you want index='index', 'index' must
        also appear in columns=['index'].

        Args
        ------
        params {dict}: query parameters, in the format
                        {table: (column, operator, value)}
        columns {list}: list of database table columns to filter results
        index {list}: name of database table column to set as returned
            DataFrame index
        logger {class}: logger

        Returns
        ------
        results {DataFrame}: query results
        """

        # tparams = {Pixels: {(Pixels.model_row, '==', int(pt[0]))}}

        unique_tables = set(params.keys())
        tables = list(params.keys())
        ntables = len(unique_tables)

        if ntables < 1 or ntables > 2:
            raise ValueError("params must have <= 2 unique tables")

        dbsession = sessionmaker(bind=self.engine)
        session = dbsession()

        from snowav.database.tables import Pixels, PixelsData

        n = 0
        for t in params.keys():
            print(t)
            if t == 'Pixels':
                table = Pixels
            if t == 'PixelsData':
                table = PixelsData

            for f in params[t]:
                print(f)
                raw = f

                try:
                    key, op, value = raw
                except ValueError:
                    raise Exception('Invalid filter: %s' % raw)

                column = getattr(table, key, None)

                if not column:
                    raise Exception('Invalid filter column: %s' % key)
                if op == 'in':
                    if isinstance(value, list):
                        filt = column.in_(value)
                    else:
                        filt = column.in_(value.split(','))
                else:
                    try:
                        attr = list(filter(
                            lambda e: hasattr(column, e % op),
                            ['%s', '%s_', '__%s__']
                        ))[0] % op
                    except IndexError:
                        raise Exception('Invalid filter operator: %s' % op)
                    if value == 'null':
                        value = None

                    if n == 0:
                        filt = getattr(column, attr)(value)
                        print("FIRST\n", filt)
                    else:
                        t = getattr(column, attr)(value)
                        filt = and_(filt, t)
                        print('AND\n ', filt)

                n += 1

        if ntables == 1:
            qry = session.query(table).filter(filt)
        else:
            qry = session.query(tables[0]).join(tables[1])

        # qry = qry.filter(filt)

        print('QRY\n', qry)
        results = pd.read_sql(qry.statement, qry.session.connection())
        session.close()

        print(results)
        print(x)

        """
        if index is None:
            index = []
        if columns is None:
            columns = []
        if not isinstance(params, dict):
            raise TypeError("params must be a dict")

        unique_tables = set(params.keys())
        ntables = len(unique_tables)

        if ntables < 1 or ntables > 2:
            raise ValueError("params can have no more than 2 unique tables")

        if not isinstance(columns, list):
            raise TypeError("columns must be a list")

        if not isinstance(index, list):
            raise TypeError("index must be a list")

        dbsession = sessionmaker(bind=self.engine)
        session = dbsession()
        tables = list(params.keys())
        ntables = len(tables)

        if ntables == 1:
            qry = session.query(tables[0]).filter_by(**params[tables[0]])
        else:
            qry = session.query(tables[0], tables[1])

        results = pd.read_sql(qry.statement, qry.session.connection())
        session.close()

        # apply columns
        applied_columns = []
        for c in columns:
            if c in results.columns:
                applied_columns.append(c)

            else:
                if logger is not None:
                    logger.warning(" columns='{}' not in {}, will not "
                                   "apply".format(c, list(results.columns)))
                else:
                    print("columns='{}' not in {}, will not "
                          "apply".format(c, list(results.columns)))

        if len(applied_columns) > 0:
            results = results[applied_columns]

        # apply index
        if len(index) > 0:
            if index[0] not in list(results.columns):
                if logger is not None:
                    logger.warning(" index='{}' not in columns={}, will not "
                                   "apply".format(index[0],
                                                  list(results.columns)))
                else:
                    print("index='{}' not in columns={}, will not "
                          "apply".format(index[0], list(results.columns)))
            else:
                results.set_index(index, inplace=True)
                results.sort_index(inplace=True)

        return results
        """

    def delete(self, dbtable, kwargs, logger=None):
        """ Delete database records.

        Args
        ------
        dbtable {string}: string format of database table name (i.e., 'Pixels')
        kwargs {dict}: items for deletion {field: value}
        logger {class}: logger
        """

        if not isinstance(kwargs, dict):
            raise TypeError("kwargs must be dict")

        dbsession = sessionmaker(bind=self.engine)
        session = dbsession()

        tbl = sa.Table(dbtable, sa.MetaData(), autoload_with=self.engine)
        results = session.query(tbl).filter_by(**kwargs)
        results.delete(synchronize_session=False)

        session.commit()
        session.close()

        if logger is not None:
            logger.debug(" Deleted database records")


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
    except:
        raise Exception('Failed to make database connection with '
                        '{}'.format(connector))
    try:
        Base.metadata.create_all(engine)
    except:
        raise Exception('Failed establishing Base for sqlalchemy')

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

    for k, v in values.items():
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

    value_options = ['swe_z', 'swe_vol', 'density', 'precip', 'precip_z', 'precip_vol',
                     'rain_z', 'rain_vol', 'swe_avail', 'swe_unavail', 'coldcont',
                     'swi_z', 'swi_vol', 'depth', 'snow_line', 'mean_air_temp',
                     'evap_z', 'L_v_E', 'melt', 'lwc', 'temp_surface',
                     'temp_lower', 'temp_bulk', 'depth_lower_layer', 'h20_sat',
                     'R_n', 'H', 'G', 'M', 'delta_Q']

    method_options = ['sum', 'difference', 'daily', 'end']

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
        df = pd.DataFrame(index=['date_time'], columns=plotorder)

    else:
        df = pd.DataFrame(np.nan, index=edges, columns=plotorder)

    for bid in plotorder:
        results = query(connector, start_date, end_date, run_name, basins, bid, value)

        if results.empty:
            raise Exception('No results on database for '
                            'run_name: {}, '
                            'start_date: {}, '
                            'end_date: {}, '
                            'value: {} '.format(run_name, start_date,
                                                end_date, value))

        if method == 'daily':
            e = results[(results['elevation'] == 'total')
                        & (results['date_time'] >= start_date)
                        & (results['date_time'] <= end_date)]
            if e.empty:
                raise Exception('No results on database for '
                                'run_name: {} ',
                                'start_date: {} ,'
                                'end_date: {}, '
                                'value: {} ,'
                                'elevation: {}'.format(run_name, start_date,
                                                       end_date, value, elev))
            else:
                e = e.set_index('date_time')
                e.sort_index(inplace=True)

                if bid == plotorder[0]:
                    df = e[['value']].copy()
                    df = df.rename(columns={'value': bid})

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
                                    'elev={}, {} to {}'.format(run_name, elev,
                                                               start_date, end_date))
                else:
                    if e['value'].values[0] is None:
                        df.loc[elev, bid] = np.nan
                    else:

                        # The database can occasionally get multiple values
                        # for the same record if it doesn't exit cleanly
                        if len(e['value'].values) == 1:
                            df.loc[elev, bid] = e['value'].values
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
                                    'elev={}, {} to {}'.format(run_name, elev,
                                                               start_date, end_date))
                else:
                    df.loc[elev, bid] = np.nansum(e['value'].values - s['value'].values)

        if method == 'sum':
            for elev in edges:
                e = results[(results['elevation'] == str(elev))
                            & (results['date_time'] >= start_date)
                            & (results['date_time'] <= end_date)]
                if e.empty:
                    raise Exception('Empty results for database query in '
                                    'database.collect(), for run_name={}, '
                                    'elev={}, {} to {}'.format(run_name, elev,
                                                               start_date, end_date))
                else:
                    df.loc[elev, bid] = e['value'].sum(skipna=False)

    df = df[plotorder]

    return df


def query(connector, start_date, end_date, run_name, basins, bid=None,
          value=None, rid=None):
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
    wval = Watershed(watershed_id=1, watershed_name=plotorder[0])
    session.add(wval)
    session.commit()

    # Initialize basins within the watershed
    for i, name in enumerate(plotorder):
        logger.append(' Adding basin {} to database'.format(name))

        bval = Basin(watershed_id=1,
                     basin_id=i + 1,
                     basin_name=name)

        session.add(bval)

    session.commit()
    session.close()

    logger.append(' Completed initializing basins and watersheds')

    return logger


def run_metadata(cfg):
    '''
    Create database RunMetadata for each snowav run.

    Args
    ------
        Variables: class, from snowav/database/models.py
        run_name : str

    '''

    watershed_id = int(cfg.basins[cfg.plotorder[0]]['watershed_id'])

    session = make_session(cfg.connector)
    qry = session.query(RunMetadata)
    df = pd.read_sql(qry.statement, qry.session.connection())
    session.close()

    # Increase each runid by 1
    if not df.run_id.values.size:
        cfg.run_id = 1
    else:
        r = df['run_id'].values
        cfg.run_id = next(filterfalse(set(r).__contains__, count(1)))

    # in most cases it's not practical to save all run_dirs
    trdir = ','.join(cfg.run_dirs)
    if len(trdir) > 1200:
        trdir = cfg.run_dirs[0]

    values = {'run_id': int(cfg.run_id),
              'run_name': cfg.run_name,
              'watershed_id': int(watershed_id),
              'pixel': int(cfg.pixel),
              'description': '',
              'smrf_version': 'smrf' + smrf_version,
              'awsm_version': 'aswm' + awsm_version,
              'snowav_version': 'snowav' + cfg.snowav_version,
              'data_type': '',
              'data_location': trdir,
              'file_type': '',
              'config_file': cfg.config_file,
              'proc_time': datetime.now()}

    insert(cfg.connector, 'RunMetadata', values)

    cfg.vid = {}
    for v in [*cfg.variables.variables.keys()]:
        u = cfg.variables.variables[v]['units']
        if (('vol' in v) or ('avail' in v)):
            u = cfg.units
        if (('z' in v) or (v == 'depth')):
            u = cfg.depthlbl
        if v == 'density':
            u = 'kg m^-3'
        if v == 'coldcont':
            u = 'MJ'

        variables = {'run_id': cfg.run_id,
                     'variable': v,
                     'unit': u,
                     'name': cfg.variables.variables[v]['description']}

        insert(cfg.connector, 'VariableUnits', variables)

        # After inserting into VariableUnits, get the existing id
        session = make_session(cfg.connector)
        qry = session.query(VariableUnits).filter(
            VariableUnits.run_id == cfg.run_id)
        session.close()
        df = pd.read_sql(qry.statement, qry.session.connection())
        cfg.vid[v] = df[df['variable'] == v]['id'].values[0]

    # snow_line
    variables = {'run_id': cfg.run_id,
                 'variable': 'snow_line',
                 'unit': cfg.elevlbl,
                 'name': 'snow_line'}

    insert(cfg.connector, 'VariableUnits', variables)
    session = make_session(cfg.connector)
    qry = session.query(VariableUnits).filter(
        VariableUnits.run_id == cfg.run_id)
    session.close()
    df = pd.read_sql(qry.statement, qry.session.connection())
    cfg.vid['snow_line'] = df['id'].values[0]


def connect(sqlite=None, sql=None, plotorder=None, user=None,
            password=None, host=None, port=None, convert=False,
            add=False):
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
                  'section'.format(user, password, host, port))
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
                          'does not exist, it is being created...')
            query = ("CREATE DATABASE {};".format(sql))
            cursor.execute(query)
            log = create_tables(db_engine, plotorder)

            engine = create_engine(db_engine)
            Base.metadata.create_all(engine)
            DBSession = sessionmaker(bind=engine)
            session = DBSession()

            for out in log:
                logger.append(out)

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
                logger.append(' Using database connection {}@{}:{} "{}" for '
                              'results'.format(user, host, port, sql))
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

        wval = Watershed(watershed_id=int(wid),
                         watershed_name=plotorder[0])
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
    watersheds = {plotorder[0]: {'watershed_id': wid, 'watershed_name': '',
                                 'basins': '', 'shapefile': ''}}

    # Initialize basins within the watershed
    for i, name in enumerate(plotorder):

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

            bval = Basin(watershed_id=int(wid),
                         basin_id=int(bid + i),
                         basin_name=name)

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

        basins[name] = {'watershed_id': wid, 'basin_id': bid}

    session.commit()
    session.close()

    return basins, connector, logger


def convert_watershed_names(name):
    '''
    Convert watershed from topo.nc mask to existing snowav database version.

    '''
    convert_list = ['Boise River Basin', 'Lakes Basin', 'Merced River Basin',
                    'Kaweah River Basin', 'Kings River Basin']

    if name in convert_list:
        watershed = name.split(' ')[0]

    elif name == 'San Joaquin River Basin':
        watershed = 'San Joaquin'

    else:
        watershed = name

    return watershed


def package(connector, basins, df, run_id, vid, output, dtime):
    """ Put process() results on the database.

    Args
    ------
    connector {str}: database connector
    basins {dict}: snowav basins identifier
    run_id {int}: run_id for database
    vid {int}: variable id
    df {DataFrame}: results dataframe
    output {str}: output variable (i.e. 'swe_z')
    dtime {datetime}: datetime
    """

    for basin in df:

        for iters, val in enumerate(df[basin].values):
            if np.isnan(val):
                val = None
            else:
                val = float(val)

            values = {'basin_id': int(basins[basin]['basin_id']),
                      'run_id': int(run_id),
                      'date_time': dtime,
                      'variable': output,
                      'variable_id': int(vid[output]),
                      'value': val,
                      'elevation': str(df[basin].index[iters])}

            insert(connector, 'Results', values)
