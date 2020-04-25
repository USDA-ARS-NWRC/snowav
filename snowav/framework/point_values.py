import coloredlogs
from datetime import datetime
import logging
import netCDF4 as nc
import numpy as np
import os
import pandas as pd
from sys import exit
import utm

from inicheck.config import MasterConfig
from inicheck.output import print_config_report
from inicheck.tools import get_user_config, check_config, cast_all_variables
import snowav
from snowav.database.models import AwsmInputsOutputs
from snowav.database.database import Database
from snowav.database.tables import Pixels, PixelsData


class PointValues(object):
    """ Point values processing to pull time-series model results from
    specific pixels.

    If config files aren't passed, the class contains only hard-coded string
    definitions.

    NOTE: Class attributes are automatically assigned from the master config
    file.

    Args
    ------
    master_config: {string} path to inicheck master config
    config_file: {sring} path to config file
    """

    def __init__(self, master_config=None, config_file=None):

        # run processing if both are passed, otherwise PointValues
        # contains only the hard-coded string definitions
        if sum([x is not None for x in [master_config, config_file]]) == 2:

            # inicheck the config file
            mcfg = MasterConfig(path=master_config)
            ucfg = get_user_config(config_file, mcfg=mcfg)
            ucfg.apply_recipes()
            ucfg = cast_all_variables(ucfg, ucfg.mcfg)

            # check config, warn and exit if necessary
            warnings, errors = check_config(ucfg)
            if len(errors) > 0 or len(warnings) > 0:
                print_config_report(warnings, errors)
                if len(errors) > 0:
                    exit()

            # get a list of master config fields to check assignment of fields
            # in PointValues class
            mcfg_fields = []
            for section, d in ucfg.mcfg.cfg.items():
                mcfg_fields += list(d.keys())

            cfg = ucfg.cfg

            # Add in class attributes from config file
            for section, d in cfg.items():
                for item, v in d.items():

                    # assign if it exists in the master config
                    if item in mcfg_fields:
                        setattr(self, item, v)

            self.create_log()
            self.logger.info(' Starting point values processing...')

            # additional checks for user config options below...

            # check sqlite base path and file extension
            if self.database is not None:
                if not os.path.exists(os.path.dirname(self.database)):
                    self.logger.error(" [point_values] database: '{}' base "
                                      "path not valid".format(self.database))
                    raise OSError("invalid path -> {}".format(self.database))

                if os.path.splitext(self.database)[1] != '.db':
                    self.logger.error(" [point_values] database: '{}' "
                                      "must be '.db' "
                                      "file".format(self.database))
                    raise OSError("file type must be '.db'")
                self.db_type = 'sqlite'

            # if sql database, make sure all credentials exist
            if self.database is None:
                input_list = [self.user, self.password, self.host, self.port]
                if sum([x is not None for x in input_list]) != 4:
                    self.logger.warning(" User has not supplied all database "
                                        "credentials in [point_values]")
                self.db_type = 'sql'

            # csv saving requires directory
            if self.csv_output and self.csv_output_dir is None:
                self.logger.error(" [point_values] csv_output: True, must "
                                  "supply csv_output_dir field")
                exit()

            self.start_date_str = self.start_date.date().strftime("%Y%m%d")
            self.end_date_str = self.end_date.date().strftime("%Y%m%d")
            self.begin_time = datetime.now()

            snowav_properties = AwsmInputsOutputs()
            all_properties = (snowav_properties.smrf_variables +
                              snowav_properties.awsm_variables)

            # check and assign properties
            self.properties_lookup = {}

            if self.properties == ['all']:
                self.properties = all_properties
                if 'precip_z' in self.properties:
                    self.properties.remove('precip_z')
            else:
                for p in self.properties:
                    if p not in all_properties:
                        self.logger.error(' Provided [basin] properties: '
                                          '"{}" not a valid property, must be '
                                          'in {}'.format(p, all_properties))

            self.bandsmap = {}
            for p in all_properties:
                if p != 'precip_z':
                    self.bandsmap[p] = snowav_properties.vars[p]['nc_name']

            # make a lookup dict that is {file: [variables]} so we know what to
            # load for each file
            files = []
            var_list = []
            for p in snowav_properties.vars.keys():
                file = snowav_properties.vars[p]['nc_name']
                if (file is not None and
                        snowav_properties.vars[p]['file'] not in files):
                    files.append(snowav_properties.vars[p]['file'])

            for file in files:
                if file not in self.properties_lookup.keys():
                    self.properties_lookup[file] = {}
                    self.properties_lookup[file]['variables'] = []
                    self.properties_lookup[file]['bands'] = {}

            self.unit_map = {}
            for f in self.properties_lookup.keys():
                for v in self.properties:
                    variable = self.bandsmap[v]

                    self.unit_map[variable] = \
                        snowav_properties.vars[v]['units']

                    if f == snowav_properties.vars[v]['file']:
                        self.properties_lookup[f]['variables'].append(v)
                        band = snowav_properties.vars[v]['band']
                        self.properties_lookup[f]['bands'][v] = band

                        # make list for the data dfs - these are the headings

                        # The variable names snow.nc and snow_density.nc are
                        # both 'snow_density', this is an issue!
                        # This is also addressed in load_data() and save_csv()
                        if v == 'density':
                            var_list.append('density')
                        else:
                            var_list.append(variable)

            self.logger.debug(" Variables: {}".format(var_list))

            # these will be converted to dicts with location keys later
            self.tvar = pd.DataFrame(columns=var_list)

        # hard coded strings relating to AWSM and PointValues csv file
        # NOTE that database Pixels.location field is created by:
        #    os.path.abspath(os.path.join(pv.run_dirs[0], '..', '..'))
        # in load_data()
        self.smrf_dir = 'smrfOutputs'
        self.data_str = '/data'
        self.datadata_str = '/data/data'
        self.runs_str = '/runs'
        self.runsrun_str = '/runs/run'
        self.name_col = 'name'
        self.description_col = 'description'
        self.model_xind = 'model xind'
        self.model_yind = 'model yind'
        self.index_columns = [self.model_xind, self.model_yind]
        self.location_columns = ['latitude', 'longitude']
        self.csv_cols = [self.name_col,
                         self.description_col,
                         self.model_xind,
                         self.model_yind] + self.location_columns

    def assign_vars(self, kwargs):
        """ Assign PointValues attributes.

        Args
        ------
        kwargs {dict}: {item, value}
        """

        if not isinstance(kwargs, dict):
            raise TypeError("kwargs must be dict")

        for item, v in kwargs.items():
            setattr(self, item, v)

    def create_log(self):
        """ Create logger. """

        level_styles = {'info': {'color': 'white'},
                        'error': {'color': 'red'},
                        'debug': {'color': 'green'},
                        'warning': {'color': 'yellow'}}

        field_styles = {'hostname': {'color': 'magenta'},
                        'programname': {'color': 'cyan'},
                        'name': {'color': 'white'},
                        'levelname': {'color': 'white', 'bold': True},
                        'asctime': {'color': 'green'}}

        # start logging
        log_level = self.log_level.upper()

        numeric_level = getattr(logging, log_level, None)
        if not isinstance(numeric_level, int):
            raise ValueError('Invalid log level: %s' % log_level)

        fmt = '%(levelname)s:%(name)s:%(message)s'
        logging.basicConfig(level=numeric_level)
        coloredlogs.install(level=numeric_level,
                            fmt=fmt,
                            level_styles=level_styles,
                            field_styles=field_styles)

        self.assign_vars({'logger': logging.getLogger(__name__)})


def get_paths(pv):
    """ Get the 'data' and 'runs' directories for inputs and outputs.

        ../runs/runYYYYMMDD
        ../data/dataYYYYMMDD/smrfOutputs

    Assigns:
    pv.data_dirs
    pv.runs_dirs

    Args
    ------
    pv {class} PointValues class
    """

    # get all run dirs
    run_path = pv.directory
    r = os.listdir(run_path)
    run_dirs = [os.path.join(run_path, x) for x in r]
    run_dirs.sort()

    for dd in run_dirs:
        pv.logger.debug(" Checking run dir: {}".format(dd))

    # get all data dirs
    data_path = run_path.replace(pv.runs_str, pv.data_str)
    d = os.listdir(data_path)
    d_dirs = [os.path.join(x, pv.smrf_dir) for x in d]
    data_dirs = [os.path.join(data_path, x) for x in d_dirs]
    data_dirs.sort()

    for dd in run_dirs:
        pv.logger.debug(" Checking data dir: {}".format(dd))

    if len(run_dirs) != len(data_dirs):
        pv.logger.error(' There are mismatched number of subdirectories '
                        'in {} and {}'.format(run_path, data_path))

    # filter do the ones in requested date range
    if sum([x is not None for x in [pv.start_date, pv.end_date]]) == 2:
        pv.logger.debug(' Parsing folder dates...')
        pv.data_dirs = []
        pv.run_dirs = []

        for d in data_dirs:
            s = d.split(pv.datadata_str)[1].split('/' + pv.smrf_dir)[0]
            dt = datetime(int(s[:4]), int(s[4:6]), int(s[-2::]))

            if pv.start_date <= dt <= pv.end_date:
                pv.data_dirs.append(d)
                pv.logger.info(' Using data dir: {}'.format(d))

        for r in run_dirs:
            s = r.split(pv.runsrun_str)[1]
            dt = datetime(int(s[:4]), int(s[4:6]), int(s[-2::]))
            if pv.start_date <= dt <= pv.end_date:
                pv.run_dirs.append(r)
                pv.logger.info(' Using run dir: {}'.format(r))

        # check again in case the parsing changed anything
        if len(run_dirs) != len(data_dirs):
            pv.logger.error(' There are mismatched number of subdirectories '
                            'in {} and {}'.format(run_path, data_path))


def get_basin_info(pv):
    """ Load topo.nc file and get basin info.

    Args
    ------
    pv {class} PointValues class
    """

    pv.logger.debug(" Loading basin info from {}".format(pv.topo))

    # load topo.nc file
    try:
        data = nc.Dataset(pv.topo)
    except Exception as e:
        print(e)

    basin_name = data.variables['mask'].long_name
    pv.basin_name = ''.join(basin_name.split())
    pv.basin_x = len(data.variables['x'])
    pv.basin_y = len(data.variables['y'])
    pv.utm_x = data.variables['x'][:]
    pv.utm_y = data.variables['y'][:]
    pv.dem = data.variables['dem'][:]
    data.close()

    pv.logger.debug(" Basin y: {}".format(pv.basin_y))
    pv.logger.debug(" Basin x: {}".format(pv.basin_x))


def get_locations(pv):
    """ Get location data from csv file. If using latitude and longitude,
    convert to model x and y indices.

    Assigns:
    pv.data

    Args
    ------
    pv {class} PointValues class
    """

    pv.logger.debug(" Loading data for locations from "
                    "{}".format(pv.locations_csv))

    # load csv
    try:
        data = pd.read_csv(pv.locations_csv)
    except Exception as e:
        print(e)

    if data.empty:
        pv.logger.error(" {} is empty ".format(pv.locations_csv))

    columns = list(data.columns)

    if pv.name_col not in columns:
        pv.logger.warning(" '{}' column is missing".format(pv.name_col))

    if pv.description_col not in columns:
        pv.logger.warning(" '{}' column is missing".format(pv.description_col))

    # check if df contains complete set of either headings
    check = True

    # index
    if all(e in columns for e in pv.index_columns):
        pv.logger.debug(" Checking file {} "
                        "headings".format(pv.index_columns))

        if np.isnan(data[pv.index_columns].values).any() > 0:
            pv.logger.debug(" {} may contain empty "
                            "fields".format(pv.index_columns))
        else:
            data = data[pv.index_columns + [pv.name_col] + [pv.description_col]]
            check = False
            pv.logger.debug(" Using data in "
                            "{}".format(pv.index_columns))

            # check that indices are within basin dimensions
            if not all(x >= 0 for x in list(data[pv.model_xind].values)):
                pv.logger.error(" csv file contains index value <0")
                raise ValueError("Index out of bounds")

            if not all(x <= pv.basin_x for x in
                       list(data[pv.model_xind].values)):
                pv.logger.error(" csv file contains index value "
                                ">{}".format(pv.basin_x))
                raise ValueError("Index out of bounds")

            if not all(x >= 0 for x in list(data[pv.model_yind].values)):
                pv.logger.error(" csv file contains index value <0")
                raise ValueError("Index out of bounds")

            if not all(x <= pv.basin_y for x in
                       list(data[pv.model_yind].values)):
                pv.logger.error(" csv file contains index value "
                                ">{}".format(pv.basin_y))
                raise ValueError("Index out of bounds")

    # latitude and longitude
    if all(e in columns for e in pv.location_columns) and check:
        pv.logger.debug(" Checking file {} "
                        "headings".format(pv.location_columns))

        if np.isnan(data[pv.location_columns].values).any() > 0:
            pv.logger.warning(" {} may contain empty "
                              "fields".format(pv.location_columns))
        else:
            check = False
            pv.logger.debug(" Using data in "
                            "{}".format(pv.location_columns))

        # convert to indices and assign
        for i, (index, row) in enumerate(data.iterrows()):
            ll = utm.from_latlon(row['latitude'], row['longitude'])
            xind = np.where(abs(pv.utm_x - ll[0]) ==
                            min(abs(pv.utm_x - ll[0])))[0]
            yind = np.where(abs(pv.utm_y - ll[1]) ==
                            min(abs(pv.utm_y - ll[1])))[0]

            if xind == 0 or xind == pv.basin_x - 1:
                pv.logger.warning(" For {} xind = {}, "
                                  "may be out of basin".format(ll, xind))

            if yind == 0 or yind == pv.basin_x - 1:
                pv.logger.warning(" For {} yind = {}, "
                                  "may be out of basin".format(ll, yind))

            # assign x, y, name, and description
            data.loc[i, pv.model_xind] = int(xind)
            data.loc[i, pv.model_yind] = int(yind)

    if check:
        pv.logger.error(" csv must contain set of either {} or {} "
                        "columns".format(pv.index_columns,
                                         pv.location_columns))

    pv.data = data.copy()


def save_csv(pv):
    """ Save results in csv format.

    Args
    ------
    pv {class} PointValues class
    """

    now_str = datetime.now().date().strftime("%Y-%m-%d")
    date_col = 'Date generated: {}'.format(now_str)

    header = ['USDA Agicultural Research Service AWSM Model Data',
              'Data provided are SMRF inputs (hourly) and iSnobal outputs (daily)',
              date_col,
              'Contact: Scott Havens <scott.havens@usda.gov>\n']

    for loc in pv.var_dict.keys():
        df = pv.var_dict[loc]['data']
        df.sort_index(inplace=True)
        df.index.name = 'date_time'

        # add units to the headings
        columns = list(df.columns)

        for c in columns:
            # prevent renaming to snow_density
            if c == 'density':
                rdict = {c: "density [kg/m^3]"}
            else:
                rdict = {c: "{} [{}]".format(c, pv.unit_map[c])}
            df = df.rename(columns=rdict)

        loc_str = "x{}_y{}".format(str(int(loc[1])), str(int(loc[0])))
        name = "{}_{}_{}_variables_{}.csv".format(pv.basin_name,
                                                  pv.start_date_str,
                                                  pv.end_date_str,
                                                  loc_str)

        # df = df.sort_index(inplace=True)

        filename = os.path.join(pv.csv_output_dir, name)

        with open(filename, mode='w', encoding='utf-8') as f:
            f.write('\n'.join(header))

        df.to_csv(filename, encoding='utf-8', mode='a')

        pv.logger.info(' Saved: {}'.format(filename))


def load_data(pv):
    """ Load smrf and iSnobal outputs from pv.run_dirs and pv.data_dirs.

    Assigns:
    pv.var_dict {dict}: {(x, y): 'data': df,
                                 'name': name,
                                 'description': description,
                                 'location': location}

    where df.columns = ['air_temp','vapor_pressure',...]

    Args
    ------
    pv {class}: PointValues class
    """

    base_location = os.path.abspath(os.path.join(pv.run_dirs[0], '..', '..'))

    pv.var_dict = {}
    for index, row in pv.data.iterrows():
        loc = (int(row[pv.model_xind]), int(row[pv.model_yind]))
        pv.var_dict[loc] = {'data': pv.tvar.copy(),
                            'name': row[pv.name_col],
                            'description': row[pv.description_col],
                            'location': base_location}

    # each runs/ and data/ directory
    for n, d in enumerate(pv.run_dirs + pv.data_dirs):
        pv.logger.debug(" Working in {}".format(d))

        # each file, i.e. snow.nc
        for file in os.listdir(d):
            filepath = os.path.join(d, file)

            # if it's a valid file and one we need
            if os.path.isfile(filepath) and \
                    file in list(pv.properties_lookup.keys()):

                try:
                    data = nc.Dataset(filepath)
                    date_time = nc.num2date(data.variables['time'][:],
                                            data.variables['time'].units)
                except Exception as e:
                    print(e)

                pv.logger.debug(" Opened {}".format(file))

                # ensure that outputs are hourly
                for i, dt in enumerate(date_time):
                    date_time[i] = date_time[i].replace(minute=0,
                                                        second=0,
                                                        microsecond=0)

                # for each band
                for v, b in pv.properties_lookup[file]['bands'].items():

                    # for each location
                    for item, row in pv.data.iterrows():
                        loc = (int(row[pv.model_xind]), int(row[pv.model_yind]))

                        # snowav --> .nc
                        if v in pv.bandsmap.keys() and b is not None:
                            v = pv.bandsmap[v]

                        try:
                            # nans in some of the netcdf files give runtime
                            # warnings when being loaded. As far as I can tell
                            # the only fix is fill_value when creating.
                            img = data.variables[v][:]
                        except KeyError as e:
                            print('{} not a variable in {}'.format(e, file))

                        # check dimensions of the first image opened
                        if n == 0 and (pv.basin_y, pv.basin_x) != img[0].shape:
                            pv.logger.error(" output dimensions {} "
                                            "do not match topo "
                                            "({}, {})".format(img[0].shape,
                                                              pv.basin_y,
                                                              pv.basin_x))

                        for ih in range(0, len(date_time)):
                            if v == 'snow_density' and file == 'snow.nc':
                                vt = 'density'
                            else:
                                vt = v

                            if ih < 24:
                                dth = date_time[ih]
                                pv.var_dict[loc]['data'].loc[dth, vt] = \
                                    img[ih, loc[1], loc[0]]

                                pv.logger.debug(" Write to df: {}, {}, "
                                                "{}".format(v, loc, dth))

                data.close()

    pv.var_dict[loc]['data'].sort_index(inplace=True)


def put_on_database(db, pv):
    """ Database stream for checking for existing results, metdata, and
    results.

    Args
    ------
    db {class}: Database class
    pv {class}: PointValues class
    """

    # for each location
    for pt in pv.var_dict:

        df = pv.var_dict[pt]['data']
        start_date = min(df.index)
        # end_date = max(df.index)

        metadata = {
            'location': str(pv.var_dict[pt]['location']),
            'model_row': int(pt[0]),
            'model_col': int(pt[1]),
            'description': str(pv.var_dict[pt]['description']),
            'name': str(pv.var_dict[pt]['name']),
            'utm_x': float(pv.utm_x[pt[0]]),
            'utm_y': float(pv.utm_y[pt[1]]),
            'elevation': float(pv.dem[pt[1], pt[0]])
        }

        # currently query() is only grabbing keys from params and pixelsparams
        params = {Pixels: {'model_row': int(pt[0]),
                           'model_col': int(pt[1]),
                           'name': str(pv.var_dict[pt]['name']),
                           'location': str(pv.var_dict[pt]['location']),
                           'description': str(pv.var_dict[pt]['description'])
                           },
                  PixelsData: {'date_time': start_date}
                  }
        pixelsparams = {Pixels: {'model_row': int(pt[0]),
                                 'model_col': int(pt[1]),
                                 'name': str(pv.var_dict[pt]['name']),
                                 'location': str(pv.var_dict[pt]['location']),
                                 'description': str(pv.var_dict[pt]['description'])
                                 }
                        }

        # from sqlalchemy import create_engine, and_
        # from sqlalchemy.orm import sessionmaker
        #
        # tparams = {Pixels: {(Pixels.model_row, '==', int(pt[0]))}}
        #
        # dbsession = sessionmaker(bind=db.engine)
        # session = dbsession()
        # tables = list(tparams.keys())
        # raw = ('model_row', 'ge', 1)
        # try:
        #     key, op, value = raw
        # except ValueError:
        #     raise Exception('Invalid filter: %s' % raw)
        #
        # column = getattr(Pixels, key, None)
        #
        # if not column:
        #     raise Exception('Invalid filter column: %s' % key)
        # if op == 'in':
        #     if isinstance(value, list):
        #         filt = column.in_(value)
        #     else:
        #         filt = column.in_(value.split(','))
        # else:
        #     try:
        #         attr = list(filter(
        #             lambda e: hasattr(column, e % op),
        #             ['%s', '%s_', '__%s__']
        #         ))[0] % op
        #     except IndexError:
        #         raise Exception('Invalid filter operator: %s' % op)
        #     if value == 'null':
        #         value = None
        #     filt = getattr(column, attr)(value)
        #
        # qry = session.query(Pixels, PixelsData)
        # qry = qry.filter(filt)
        #
        # results = pd.read_sql(qry.statement, qry.session.connection())
        # session.close()
        #
        # print(qry)
        # print(results)
        #
        # print(x)

        params = {'Pixels': [
                            ('model_row', 'eq', int(pt[0])),
                            ('name', 'eq', 'MAM')
        ],
                  # 'PixelsData': [('date_time', 'ge', datetime(2020, 1, 1))]
                  }

        # first, see if records already exist with the same values
        results = db.query(params, logger=pv.logger)

        # filter outside of query for now
        results = results[(results.model_row == pt[0]) &
                          (results.model_col == pt[1]) &
                          (results.name == str(pv.var_dict[pt]['name'])) &
                          (results.date_time >= pv.start_date) &
                          (results.date_time <= pv.end_date)]

        # if no existing records, put on database
        if results.empty:

            # insert metadata
            db.insert('Pixels', metadata, logger=pv.logger)
            pv.logger.debug(" Metadata to database for point: "
                            "({},{})".format(metadata['model_row'],
                                             metadata['model_col']))

            # get associated Pixels.id metadata for PixelsData
            results = db.query(pixelsparams, logger=pv.logger)

            # for this case we can't filter by 'date_time', because there
            # aren't any records yet
            results = results[(results.model_row == pt[0]) &
                              (results.model_col == pt[1]) &
                              (results.name == str(pv.var_dict[pt]['name']))]

            max_record_id = int(max(pd.unique(results['id'])))

            for idx, row in df.iterrows():
                row = row.astype('float')
                p = {'pixel_id': max_record_id,
                     'date_time': idx}
                res = row.to_dict()
                res = {k: res[k] for k in res if not np.isnan(res[k])}
                pixeldata = {**p, **res}

                # insert data
                db.insert('PixelsData', pixeldata, logger=pv.logger)

            # log just one
            pv.logger.debug(" Data to database for point: "
                            "({},{})".format(metadata['model_row'],
                                             metadata['model_col']))

        # if there are existing record, check overwrite
        else:
            existing_record_id = pd.unique(results['id'])

            if pv.overwrite:
                # first, delete the existing metadata and records based on the
                # existing Pixels.id
                print(existing_record_id)
                for rec in existing_record_id:

                    # first, delete the data
                    dparams = {'pixel_id': int(rec)}
                    db.delete('PixelsData', dparams, logger=pv.logger)

                    results = db.query(params, logger=pv.logger)
                    results = results[(results.model_row == pt[0]) &
                                      (results.model_col == pt[1]) &
                                      (results.name == str(pv.var_dict[pt]['name'])) &
                                      (results.date_time >= pv.start_date) &
                                      (results.date_time <= pv.end_date)]

                    # if that removed everything, also delete the metadata
                    if results.empty:
                        db.delete('Pixels',
                                  {'id': int(rec)},
                                  logger=pv.logger)

                # insert metadata
                db.insert('Pixels', metadata, logger=pv.logger)
                pv.logger.debug(" Metadata to database for point: "
                                "({},{})".format(metadata['model_row'],
                                                 metadata['model_col']))

                # get associated Pixels.id metadata for PixelsData with the
                # same query params
                results = db.query(pixelsparams, logger=pv.logger)

                # filter outside of query for now
                results = results[(results.model_row == pt[0]) &
                                  (results.model_col == pt[1]) &
                                  (results.name == str(pv.var_dict[pt]['name']))]
                record_id = int(max(pd.unique(results['id'])))

                # prepare data to put on database
                df = pv.var_dict[pt]['data']

                for idx, row in df.iterrows():
                    row = row.astype('float')
                    p = {'pixel_id': int(record_id),
                         'date_time': idx}
                    res = row.to_dict()
                    res = {k: res[k] for k in res if not np.isnan(res[k])}
                    pixeldata = {**p, **res}

                    # insert data
                    db.insert('PixelsData', pixeldata, logger=pv.logger)

                # log just one
                pv.logger.debug(" Data to database for point: "
                                "({},{})".format(metadata['model_row'],
                                                 metadata['model_col']))

            else:
                pv.logger.info(" Database records exist for {}, ({},{}) "
                               "and [point_values] overwrite: False, "
                               "values not put on "
                               "database".format(metadata['location'],
                                                 metadata['model_row'],
                                                 metadata['model_col']))


def run_point_values(point_values_config, master_config, blank):
    """ Run point values processing.

    This does a bulk records delete for matching 'location', 'model_row',
    and 'model_col' on the database if existing records are present.

    Args
    ------
    point_values_config {string}: config file
    master_config {string}: point values master config
    blank {bool}: make a blank locations_csv file
    """

    # use default CoreConfig.ini if one is not supplied
    if master_config is None:
        master_config = os.path.abspath(
            os.path.join(snowav.__path__[0],
                         "./config/PointValuesCoreConfig.ini")
        )

    if point_values_config is not None:
        if not os.path.isfile(os.path.abspath(point_values_config)):
            raise Exception('{} not a valid file'.format(master_config))

        if not os.path.isfile(os.path.abspath(master_config)):
            raise Exception('{} not a valid file'.format(master_config))

        # initialize class
        pv = PointValues(master_config, point_values_config)

        # initialize Database class
        db = Database(user=pv.user, password=pv.password, host=pv.host,
                      port=pv.port, db_type=pv.db_type, database=pv.database)

        # make and check database connection
        db.make_connection(logger=pv.logger)

        # check database tables
        db.check_tables(logger=pv.logger)

        # load topo.nc and get basin info
        get_basin_info(pv)

        # get and parse locations from the csv file
        get_locations(pv)

        # get data paths
        get_paths(pv)

        # load data
        load_data(pv)

        # database stream
        put_on_database(db, pv)

        # save to csv
        if pv.csv_output:
            save_csv(pv)

        elapsed = str(datetime.now() - pv.begin_time)
        pv.logger.info(" Completed point values processing, elapsed "
                       "time: {}".format(elapsed))

    # make a blank locations_csv file
    if blank:
        pv = PointValues()
        blank = pd.DataFrame(columns=pv.csv_cols)
        blank.to_csv('basin_locations.csv', index=False)
        print('\nCreated a blank locations_csv file:\n '
              '    basin_locations.csv\n\n'
              'Fill in either ["model xind", "model yind"] '
              'or ["latitude","longitude"] columns\n')
