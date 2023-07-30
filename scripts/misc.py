"""
Miscellaneous model inputs.
"""
import os
import configparser
import glob
import pandas as pd
import geopandas as gpd
# from shapely.ops import transform
from shapely.geometry import shape, Point, mapping, LineString, MultiPolygon

CONFIG = configparser.ConfigParser()
CONFIG.read(os.path.join(os.path.dirname(__file__), 'script_config.ini'))
BASE_PATH = CONFIG['file_locations']['base_path']

DATA_RAW = os.path.join(BASE_PATH, 'raw', '..', '..', '..', 'data_raw')
DATA_PROCESSED = os.path.join(BASE_PATH, 'processed')


technologies = [
    'GSM',
    'UMTS',
    'LTE',
    'NR',
]


def get_countries():
    """
    Get all countries.
    """
    filename = "countries.csv"
    path = os.path.join(DATA_RAW, filename)
    
    countries = pd.read_csv(path, encoding='latin-1')
    countries = countries[countries.Exclude == 0]

    return countries


def get_regions(country, region_type):
    """
    Get region information.
    """
    if region_type == 'use_csv':
        filename = 'regions_{}_{}.shp'.format(
            country['lowest'],
            country['iso3']
        )
        folder = os.path.join(DATA_PROCESSED, country['iso3'], 'regions')
    elif region_type in [1, 2]:
        filename = 'regions_{}_{}.shp'.format(
            region_type,
            country['iso3']
        )
        folder = os.path.join(DATA_PROCESSED, country['iso3'], 'regions')
    elif region_type == 0:
        filename = 'national_outline.shp'
        folder = os.path.join(DATA_PROCESSED, country['iso3'])
    else:
        print("Did not recognize region_type arg provided to get_regions()")

    path = os.path.join(folder, filename)

    if not os.path.exists(path):
        print('This path did not exist/load: {}'.format(path))
        return []

    regions = gpd.read_file(path, crs='epsg:4326')#[:1]

    return regions


def process_country_shapes(iso3):
    """
    Creates a single national boundary for the desired country.
    Parameters
    ----------
    country : dict
        Contains all desired country information.
    """
    path = os.path.join(DATA_PROCESSED, iso3)

    if os.path.exists(os.path.join(path, 'national_outline.shp')):
        return 'Completed national outline processing'

    print('Processing country shapes')

    if not os.path.exists(path):
        os.makedirs(path)

    shape_path = os.path.join(path, 'national_outline.shp')

    path = os.path.join(DATA_RAW, 'gadm36_levels_shp', 'gadm36_0.shp')

    countries = gpd.read_file(path)

    single_country = countries[countries.GID_0 == iso3].reset_index()

    single_country = single_country.copy()
    single_country["geometry"] = single_country.geometry.simplify(
        tolerance=0.01, preserve_topology=True)

    single_country['geometry'] = single_country.apply(
        remove_small_shapes, axis=1)

    glob_info_path = os.path.join(DATA_RAW, 'countries.csv')
    load_glob_info = pd.read_csv(glob_info_path, encoding = "ISO-8859-1",
        keep_default_na=False)
    single_country = single_country.merge(
        load_glob_info, left_on='GID_0', right_on='iso3')

    single_country.to_file(shape_path, driver='ESRI Shapefile')

    return


def remove_small_shapes(x):
    """
    Remove small multipolygon shapes.
    Parameters
    ---------
    x : polygon
        Feature to simplify.
    Returns
    -------
    MultiPolygon : MultiPolygon
        Shapely MultiPolygon geometry without tiny shapes.
    """
    if x.geometry.type == 'Polygon':
        return x.geometry

    elif x.geometry.type == 'MultiPolygon':

        area1 = 0.01
        area2 = 50

        if x.geometry.area < area1:
            return x.geometry

        if x['GID_0'] in ['CHL','IDN']:
            threshold = 0.01
        elif x['GID_0'] in ['RUS','GRL','CAN','USA']:
            threshold = 0.01
        elif x.geometry.area > area2:
            threshold = 0.1
        else:
            threshold = 0.001

        new_geom = []
        for y in list(x['geometry'].geoms):
            if y.area > threshold:
                new_geom.append(y)

        return MultiPolygon(new_geom)


def process_regions(iso3, level):
    """
    Function for processing the lowest desired subnational
    regions for the chosen country.
    Parameters
    ----------
    country : dict
        Contains all desired country information.
    """
    regions = []

    for regional_level in range(1, int(level) + 1):

        filename = 'regions_{}_{}.shp'.format(regional_level, iso3)
        folder = os.path.join(DATA_PROCESSED, iso3, 'regions')
        path_processed = os.path.join(folder, filename)

        if os.path.exists(path_processed):
            continue

        print('Processing GID_{} region shapes'.format(regional_level))

        if not os.path.exists(folder):
            os.mkdir(folder)

        filename = 'gadm36_{}.shp'.format(regional_level)
        path_regions = os.path.join(DATA_RAW, 'gadm36_levels_shp', filename)
        regions = gpd.read_file(path_regions)

        regions = regions[regions.GID_0 == iso3]

        regions = regions.copy()
        regions["geometry"] = regions.geometry.simplify(
            tolerance=0.005, preserve_topology=True)

        regions['geometry'] = regions.apply(remove_small_shapes, axis=1)

        try:
            regions.to_file(path_processed, driver='ESRI Shapefile')
        except:
            print('Unable to write {}'.format(filename))
            pass

    return


income_lut = {
    'BFA': {
        ".6_1": 1773,
        ".2_.6": 637,
        "-.2_.2": 398,
        "-.6_-.2": 281,
        "-1_-.6": 180,
    },
    'MLI': {
        ".6_1": 1563,
        ".2_.6": 787,
        "-.2_.2": 545,
        "-.6_-.2": 399,
        "-1_-.6": 267,
    },
    'NER': {
        ".6_1": 1266,
        ".2_.6": 565,
        "-.2_.2": 415,
        "-.6_-.2": 316,
        "-1_-.6": 210,
    },
}


parameters = {
    'restoration_days': 21,
    'dependence_percent': 20,
    'annual_protection_costs': 800, 
    'time_period': 5,
    'total_sites_BFA': 1700, #towerxchange
    'total_sites_MLI': 1900, #towerxchange
    'total_sites_NER': 1800, #towerxchange
    'sites_to_protect': 5,
}


if __name__ == '__main__':

    countries = get_countries()
    for idx, country in countries.iterrows():
       print(country['country'])

