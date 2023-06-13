"""
Preprocess sites data.

Ed Oughton

February 2022

"""
import sys
import os
import configparser
import json
import pandas as pd
import geopandas as gpd
import pyproj
from shapely.ops import transform
from shapely.geometry import Point, box
import rasterio
from rasterio.mask import mask
from tqdm import tqdm

from misc import get_countries, process_country_shapes, process_regions, get_regions 

CONFIG = configparser.ConfigParser()
CONFIG.read(os.path.join(os.path.dirname(__file__),'..', 'scripts', 'script_config.ini'))
BASE_PATH = CONFIG['file_locations']['base_path']

DATA_RAW = os.path.join(BASE_PATH, 'raw')
DATA_PROCESSED = os.path.join(BASE_PATH, 'processed')


def run_preprocessing(iso3):
    """
    Meta function for running preprocessing.

    """
    filename = "countries.csv"
    path = os.path.join(DATA_RAW, filename)

    countries = pd.read_csv(path, encoding='latin-1')
    country = countries[countries.iso3 == iso3]
    country = country.to_records('dicts')[0]
    regional_level = int(country['gid_region'])

    print('Working on create_national_sites_csv')
    create_national_sites_csv(country)

    print('Working on process_country_shapes')
    process_country_shapes(iso3)

    print('Working on process_regions')
    process_regions(iso3, regional_level)

    print('Working on create_national_sites_shp')
    create_national_sites_shp(iso3)

    print('Working on process_acled_layer')
    process_acled_layer(iso3)

    print('Working on subset_acled_telecom')
    intersect_acled(iso3)

    print('Working on process_scdi_layer')
    process_scdi_layer(iso3)
    
    return


def create_national_sites_csv(country):
    """
    Create a national sites csv layer for a selected country.

    """
    iso3 = country['iso3']#.values[0]

    filename = "mobile_codes.csv"
    path = os.path.join(DATA_RAW, filename)
    mobile_codes = pd.read_csv(path)
    mobile_codes = mobile_codes[['iso3', 'mcc', 'mnc']].drop_duplicates()
    all_mobile_codes = mobile_codes[mobile_codes['iso3'] == iso3]
    all_mobile_codes = all_mobile_codes.to_dict('records')

    output = []

    filename = '{}.csv'.format(iso3)
    folder = os.path.join(DATA_PROCESSED, iso3, 'sites')
    path_csv = os.path.join(folder, filename)

    ### Produce national sites data layers
    if os.path.exists(path_csv):
        return

    print('-site.csv data does not exist')
    print('-Subsetting site data for {}'.format(iso3))

    if not os.path.exists(folder):
        os.makedirs(folder)

    filename = "cell_towers_2022-12-24.csv"
    path = os.path.join(DATA_RAW, filename)

    for row in all_mobile_codes:

        # if not row['mnc'] in [10,2,11,33,34,20,94,30,31,32,27,15,91,89]:
        #     continue

        mcc = row['mcc']
        seen = set()
        chunksize = 10 ** 6
        for idx, chunk in enumerate(pd.read_csv(path, chunksize=chunksize)):

            country_data = chunk.loc[chunk['mcc'] == mcc]#[:1]

            country_data = country_data.to_dict('records')

            for site in country_data:

                # if not -4 > site['lon'] > -6:
                #     continue

                # if not 49.8 < site['lat'] < 52:
                #     continue

                if site['cell'] in seen:
                    continue

                seen.add(site['cell'])

                output.append({
                    'radio': site['radio'],
                    'mcc': site['mcc'],
                    'net': site['net'],
                    'area': site['area'],
                    'cell': site['cell'],
                    'unit': site['unit'],
                    'lon': site['lon'],
                    'lat': site['lat'],
                    # 'range': site['range'],
                    # 'samples': site['samples'],
                    # 'changeable': site['changeable'],
                    # 'created': site['created'],
                    # 'updated': site['updated'],
                    # 'averageSignal': site['averageSignal']
                })
            # if len(output) > 0:
            #     break

    if len(output) == 0:
        return

    output = pd.DataFrame(output)
    output.to_csv(path_csv, index=False)

    return


def create_national_sites_shp(iso3):
    """
    Create a national sites csv layer for a selected country.

    """
    filename = '{}.csv'.format(iso3)
    folder = os.path.join(DATA_PROCESSED, iso3, 'sites')
    path_csv = os.path.join(folder, filename)

    filename = '{}.shp'.format(iso3)
    path_shp = os.path.join(folder, filename)

    if not os.path.exists(path_shp):

        # print('-Writing site shapefile data for {}'.format(iso3))

        country_data = pd.read_csv(path_csv)#[:10]

        output = []

        for idx, row in country_data.iterrows():
            output.append({
                'type': 'Feature',
                'geometry': {
                    'type': 'Point',
                    'coordinates': [row['lon'],row['lat']]
                },
                'properties': {
                    'radio': row['radio'],
                    'mcc': row['mcc'],
                    'net': row['net'],
                    'area': row['area'],
                    'cell': row['cell'],
                }
            })

        output = gpd.GeoDataFrame.from_features(output, crs='epsg:4326')

        output.to_file(path_shp)


def process_acled_layer(iso3):
    """
    Process the acled data by country.

    """
    filename = 'acled_{}.shp'.format(iso3)
    folder = os.path.join(DATA_PROCESSED, iso3, 'acled')
    path_out = os.path.join(folder, filename)

    # if os.path.exists(path_out):
    #     return

    if not os.path.exists(folder):
        os.makedirs(folder)

    filename = 'national_outline.shp'
    folder = os.path.join(DATA_PROCESSED, iso3)
    path_in = os.path.join(folder, filename)
    national_outline = gpd.read_file(path_in, crs='epsg:4326')

    filename = 'D1_Infrastructure_Database.shp'
    folder = os.path.join(DATA_RAW, "acled")
    path_in = os.path.join(folder, filename)
    data = gpd.read_file(path_in, crs='epsg:4326')

    data = gpd.overlay(data, national_outline, how='intersection')

    data = data[data['infra_cate'] == 'communications']

    data.to_file(path_out, crs='epsg:4326')

    return 


def intersect_acled(iso3):
    """
    Subset acled telecom incidents. 
    
    """
    output = []

    filename = 'acled_{}.shp'.format(iso3)
    folder = os.path.join(DATA_PROCESSED, iso3, 'acled')
    path_in = os.path.join(folder, filename)
    data = gpd.read_file(path_in, crs='epsg:4326')#[:10]
    data = data.to_crs(3857)
    data['geometry'] = data['geometry'].buffer(10000)

    filename = '{}.shp'.format(iso3)
    folder = os.path.join(DATA_PROCESSED, iso3, 'sites')
    path_in = os.path.join(folder, filename)
    cells = gpd.read_file(path_in, crs='epsg:4326') 
    cells = cells.to_crs(3857)

    output = []

    for idx, item in data.iterrows():
        for idx, cell in cells.iterrows():
            if item['geometry'].intersects(cell['geometry']):
                output.append({
                    'year': item['year'],
                    'month': item['month'],
                    'country_1': item['country_1'],
                    'admin1': item['admin1'],
                    'location': item['location'],
                    'infra_cate': item['infra_cate'],
                    'infra': item['infra'],
                    'sub_event_': item['sub_event_'],
                    'damage_inf': item['damage_inf'],
                    'severity_i': item['severity_i'],
                    'notes': item['notes'],
                    'fatalities': item['fatalities'],
                    'actor1': item['actor1'],
                    'assoc_acto': item['assoc_acto'],
                    'inter1': item['inter1'],
                    'actor2': item['actor2'],
                    'assoc_ac_1': item['assoc_ac_1'],
                    'inter2': item['inter2'],
                    'latitude': item['latitude'],
                    'longitude': item['longitude'],
                    'coordinate': item['coordinate'],
                    'source': item['source'],
                    'radio': cell['radio'],
                    'mcc': cell['mcc'],
                    'net': cell['net'],
                    'area': cell['area'],
                })                

    if len(output) == 0:
        return

    output = pd.DataFrame.from_records(output)

    filename = 'acled_{}.csv'.format(iso3)
    folder = os.path.join(BASE_PATH, '..', 'results', iso3)
    if not os.path.exists(folder):
        os.makedirs(folder)
    path_out = os.path.join(folder, filename)

    output.to_csv(path_out, index=False)

    return


def process_scdi_layer(iso3):
    """
    Process the scdi data by country.

    """
    filename = 'scdi_{}.shp'.format(iso3)
    folder = os.path.join(DATA_PROCESSED, iso3, 'scdi')
    path_out = os.path.join(folder, filename)

    if os.path.exists(path_out):
        return

    if not os.path.exists(folder):
        os.makedirs(folder)

    filename = 'national_outline.shp'
    folder = os.path.join(DATA_PROCESSED, iso3)
    path_in = os.path.join(folder, filename)
    national_outline = gpd.read_file(path_in, crs='epsg:4326')

    filename = 'SCDiDec2022.shp'
    folder = os.path.join(DATA_RAW, "OECD-SCDi-main", "data")
    path_in = os.path.join(folder, filename)
    data = gpd.read_file(path_in, crs='epsg:102024')
    data = data.to_crs(4326)

    data = gpd.overlay(data, national_outline, how='intersection')

    data.to_file(path_out, crs='epsg:4326')

    return 


if __name__ == "__main__":

    countries = get_countries()

    failures = []
    for idx, country in countries.iterrows():

        if not country['iso3'] in ['MLI', 'NER', 'BFA']:#, ]:
           continue

        print('Working on {}'.format(country['iso3']))

        run_preprocessing(country['iso3'])



