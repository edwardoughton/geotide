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
from shapely.geometry import Point, box, LineString, mapping
import rasterio
from rasterio.mask import mask
from rasterstats import zonal_stats
import numpy as np
from shapely.ops import cascaded_union
from geovoronoi import voronoi_regions_from_coords, points_to_coords

from tqdm import tqdm

from misc import get_countries, process_country_shapes, process_regions, get_regions 

CONFIG = configparser.ConfigParser()
CONFIG.read(os.path.join(os.path.dirname(__file__),'..', 'scripts', 'script_config.ini'))
BASE_PATH = CONFIG['file_locations']['base_path']

DATA_RAW = os.path.join(BASE_PATH, 'raw')
DATA_PROCESSED = os.path.join(BASE_PATH, 'processed')
RESULTS = os.path.join(BASE_PATH, '..', 'results')


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
    path = os.path.join(BASE_PATH, '..', '..', 'data_raw', filename)

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
    national_outline = national_outline[['geometry']]#,'iso3','iso2','gid_region']]

    filename = 'regions_{}_{}.shp'.format(2, iso3)
    folder = os.path.join(DATA_PROCESSED, iso3, 'regions')
    path_in = os.path.join(folder, filename)
    regions = gpd.read_file(path_in, crs='epsg:4326')#[:1]

    filename = 'D1_Infrastructure_Database.shp'
    folder = os.path.join(DATA_RAW, "acled")
    path_in = os.path.join(folder, filename)
    data = gpd.read_file(path_in, crs='epsg:4326')

    data = gpd.overlay(data, national_outline, how='intersection')

    data = data[data['infra_cate'] == 'communications']

    data = gpd.overlay(data, regions, how='intersection')

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
    data = gpd.read_file(path_in, crs='epsg:4326')#[:5]
    data = data.to_crs(3857)
    data['geometry'] = data['geometry'].buffer(10000)

    filename = '{}.shp'.format(iso3)
    folder = os.path.join(DATA_PROCESSED, iso3, 'sites')
    path_in = os.path.join(folder, filename)
    cells = gpd.read_file(path_in, crs='epsg:4326') 
    cells = cells.to_crs(3857)

    interim = []
    seen = set()

    for idx, item in data.iterrows():

        for idx, cell in cells.iterrows():
            if item['geometry'].intersects(cell['geometry']):

                point1 = item['geometry'].centroid
                point2 = cell['geometry'].centroid
                distance = LineString([point1, point2]).length

                interim.append({
                    'GID_0': item['GID_0'],
                    'NAME_0': item['NAME_0'],
                    'GID_2': item['GID_2'],
                    'NAME_1': item['NAME_1'],
                    'year': item['year'],
                    'month': item['month'],
                    'day': item['day'],
                    'country': item['country'],
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
                    'distance_to_event_m': distance,
                    'space_time_id': "{}_{}_{}_{}_{}".format(
                        item['year'],
                        item['month'], 
                        item['day'], 
                        item['longitude'], 
                        item['latitude'] 
                    ),
                })
                seen.add(item['coordinate'])
        if not item['coordinate'] in list(seen):
            interim.append({
                    'GID_0': item['GID_0'],
                    'NAME_0': item['NAME_0'],
                    'GID_2': item['GID_2'],
                    'NAME_1': item['NAME_1'],
                    'year': item['year'],
                    'month': item['month'],
                    'day': item['day'],
                    'country': item['country'],
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
                    'radio': 'NA',
                    'mcc': 'NA',
                    'net': 'NA',
                    'area': 'NA',
                    'distance_to_event_m': 'NA',
                    'space_time_id': "{}_{}_{}_{}_{}".format(
                        item['year'], 
                        item['month'], 
                        item['day'], 
                        item['longitude'], 
                        item['latitude'] 
                    ),
                })

    if len(interim) == 0:
        return

    output = remove_duplicates(interim)

    output = pd.DataFrame.from_records(output)

    filename = 'acled_{}.csv'.format(iso3)
    folder = os.path.join(BASE_PATH, '..', 'results', iso3)
    if not os.path.exists(folder):
        os.makedirs(folder)
    path_out = os.path.join(folder, filename)

    output.to_csv(path_out, index=False)

    return


def remove_duplicates(interim):
    """
    This function removes duplicated cells.
    
    Only the 3 nearest cells are selected for each radio, using the distance estimate.  
    
    """
    output = []

    coordinates = set()
    for item in interim:
        coordinates.add(item['coordinate'])

    radios = ['GSM', 'UMTS', 'LTE'] 

    for coordinate in list(coordinates):
        
        cells = []
        
        #first get cells matching coordinate
        for item in interim:
            if coordinate == item['coordinate']:
                cells.append(item)

        #then subset by radio type
        for radio in radios:

            cells_by_radio = []
            
            for item in cells:
                if radio == item['radio']:
                    cells_by_radio.append(item)

            if not len(cells_by_radio) > 0:
                continue

            cells_by_radio = pd.DataFrame(cells_by_radio)
            cells_by_radio = cells_by_radio.sort_values(by=['distance_to_event_m'], ascending=True)
            #now get 3 closest cells
            cells_by_radio = cells_by_radio[:3]
            cells_by_radio = cells_by_radio.to_dict('records')
            output = output + cells_by_radio

        for item in cells:
            if item['distance_to_event_m'] == 'NA':
                output.append(item)
                output.append(item)
                output.append(item)

    return output


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


def count_cells(country):
    """
    Add up cells by technology.
    
    """
    regional_level = country['gid_region']
    gid_level = 'gid_{}'.format(regional_level) #regional_level
    folder = os.path.join(DATA_PROCESSED, country['iso3'], 'sites', gid_level)

    regions_df = get_regions(country, regional_level)

    if len(regions_df) == 0:
        print("len(regions_df) == 0")
        return

    output = []

    for idx, region in regions_df.iterrows():

        region = region["GID_{}".format(regional_level)]
        
        filename = "{}.csv".format(region)
        path_in = os.path.join(folder, filename)
        
        if not os.path.exists(path_in):
            print("path_in did not exist: {}".format(path_in))
            continue

        data = pd.read_csv(path_in)

        cells_2g = 0
        cells_3g = 0
        cells_4g = 0
        cells_5g = 0

        for idx, row in data.iterrows():
            
            if row['radio'] == 'GSM':
                cells_2g += 1
            elif row['radio'] == 'UMTS':
                cells_3g += 1
            elif row['radio'] == 'LTE':
                cells_4g += 1
            elif row['radio'] == 'NR':
                cells_5g += 1

        output.append({
            "iso3": country['iso3'],
            "gid_level": gid_level,
            "gid_id": region,
            "cells_2g": cells_2g,
            "cells_3g": cells_3g,
            "cells_4g": cells_4g,
            "cells_5g": cells_5g,
        })
        
    if not len(output) > 0:
        return

    output = pd.DataFrame(output)

    filename = "cell_count.csv"
    folder_out = os.path.join(DATA_PROCESSED, country['iso3'], 'sites')
    path_out = os.path.join(folder_out, filename)
    output.to_csv(path_out, index=False)

    return


def collect_cells(countries):
    """
    Collect all cells.
    """
    output = []

    for country in countries:

        filename = "cell_count.csv"
        folder_in = os.path.join(DATA_PROCESSED, country['iso3'], 'sites')
        path_in = os.path.join(folder_in, filename)

        if not os.path.exists(path_in):
            continue
        print(path_in)
        print("Collecting cells for {}".format(country['iso3']))

        data = pd.read_csv(path_in)
        data = data.to_dict('records')
        output = output + data

    if not len(output) > 0:
        return

    output = pd.DataFrame(output)

    filename = "cell_count_regional.csv"
    folder_out = os.path.join(DATA_PROCESSED, 'results', 'sites')
    if not os.path.exists(folder_out):
        os.makedirs(folder_out)
    path_out = os.path.join(folder_out, filename)
    output.to_csv(path_out, index=False)


def process_settlement_layer(country):
    """
    Clip the settlement layer to the chosen country boundary
    and place in desired country folder.

    Parameters
    ----------
    country : dict
        Contains all desired country information.

    """
    iso3 = country['iso3']
    regional_level = country['gid_region']

    filename = 'ppp_2020_1km_Aggregated.tif'
    folder = os.path.join(BASE_PATH,'..', '..', 'data_raw', 'settlement_layer')
    path_settlements = os.path.join(folder, filename)

    settlements = rasterio.open(path_settlements, 'r+')
    settlements.nodata = 255
    settlements.crs.from_epsg(4326)

    iso3 = country['iso3']
    path_country = os.path.join(DATA_PROCESSED, iso3,
        'national_outline.shp')

    if os.path.exists(path_country):
        country = gpd.read_file(path_country)
    else:
        print('Must generate national_outline.shp first' )

    path_country = os.path.join(DATA_PROCESSED, iso3)
    shape_path = os.path.join(path_country, 'settlements.tif')

    if os.path.exists(shape_path):
        return #print('Completed settlement layer processing')

    # print('----')
    # print('Working on {} level {}'.format(iso3, regional_level))

    # bbox = country.envelope
    geo = gpd.GeoDataFrame()

    geo = gpd.GeoDataFrame({'geometry': country['geometry']})

    coords = [json.loads(geo.to_json())['features'][0]['geometry']]

    out_img, out_transform = mask(settlements, coords, crop=True)

    out_meta = settlements.meta.copy()

    out_meta.update({"driver": "GTiff",
                    "height": out_img.shape[1],
                    "width": out_img.shape[2],
                    "transform": out_transform,
                    "crs": 'epsg:4326'})

    with rasterio.open(shape_path, "w", **out_meta) as dest:
            dest.write(out_img)

    return #print('Completed processing of settlement layer')


def generate_agglomerations(country):
    """
    Generate a lookup table of agglomerations.

    """
    iso3 = country['iso3']
    regional_level = country['gid_region']
    GID_level = 'GID_{}'.format(regional_level)

    folder = os.path.join(DATA_PROCESSED, iso3, 'agglomerations')
    if not os.path.exists(folder):
        os.makedirs(folder)
    path_output = os.path.join(folder, 'agglomerations.shp')

    if os.path.exists(path_output):
        return print('Agglomeration processing has already completed')

    print('Working on {} agglomeration lookup table'.format(iso3))

    filename = 'regions_{}_{}.shp'.format(regional_level, iso3)
    folder = os.path.join(DATA_PROCESSED, iso3, 'regions')
    path = os.path.join(folder, filename)
    regions = gpd.read_file(path, crs="epsg:4326")

    path_settlements = os.path.join(DATA_PROCESSED, iso3, 'settlements.tif')
    settlements = rasterio.open(path_settlements, 'r+')
    settlements.nodata = 255
    settlements.crs = {"init": "epsg:4326"}

    folder_tifs = os.path.join(DATA_PROCESSED, iso3, 'agglomerations', 'tifs')
    if not os.path.exists(folder_tifs):
        os.makedirs(folder_tifs)

    for idx, region in regions.iterrows():

        # bbox = region['geometry'].envelope
        geo = gpd.GeoDataFrame()
        geo = gpd.GeoDataFrame({'geometry': region['geometry']}, index=[idx])
        coords = [json.loads(geo.to_json())['features'][0]['geometry']]

        #chop on coords
        out_img, out_transform = mask(settlements, coords, crop=True)

        # Copy the metadata
        out_meta = settlements.meta.copy()

        out_meta.update({"driver": "GTiff",
                        "height": out_img.shape[1],
                        "width": out_img.shape[2],
                        "transform": out_transform,
                        "crs": 'epsg:4326'})

        path_output = os.path.join(folder_tifs, region[GID_level] + '.tif')

        with rasterio.open(path_output, "w", **out_meta) as dest:
                dest.write(out_img)

    print('Completed settlement.tif regional segmentation')

    nodes, missing_nodes = find_nodes(country, regions)

    missing_nodes = get_missing_nodes(country, regions, missing_nodes, 10, 10)

    nodes = nodes + missing_nodes

    nodes = gpd.GeoDataFrame.from_features(nodes, crs='epsg:4326')

    bool_list = nodes.intersects(regions['geometry'].unary_union)
    nodes = pd.concat([nodes, bool_list], axis=1)
    nodes = nodes[nodes[0] == True].drop(columns=0)

    agglomerations = []

    print('Identifying agglomerations')
    for idx1, region in regions.iterrows():
        seen = set()
        for idx2, node in nodes.iterrows():
            if node['geometry'].intersects(region['geometry']):
                agglomerations.append({
                    'type': 'Feature',
                    'geometry': mapping(node['geometry']),
                    'properties': {
                        'id': idx1,
                        'GID_0': region['GID_0'],
                        GID_level: region[GID_level],
                        'population': node['sum'],
                    }
                })
                seen.add(region[GID_level])
        if len(seen) == 0:
            agglomerations.append({
                    'type': 'Feature',
                    'geometry': mapping(region['geometry'].centroid),
                    'properties': {
                        'id': 'regional_node',
                        'GID_0': region['GID_0'],
                        GID_level: region[GID_level],
                        'population': 1,
                    }
                })

    agglomerations = gpd.GeoDataFrame.from_features(
            [
                {
                    'geometry': item['geometry'],
                    'properties': {
                        'id': item['properties']['id'],
                        'GID_0':item['properties']['GID_0'],
                        GID_level: item['properties'][GID_level],
                        'population': item['properties']['population'],
                    }
                }
                for item in agglomerations
            ],
            crs='epsg:4326'
        )

    folder = os.path.join(DATA_PROCESSED, iso3, 'agglomerations')
    path_output = os.path.join(folder, 'agglomerations' + '.shp')

    agglomerations.to_file(path_output)

    agglomerations['lon'] = agglomerations['geometry'].x
    agglomerations['lat'] = agglomerations['geometry'].y
    agglomerations = agglomerations[['lon', 'lat', GID_level, 'population']]
    agglomerations.to_csv(os.path.join(folder, 'agglomerations.csv'), index=False)

    return print('Agglomerations layer complete')


def find_nodes(country, regions):
    """
    Find key nodes.

    """
    iso3 = country['iso3']
    regional_level = country['gid_region']
    GID_level = 'GID_{}'.format(regional_level)

    threshold = country['pop_density_km2']
    settlement_size = country['settlement_size']

    folder_tifs = os.path.join(DATA_PROCESSED, iso3, 'agglomerations', 'tifs')

    interim = []
    missing_nodes = set()

    print('Working on gathering data from regional rasters')
    for idx, region in regions.iterrows():

        path = os.path.join(folder_tifs, region[GID_level] + '.tif')

        with rasterio.open(path) as src:
            data = src.read()
            data[data < threshold] = 0
            data[data >= threshold] = 1
            polygons = rasterio.features.shapes(data, transform=src.transform)
            shapes_df = gpd.GeoDataFrame.from_features(
                [
                    {'geometry': poly, 'properties':{'value':value}}
                    for poly, value in polygons
                    if value > 0
                ],
                crs='epsg:4326'
            )

        geojson_region = [
            {
                'geometry': region['geometry'],
                'properties': {
                    GID_level: region[GID_level]
                }
            }
        ]

        gpd_region = gpd.GeoDataFrame.from_features(
                [
                    {'geometry': poly['geometry'],
                    'properties':{
                        GID_level: poly['properties'][GID_level]
                        }}
                    for poly in geojson_region
                ], crs='epsg:4326'
            )

        if len(shapes_df) == 0:
            continue

        nodes = gpd.overlay(shapes_df, gpd_region, how='intersection')

        stats = zonal_stats(shapes_df['geometry'], path, stats=['count', 'sum'])

        stats_df = pd.DataFrame(stats)

        nodes = pd.concat([shapes_df, stats_df], axis=1).drop(columns='value')

        nodes_subset = nodes[nodes['sum'] >= settlement_size]

        if len(nodes_subset) == 0:
            missing_nodes.add(region[GID_level])

        for idx, item in nodes_subset.iterrows():
            interim.append({
                    'geometry': item['geometry'].centroid,
                    'properties': {
                        GID_level: region[GID_level],
                        'count': item['count'],
                        'sum': item['sum']
                    }
            })

    return interim, missing_nodes


def get_missing_nodes(country, regions, 
    missing_nodes, threshold, settlement_size):
    """
    Find any missing nodes.

    """
    iso3 = country['iso3']
    regional_level = country['gid_region']
    GID_level = 'GID_{}'.format(regional_level)

    folder_tifs = os.path.join(DATA_PROCESSED, iso3, 'agglomerations', 'tifs')

    interim = []

    for idx, region in regions.iterrows():

        if not region[GID_level] in list(missing_nodes):
            continue

        path = os.path.join(folder_tifs, region[GID_level] + '.tif')

        with rasterio.open(path) as src:
            data = src.read()
            data[data < threshold] = 0
            data[data >= threshold] = 1
            polygons = rasterio.features.shapes(data, transform=src.transform)
            shapes_df = gpd.GeoDataFrame.from_features(
                [
                    {'geometry': poly, 'properties':{'value':value}}
                    for poly, value in polygons
                    if value > 0
                ],
                crs='epsg:4326'
            )

        geojson_region = [
            {
                'geometry': region['geometry'],
                'properties': {
                    GID_level: region[GID_level]
                }
            }
        ]

        gpd_region = gpd.GeoDataFrame.from_features(
                [
                    {'geometry': poly['geometry'],
                    'properties':{
                        GID_level: poly['properties'][GID_level]
                        }}
                    for poly in geojson_region
                ], crs='epsg:4326'
            )

        nodes = gpd.overlay(shapes_df, gpd_region, how='intersection')

        stats = zonal_stats(shapes_df['geometry'], path, stats=['count', 'sum'])

        stats_df = pd.DataFrame(stats)

        nodes = pd.concat([shapes_df, stats_df], axis=1).drop(columns='value')

        max_sum = nodes['sum'].max()

        nodes = nodes[nodes['sum'] > max_sum - 1]

        for idx, item in nodes.iterrows():
            interim.append({
                    'geometry': item['geometry'].centroid,
                    'properties': {
                        GID_level: region[GID_level],
                        'count': item['count'],
                        'sum': item['sum']
                    }
            })

    return interim


def generate_distance_lut(country):
    """
    Generate distance lut from each site to the nearest agglomeration. 
    
    """
    filename = '{}.shp'.format(country['iso3'])
    folder = os.path.join(DATA_PROCESSED, country['iso3'], 'sites')
    path = os.path.join(folder, filename)
    sites = gpd.read_file(path, crs='epsg:4326')
    sites = sites.to_crs(3857)#[:1]

    filename = 'agglomerations.shp'
    folder = os.path.join(DATA_PROCESSED, country['iso3'], 'agglomerations')
    path = os.path.join(folder, filename)
    agglomerations = gpd.read_file(path, crs='epsg:4326')
    agglomerations = agglomerations.to_crs(3857)

    output = []

    for idx, site in sites.iterrows():
        
        geom_site = site['geometry']
        
        distance_lut = []

        for idx, agglomeration in agglomerations.iterrows():

            geom_agglomeration = agglomeration['geometry']

            line = LineString([
                geom_site,
                geom_agglomeration
            ])

            distance_lut.append({
                'gid_id': agglomeration['GID_1'],
                'radio': site['radio'],
                'mcc': site['mcc'],
                'net': site['net'],
                'area': site['area'],
                'cell': site['cell'],
                'coords_site': '{}_{}'.format(geom_site.x, geom_site.y),
                'coords_agglomeration': '{}_{}'.format(geom_agglomeration.x, geom_agglomeration.y),
                'distance_m': line.length,
            })

        distance_lut = pd.DataFrame(distance_lut)
        distance_lut = distance_lut.sort_values(by=['distance_m'], ascending=True)
        distance_lut = distance_lut[:1]
        distance_lut = distance_lut.to_dict('records')

        output = output + distance_lut

    output = pd.DataFrame(output)

    filename = 'site_to_agglomeration_lut.csv'
    folder = os.path.join(DATA_PROCESSED, country['iso3'], 'agglomerations')
    path = os.path.join(folder, filename)
    output.to_csv(path, index=False)

    return


def concat_acled_results():
    """
    Generate single .csv file. 

    """
    output = []

    for subfolder in os.listdir(RESULTS):

        if subfolder.endswith(".csv"):
            continue
        if subfolder.endswith("xlsx"):
            continue
        path = os.path.join(RESULTS, subfolder, "acled_{}.csv".format(subfolder))
        data = pd.read_csv(path)
        data = data.to_dict('records')
        output = output + data

    output = pd.DataFrame(output)
    output.to_csv(os.path.join(RESULTS, "all_results.csv"), index=False)

    return


def generate_site_coverage_areas(iso3):
    """
    Generate initial site coverage areas events layer.  

    """
    filename = "acled_{}.shp".format(iso3)
    folder = os.path.join(DATA_PROCESSED, iso3, 'acled')
    path_in = os.path.join(folder, filename)
    data = gpd.read_file(path_in, crs='epsg:4326')
    data = data.to_crs(3857)
    data['geometry'] = data['geometry'].buffer(10000)
    data = data.to_crs(4326)

    filename = "initial_coverage_polygons.shp"
    folder = os.path.join(DATA_PROCESSED, iso3, 'coverage_polygons', 'interim')
    if not os.path.exists(folder):
        os.makedirs(folder)
    path_out = os.path.join(folder, filename)
    data.to_file(path_out, crs='epsg:4326') 

    return


def clean_coverage_areas(iso3):
    """
    Clean site coverage areas events layer.  

    """
    
    filename = "coverage_polygon_boundaries.shp"
    folder = os.path.join(DATA_PROCESSED, iso3, 'coverage_polygons', 'interim')
    path_out = os.path.join(folder, filename)

    # if not os.path.exists(path):
    filename = "initial_coverage_polygons.shp"
    folder = os.path.join(DATA_PROCESSED, iso3, 'coverage_polygons', 'interim')
    path = os.path.join(folder, filename)
    boundaries = gpd.read_file(path, crs='epsg:4326') 
    boundaries = boundaries[['geometry', 'fid']]
    union_boundaries = boundaries['geometry'].unary_union

    output = []

    for union_boundary in union_boundaries:
        for idx, boundary in boundaries.iterrows():
            if union_boundary.intersects(boundary['geometry']):
                output.append({
                    'geometry': union_boundary,
                    'properties': {
                        'fid': boundary['fid']
                    }
                })

    output = gpd.GeoDataFrame.from_features(output, crs='epsg:4326')
    output.to_file(path_out, crs='epsg:4326')
    boundaries = output.to_dict('records')

    # # else:
    # #     boundaries = gpd.read_file(path, crs='epsg:4326')
    # #     boundaries = boundaries.to_dict('records')

    filename = "acled_{}.shp".format(iso3)
    folder = os.path.join(DATA_PROCESSED, iso3, 'acled')
    path_in = os.path.join(folder, filename)
    points = gpd.read_file(path_in, crs='epsg:4326')
    points = points.drop_duplicates(subset='geometry')
    points = points.to_dict('records')

    interim = []

    for boundary in boundaries:

        # if not boundary['FID'] == 38:
        #     continue

        ids = []
        geoms = []

        for point in points:
            if point['geometry'].intersects(boundary['geometry']):
                ids.append(point['id'])
                geoms.append(point['geometry'])

        geoms = gpd.GeoDataFrame([geom for geom in geoms]).set_geometry(0)

        coords = points_to_coords(geoms.geometry)

        only_2_coords = []

        if len(coords) > 2:
            polys, pts = voronoi_regions_from_coords(coords, boundary['geometry'])
            for key, value in polys.items():
                interim.append(value)
        # elif len(coords) == 2:
        #     only_2_coords.append(boundary['FID'])
        elif len(coords) == 1: 
            interim.append(boundary['geometry'])

    output = []
    seen = []

    for boundary in interim:
        for point in points:
            if point['geometry'].intersects(boundary):
                output.append({
                    'geometry': boundary,
                    'properties': {
                        'GID_0': point['GID_0'],
                        'NAME_0': point['NAME_0'],
                        'GID_2': 'GID_2',
                        'NAME_1': 'NAME_1',
                        'year': point['year'],
                        'month': point['month'],
                        'day': point['day'],
                        'country': point['country'],
                        'admin1': point['admin1'],
                        'location': point['location'],
                        'infra_cate': point['infra_cate'],
                        'infra': point['infra'],
                        'sub_event_': point['sub_event_'],
                        'damage_inf': point['damage_inf'],
                        'severity_i': point['severity_i'],
                        'notes': point['notes'],
                        'fatalities': point['fatalities'],
                        'actor1': point['actor1'],
                        'assoc_acto': point['assoc_acto'],
                        'inter1': point['inter1'],
                        'actor2': point['actor2'],
                        'assoc_ac_1': point['assoc_ac_1'],
                        'inter2': point['inter2'],
                        'latitude': point['latitude'],
                        'longitude': point['longitude'],
                        'coordinate': point['coordinate'],
                        'source': point['source'],
                    }
                })
                seen.append(point['id'])

    # # add back in any with only 2 points, as the Voronoi tesselation 
    # # will not work with <3 points (requirement of delaunay triangulation)
    # filename = "initial_coverage_polygons.shp"
    # folder = os.path.join(DATA_PROCESSED, iso3, 'coverage_polygons', 'interim')
    # path = os.path.join(folder, filename)
    # boundaries = gpd.read_file(path, crs='epsg:4326')
    # boundaries = boundaries.to_dict('records')

    # for item in boundaries:
    #     if not item['id'] in seen:
    #         output.append({
    #             'geometry': item['geometry'],
    #             'properties': {
    #                 'GID_0': point['GID_0'],
    #                 'NAME_0': point['NAME_0'],
    #                 'GID_2': 'GID_2',
    #                 'NAME_1': 'NAME_1',
    #                 'year': point['year'],
    #                 'month': point['month'],
    #                 'day': point['day'],
    #                 'country': point['country'],
    #                 'admin1': point['admin1'],
    #                 'location': point['location'],
    #                 'infra_cate': point['infra_cate'],
    #                 'infra': point['infra'],
    #                 'sub_event_': point['sub_event_'],
    #                 'damage_inf': point['damage_inf'],
    #                 'severity_i': point['severity_i'],
    #                 'notes': point['notes'],
    #                 'fatalities': point['fatalities'],
    #                 'actor1': point['actor1'],
    #                 'assoc_acto': point['assoc_acto'],
    #                 'inter1': point['inter1'],
    #                 'actor2': point['actor2'],
    #                 'assoc_ac_1': point['assoc_ac_1'],
    #                 'inter2': point['inter2'],
    #                 'latitude': point['latitude'],
    #                 'longitude': point['longitude'],
    #                 'coordinate': point['coordinate'],
    #                 'source': point['source'],
    #             }
    #         })

    filename = "coverage_polygons.shp"
    folder = os.path.join(DATA_PROCESSED, iso3, 'coverage_polygons')
    path = os.path.join(folder, filename)

    output = gpd.GeoDataFrame.from_features(output, crs='epsg:4326')
    output.to_file(path, crs='epsg:4326')

    return


def process_rwi_geometry(iso3):
    """
    Process rwi .csv layer to points.

    """
    filename = 'regions_{}_{}.shp'.format(2, iso3)
    path_in = os.path.join(BASE_PATH, 'processed', iso3, 'regions', filename)
    regions = gpd.read_file(path_in, crs="epsg:4326")

    filename = '{}_relative_wealth_index.csv'.format(iso3)
    folder = os.path.join(BASE_PATH, '..', '..', 'data_raw', 'relative_wealth_index')
    path_in = os.path.join(folder, filename)
    data = gpd.read_file(path_in, encoding='latin-1')

    points = gpd.GeoDataFrame(
        data, geometry=gpd.points_from_xy(data.longitude, data.latitude), crs="EPSG:4326"
    )  

    for idx, region in regions.iterrows():

        filename_out = '{}.shp'.format(region['GID_2']) #each regional file is named using the gid id
        folder_out = os.path.join(BASE_PATH, 'processed', iso3 , 'rwi', 'points')
        if not os.path.exists(folder_out):
            os.makedirs(folder_out)
        path_out = os.path.join(folder_out, filename_out)

        if os.path.exists(path_out):
            continue

        output = []

        for idx, point in points.iterrows():
            if point['geometry'].intersects(region['geometry']):
                output.append({
                    'type': 'Point',
                    'geometry': point['geometry'],
                    'properties': {
                        'rwi': point['rwi'],
                        'error': point['error'],
                    }
                })

        if len(output) == 0:
            continue

        output = gpd.GeoDataFrame.from_features(output, crs='epsg:4326')
        output.to_file(path_out,crs="epsg:4326")
    
    return


def create_rwi_voronoi(iso3):
    """
    Creates rwi voronoi polygons by region.

    """
    filename = 'regions_{}_{}.shp'.format(2, iso3)
    path_in = os.path.join(BASE_PATH, 'processed', iso3, 'regions', filename)
    regions = gpd.read_file(path_in, crs="epsg:4326")

    for idx, region in regions.iterrows():

        print('Working on vononoi polygons for: {}'.format(region['GID_2']))

        filename = '{}.shp'.format(region['GID_2'])
        path_in = os.path.join(BASE_PATH, 'processed', iso3, 'rwi', 'points', filename)
        if not os.path.exists(path_in):
            continue
        rwi_df = gpd.read_file(path_in, crs="epsg:4326")#[:500]

        boundary_shape = cascaded_union(region.geometry)
        coords = points_to_coords(rwi_df.geometry)

        region_polys, region_pts = voronoi_regions_from_coords(coords, boundary_shape)

        interim = []
        for key, value in region_polys.items():
            for idx, rwi in rwi_df.iterrows():
                if value.intersects(rwi['geometry']):
                    interim.append({ 
                        'type': 'Polygon',
                        'geometry': value,
                        'properties': {
                            "rwi": rwi['rwi'],
                            'error': rwi['error']
                        }
                    })
        output = gpd.GeoDataFrame.from_features(interim, crs="EPSG:4326")

        if len(output) == 0:
            continue

        filename = '{}.shp'.format(region['GID_2'])
        folder_out = os.path.join(BASE_PATH, 'processed', iso3, 'rwi', 'polygons')
        if not os.path.exists(folder_out):
            os.makedirs(folder_out)
        path_out = os.path.join(folder_out, filename)
        output.to_file(path_out, crs="EPSG:4326")

    return


def convert_settlement_to_polygon(iso3):
    """
    Convert to .shp.
    
    """
    filename = 'regions_{}_{}.shp'.format(2, iso3)
    path_in = os.path.join(BASE_PATH, 'processed', iso3, 'regions', filename)
    regions = gpd.read_file(path_in, crs="epsg:4326")

    for idx, region in regions.iterrows():

        filename = "{}.tif".format(region['GID_2'])
        folder = os.path.join(DATA_PROCESSED, iso3, 'settlements')
        path_in = os.path.join(folder, filename)

        if not os.path.exists(path_in):
            continue

        filename = "{}.shp".format(region['GID_2'])
        folder = os.path.join(DATA_PROCESSED, iso3, 'settlements', 'shapes')
        if not os.path.exists(folder):
            os.makedirs(folder)
        path_out = os.path.join(folder, filename)

        if os.path.exists(path_out):
            continue

        output = []

        mask = None
        with rasterio.Env():
            with rasterio.open(path_in) as src:
                image = src.read(1)#[:1] # first band
                results = (
                    {'properties': {'raster_val': v}, 'geometry': s}
                    for i, (s, v) in enumerate(
                        shapes(image, mask=mask, transform=src.transform))
                    )
                for item in results:
                    output.append({
                        'geometry': item['geometry'],
                        'properties': item['properties'],
                    })

        if len(output) == 0:
            continue

        output = gpd.GeoDataFrame.from_features(output)
        output.to_file(path_out, crs='epsg:4326')

    return


def intersect_coverage_areas_and_settlements(iso3):
    """
    Intersect the coverage area and population polygons. 
    
    """
    filename = 'regions_{}_{}.shp'.format(2, iso3)
    path_in = os.path.join(BASE_PATH, 'processed', iso3, 'regions', filename)
    regions = gpd.read_file(path_in, crs="epsg:4326")#[:1]

    for idx, region in regions.iterrows():

        if not region['GID_2'] == 'BFA.7.1_1':
            continue

        filename = '{}.shp'.format(region['GID_2'])
        folder_out = os.path.join(DATA_PROCESSED, iso3, 'estimated_income')
        if not os.path.exists(folder_out):
            os.mkdir(folder_out)
        path_out = os.path.join(folder_out, filename)
        
        # if os.path.exists(path_out):
        #     continue

        filename = "{}.shp".format(region['GID_2'])
        folder_in = os.path.join(DATA_PROCESSED, iso3, 'settlements', 'shapes')
        path_in = os.path.join(folder_in, filename)
        if not os.path.exists(path_in):
            continue
        data_pop = gpd.read_file(path_in, crs='epsg:4326')#[:100]

        filename = '{}.shp'.format(region['GID_2'])
        folder_in = os.path.join(DATA_PROCESSED, iso3, 'rwi', 'polygons')
        path_in = os.path.join(folder_in, filename)
        if not os.path.exists(path_in):
            continue
        data_rwi = gpd.read_file(path_in, crs='epsg:4326')#[:1]

        output = []

        for idx1, pop_tile in data_pop.iterrows():
            for idx2, rwi_poly in data_rwi.iterrows():
                geom = pop_tile['geometry'].centroid
                if geom.intersects(rwi_poly['geometry']):
                    output.append({
                        'geometry': pop_tile['geometry'], #.representative_point(),
                        'properties': {
                            'pop': pop_tile['raster_val'],
                            'rwi': rwi_poly['rwi'],
                            'error': rwi_poly['error'],
                        } 
                    })

        if len(output) == 0:
            continue

        output = gpd.GeoDataFrame.from_features(output)
        output.to_file(path_out, crs='epsg:4326')

    return


def estimate_customers_per_site(iso3):
    """
    Estimate the population (and their wealth) served by each cell.
    
    """
    path_settlements = os.path.join(DATA_PROCESSED, iso3,
        'settlements.tif')

    filename = 'coverage_polygons.shp'
    folder = os.path.join(DATA_PROCESSED, iso3, 'coverage_polygons')
    path = os.path.join(folder, filename)
    coverage_areas = gpd.read_file(path, crs='epsg:4326')#[:10]
    coverage_areas = coverage_areas.to_dict('records')

    output = []
    
    for coverage_area in coverage_areas:

        # if not coverage_area['coordinate'] == '12.3703,-1.5247':
        #     continue

        gid_ids = find_gid_ids(country, coverage_area)

        for gid_id in gid_ids:

            filename = '{}.shp'.format(gid_id)
            folder = os.path.join(DATA_PROCESSED, country['iso3'], 'estimated_income')
            path = os.path.join(folder, filename)
            data = gpd.read_file(path, crs='epsg:4326')

            coverage_area_df = gpd.GeoDataFrame(
                geometry=gpd.GeoSeries(coverage_area['geometry']), crs='epsg:4326')

            coverage_area_df = gpd.overlay(coverage_area_df, data, how='intersection')
            coverage_area_df = coverage_area_df.to_crs(3857)
            coverage_area_df = coverage_area_df.to_dict('records')

            for item in coverage_area_df:

                area_km2 = item['geometry'].area / 1e6

                if area_km2 > 1.5:
                    continue
                
                pop = item['pop'] * area_km2
                
                if pop < 0:
                    pop = 0
                
                try: 
                    geom = item['geometry'].representative_point()
                    lat_tile = geom.y
                    lon_tile = geom.x
                except:
                    lat_tile = 'NA'
                    lon_tile = 'NA'

                output.append({
                    'geometry': item['geometry'],
                    'properties': {
                        'GID_0': coverage_area['GID_0'],
                        'NAME_0': coverage_area['NAME_0'],
                        'GID_2': coverage_area['GID_2'],
                        'NAME_1': coverage_area['NAME_1'],
                        'year': coverage_area['year'],
                        'month': coverage_area['month'],
                        'day': coverage_area['day'],
                        'country': coverage_area['country'],
                        'admin1': coverage_area['admin1'],
                        'location': coverage_area['location'],
                        'infra_cate': coverage_area['infra_cate'],
                        'infra': coverage_area['infra'],
                        'sub_event_': coverage_area['sub_event_'],
                        'damage_inf': coverage_area['damage_inf'],
                        'severity_i': coverage_area['severity_i'],
                        'notes': coverage_area['notes'],
                        'fatalities': coverage_area['fatalities'],
                        'actor1': coverage_area['actor1'],
                        'assoc_acto': coverage_area['assoc_acto'],
                        'inter1': coverage_area['inter1'],
                        'actor2': coverage_area['actor2'],
                        'assoc_ac_1': coverage_area['assoc_ac_1'],
                        'inter2': coverage_area['inter2'],
                        'lat_site': coverage_area['latitude'],
                        'lon_site': coverage_area['longitude'],
                        'coords_site': coverage_area['coordinate'],
                        'lat_tile': lat_tile,
                        'lon_tile': lon_tile,
                        'source': coverage_area['source'],
                        'pop': pop,
                        'rwi': item['rwi'],
                        'error': item['error'],
                    }
                })

    if len(output) == 0:
        return

    output = gpd.GeoDataFrame.from_features(output, crs='epsg:3857')#[:200]

    output = output.to_crs(4326)
    filename = 'customers_per_site.shp'
    folder_out = os.path.join(DATA_PROCESSED, iso3, 'customers_per_site')
    if not os.path.exists(folder_out):
        os.makedirs(folder_out)
    path_out = os.path.join(folder_out, filename)
    output.to_file(path_out, crs='epsg:4326')

    return


def find_gid_ids(country, coverage_area):
    """
    Find gid_ids intersecting with each site.
    
    """
    coverage_bounds = coverage_area['geometry'].bounds
    c_minx, c_miny, c_maxx, c_maxy = coverage_bounds
    c_geom = box(*coverage_bounds, ccw=True)

    filename = 'regions_{}_{}.shp'.format(2, country['iso3'])
    folder = os.path.join(DATA_PROCESSED, country['iso3'], 'regions')
    path = os.path.join(folder, filename)
    regions = gpd.read_file(path, crs='epsg:4326')
    regions = regions.to_dict('records')

    output = []

    for region in regions:

        region_bounds = region['geometry'].bounds
        r_geom = box(*region_bounds, ccw=True)
        
        if  coverage_area['geometry'].intersects(region['geometry']):
            output.append(region['GID_2'])

    return output


def export_to_csv(iso3):
    """
    Export shapefiles to .csv
    
    """
    filename = 'customers_per_site.shp'
    folder = os.path.join(DATA_PROCESSED, iso3, 'customers_per_site')
    path = os.path.join(folder, filename)
    data = gpd.read_file(path, crs='epsg:4326')#[:10]
    del data['geometry']

    data['space_time_id'] = (
        data['year'].astype(str) + '_' +
        data['month'].astype(str) + '_' +
        data['day'].astype(str) + '_' +
        data['lon_site'].astype(str) + '_' +
        data['lat_site'].astype(str) 
    )

    # drop any broken tiles with pop < 0 
    data = data[data['pop'] >= 0]

    # get the pop for each site using the space_time_id
    data = data[['space_time_id', 'pop', 'rwi', 'error','lon_tile','lat_tile']]
    data = data.drop_duplicates()

    filename = 'customers_per_site.csv'
    path_out = os.path.join(RESULTS, iso3, filename)
    data.to_csv(path_out, index=False)

    return


if __name__ == "__main__":

    countries = get_countries()
    countries = countries.to_dict('records')

    for country in tqdm(countries):

        if not country['iso3'] in ['BFA']:#, 'MLI', 'NER']:#]:#, 
           continue

        print('---Working on {}'.format(country['iso3']))

        filename = "countries.csv"
        path = os.path.join(DATA_RAW, filename)

        countries = pd.read_csv(path, encoding='latin-1')
        country = countries[countries.iso3 == country['iso3']]
        country = country.to_records('dicts')[0]
        regional_level = int(country['gid_region'])

        print('Working on create_national_sites_csv')
        create_national_sites_csv(country)

        print('Working on process_country_shapes')
        process_country_shapes(country['iso3'])

        print('Working on process_regions')
        process_regions(country['iso3'], 2)

        print('Working on create_national_sites_shp')
        create_national_sites_shp(country['iso3'])

        print('Working on process_acled_layer')
        process_acled_layer(country['iso3'])

        print('Working on subset_acled_telecom')
        intersect_acled(country['iso3'])

        print('Working on process_scdi_layer')
        process_scdi_layer(country['iso3'])

        print('Processing count_cells')
        count_cells(country)

        print('Processing process_settlement_layer')
        process_settlement_layer(country)

        print('Working on generate_agglomerations')
        generate_agglomerations(country)

        print('Working on generate_distance_lut')
        generate_distance_lut(country)
        
        print('---Generate initial site coverage areas')
        generate_site_coverage_areas(country['iso3'])

        print('---Clean site coverage areas')
        clean_coverage_areas(country['iso3'])

        print('---Processing rwi geometry')
        process_rwi_geometry(country['iso3'])

        print('---Processing rwi voronoi polygons')
        create_rwi_voronoi(country['iso3'])

        print('---Convert settlement .tiff layer to .shp')
        convert_settlement_to_polygon(country['iso3'])

        print('---Processing intersect_coverage_areas_and_settlements')
        intersect_coverage_areas_and_settlements(country['iso3'])

        print('---Estimate customers per site')
        estimate_customers_per_site(country['iso3'])

        print('---Export site customer data to .csv')
        export_to_csv(country['iso3'])

    collect_cells(countries)

    concat_acled_results()


