"""
Generate demand metrics.

Written by Ed Oughton. 

June 2023. 

"""
import sys
import os
import configparser
import json
import pandas as pd
import geopandas as gpd
import pyproj
from shapely.ops import transform
from shapely.geometry import Point, box, LineString
import rasterio
from rasterio.mask import mask
from rasterio.features import shapes
from rasterstats import zonal_stats
import numpy as np
from shapely.ops import cascaded_union
from geovoronoi import voronoi_regions_from_coords, points_to_coords

# from tqdm import tqdm

from misc import get_countries#, process_country_shapes, process_regions, get_regions 

CONFIG = configparser.ConfigParser()
CONFIG.read(os.path.join(os.path.dirname(__file__),'..', 'scripts', 'script_config.ini'))
BASE_PATH = CONFIG['file_locations']['base_path']

DATA_RAW = os.path.join(BASE_PATH, 'raw')
DATA_PROCESSED = os.path.join(BASE_PATH, 'processed')
RESULTS = os.path.join(BASE_PATH, '..', 'results')


def generate_site_coverage_areas(iso3):
    """
    Generate site coverage areas events layer.  

    """
    filename = "acled_{}.shp".format(iso3)
    folder = os.path.join(DATA_PROCESSED, iso3, 'acled')
    path_in = os.path.join(folder, filename)
    data = gpd.read_file(path_in, crs='epsg:4326')
    data = data.to_crs(3857)
    data['geometry'] = data['geometry'].buffer(10000)
    data = data.to_crs(4326)

    filename = "coverage_polygons.shp"
    folder = os.path.join(DATA_PROCESSED, iso3, 'coverage_polygons')
    if not os.path.exists(folder):
        os.makedirs(folder)
    path_out = os.path.join(folder, filename)
    data.to_file(path_out, crs='epsg:4326') 

    return


def process_rwi_geometry(iso3):
    """
    Process rwi .csv layer to points.

    """
    filename = 'regions_{}_{}.shp'.format(1, iso3)
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

        filename_out = '{}.shp'.format(region['GID_1']) #each regional file is named using the gid id
        folder_out = os.path.join(BASE_PATH, 'processed', iso3 , 'rwi', 'points')
        if not os.path.exists(folder_out):
            os.makedirs(folder_out)
        path_out = os.path.join(folder_out, filename_out)

        output = gpd.GeoDataFrame.from_features(output, crs='epsg:4326')
        output.to_file(path_out,crs="epsg:4326")
    
    return


def create_rwi_voronoi(iso3):
    """
    Creates rwi voronoi polygons by region.

    """
    filename = 'regions_{}_{}.shp'.format(1, iso3)
    path_in = os.path.join(BASE_PATH, 'processed', iso3, 'regions', filename)
    regions = gpd.read_file(path_in, crs="epsg:4326")

    for idx, region in regions.iterrows():

        print('Working on vononoi polygons for: {}'.format(region['GID_1']))

        filename = '{}.shp'.format(region['GID_1'])
        path_in = os.path.join(BASE_PATH, 'processed', iso3, 'rwi', 'points', filename)
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

        filename = '{}.shp'.format(region['GID_1'])
        folder_out = os.path.join(BASE_PATH, 'processed', iso3, 'rwi', 'polygons')
        if not os.path.exists(folder_out):
            os.makedirs(folder_out)
        path_out = os.path.join(folder_out, filename)
        output.to_file(path_out, crs="EPSG:4326")

    return


def process_settlement_layer(iso3):
    """
    Clip the settlement layer to the chosen country boundary and place in
    desired country folder.

    """
    filename = 'regions_{}_{}.shp'.format(1, iso3)
    path_in = os.path.join(BASE_PATH, 'processed', iso3, 'regions', filename)
    regions = gpd.read_file(path_in, crs="epsg:4326")

    folder = os.path.join(BASE_PATH, '..', '..', 'data_raw', 'settlement_layer')
    path_settlements = os.path.join(folder, 'ppp_2020_1km_Aggregated.tif')

    settlements = rasterio.open(path_settlements, 'r+')
    settlements.nodata = 255
    settlements.crs = {"init": "epsg:4326"}

    for idx, region in regions.iterrows():

        folder_out = os.path.join(DATA_PROCESSED, iso3, 'settlements')
        if not os.path.exists(folder_out):
            os.makedirs(folder_out)
        shape_path = os.path.join(folder_out, '{}.tif'.format(region['GID_1']))

        if os.path.exists(shape_path):
            return print('Completed settlement layer processing')

        geo = gpd.GeoDataFrame()
        geo = gpd.GeoDataFrame({'geometry': region['geometry']}, index=[0])

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

    return print('Completed processing of settlement layer')


def convert_settlement_to_polygon(iso3):
    """
    Convert to .shp.
    
    """
    filename = 'regions_{}_{}.shp'.format(1, iso3)
    path_in = os.path.join(BASE_PATH, 'processed', iso3, 'regions', filename)
    regions = gpd.read_file(path_in, crs="epsg:4326")

    for idx, region in regions.iterrows():

        filename = "{}.tif".format(region['GID_1'])
        folder = os.path.join(DATA_PROCESSED, iso3, 'settlements')
        path_in = os.path.join(folder, filename)

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

        output = gpd.GeoDataFrame.from_features(output)
        filename = "{}.shp".format(region['GID_1'])
        folder = os.path.join(DATA_PROCESSED, iso3, 'settlements', 'shapes')
        if not os.path.exists(folder):
            os.makedirs(folder)
        path_out = os.path.join(folder, filename)
        output.to_file(path_out, crs='epsg:4326')

    return


def intersect_coverage_areas_and_settlements(iso3):
    """
    Intersect the coverage area and population polygons. 
    
    """
    filename = 'regions_{}_{}.shp'.format(1, iso3)
    path_in = os.path.join(BASE_PATH, 'processed', iso3, 'regions', filename)
    regions = gpd.read_file(path_in, crs="epsg:4326")#[:1]

    for idx, region in regions.iterrows():

        filename = "{}.shp".format(region['GID_1'])
        folder_in = os.path.join(DATA_PROCESSED, iso3, 'settlements', 'shapes')
        path_in = os.path.join(folder_in, filename)
        data_pop = gpd.read_file(path_in, crs='epsg:4326')#[:100]

        filename = '{}.shp'.format(region['GID_1'])
        folder_in = os.path.join(DATA_PROCESSED, iso3, 'rwi', 'polygons')
        path_in = os.path.join(folder_in, filename)
        data_rwi = gpd.read_file(path_in, crs='epsg:4326')#[:1]

        output = []

        for idx1, pop_tile in data_pop.iterrows():
            for idx2, rwi_poly in data_rwi.iterrows():
                geom = pop_tile['geometry'].centroid
                if geom.intersects(rwi_poly['geometry']):
                    output.append({
                        'geometry': pop_tile['geometry'],
                        'properties': {
                            'pop': pop_tile['raster_val'],
                            'rwi': rwi_poly['rwi'],
                            'error': rwi_poly['error'],
                        } 
                    })

        if len(output) == 0:
            continue

        output = gpd.GeoDataFrame.from_features(output)
        filename = '{}.shp'.format(region['GID_1'])
        folder_out = os.path.join(DATA_PROCESSED, iso3, 'estimated_income')
        if not os.path.exists(folder_out):
            os.mkdir(folder_out)
        path_out = os.path.join(folder_out, filename)
        # print(output)
        output.to_file(path_out, crs='epsg:4326')

    return


def estimate_customers_per_site(iso3):
    """

    """
    path_settlements = os.path.join(DATA_PROCESSED, iso3,
        'settlements.tif')

    filename = 'coverage_polygons.shp'
    folder = os.path.join(DATA_PROCESSED, iso3, 'coverage_polygons')
    path = os.path.join(folder, filename)
    coverage_areas = gpd.read_file(path)

    output = []

    for index, coverage_area in coverage_areas.iterrows():

        gid_id = coverage_area['GID_1']

        with rasterio.open(path_settlements) as src:

            affine = src.transform
            array = src.read(1)
            array[array <= 0] = 0

            population_summation = [d['sum'] for d in zonal_stats(
                coverage_area['geometry'],
                array,
                stats=['sum'],
                nodata=0,
                affine=affine)][0]

            if population_summation > 50000:
                population_summation = 50000

            output.append({
                'geometry': coverage_area['geometry'],
                'properties': {
                    'id': coverage_area['id'],
                    'country': coverage_area['country'],
                    'GID_1': coverage_area['GID_1'],
                    'NAME_1': coverage_area['NAME_1'],
                    'coordinate': coverage_area['coordinate'],
                    'population': round(population_summation),
                }
            })

    if len(output) == 0:
        return

    output = gpd.GeoDataFrame.from_features(output, crs='epsg:4326')
    filename = 'customers_per_site.shp'
    folder_out = os.path.join(DATA_PROCESSED, iso3, 'customers_per_site')
    if not os.path.exists(folder_out):
        os.makedirs(folder_out)
    path_out = os.path.join(folder_out, filename)
    output.to_file(path_out, crs='epsg:4326')

    output = output[['id','country','GID_1','NAME_1','coordinate', 'population']]
    filename = 'customers_per_site.csv'
    folder_out = os.path.join(DATA_PROCESSED, iso3, 'customers_per_site')
    if not os.path.exists(folder_out):
        os.makedirs(folder_out)
    path_out = os.path.join(folder_out, filename)
    output.to_csv(path_out, index=False)

    return


if __name__ == "__main__":

    countries = get_countries()

    failures = []
    for idx, country in countries.iterrows():

        if not country['iso3'] in ['MLI', 'NER']: #'BFA'
           continue

        print('Generate site coverage areas')
        generate_site_coverage_areas(country['iso3'])

        print('Processing rwi geometry')
        process_rwi_geometry(country['iso3'])

        print('Processing rwi voronoi polygons')
        create_rwi_voronoi(country['iso3'])

        print('Processing settlement layer')
        process_settlement_layer(country['iso3'])

        print('Convert settlement .tiff layer to .shp')
        convert_settlement_to_polygon(country['iso3'])

        print('Convert settlement .tiff layer to .shp')
        intersect_coverage_areas_and_settlements(country['iso3'])

        # print('Estimate customers per site')
        # estimate_customers_per_site(country['iso3'])

