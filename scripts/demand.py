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
from tqdm import tqdm

from misc import get_countries#, process_country_shapes, process_regions, get_regions 

CONFIG = configparser.ConfigParser()
CONFIG.read(os.path.join(os.path.dirname(__file__),'..', 'scripts', 'script_config.ini'))
BASE_PATH = CONFIG['file_locations']['base_path']

DATA_RAW = os.path.join(BASE_PATH, 'raw')
DATA_PROCESSED = os.path.join(BASE_PATH, 'processed')
RESULTS = os.path.join(BASE_PATH, '..', 'results')


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


def process_settlement_layer(iso3):
    """
    Clip the settlement layer to the chosen country boundary and place in
    desired country folder.

    """
    filename = 'regions_{}_{}.shp'.format(2, iso3)
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
        shape_path = os.path.join(folder_out, '{}.tif'.format(region['GID_2']))

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

    failures = []
    for country in tqdm(countries):

        if not country['iso3'] in ['BFA', 'MLI', 'NER']: #, #
           continue

        print('---Generate initial site coverage areas')
        generate_site_coverage_areas(country['iso3'])

        print('---Clean site coverage areas')
        clean_coverage_areas(country['iso3'])

        print('---Processing rwi geometry')
        process_rwi_geometry(country['iso3'])

        print('---Processing rwi voronoi polygons')
        create_rwi_voronoi(country['iso3'])

        print('---Processing settlement layer')
        process_settlement_layer(country['iso3'])

        print('---Convert settlement .tiff layer to .shp')
        convert_settlement_to_polygon(country['iso3'])

        print('---Processing intersect_coverage_areas_and_settlements')
        intersect_coverage_areas_and_settlements(country['iso3'])

        print('---Estimate customers per site')
        estimate_customers_per_site(country['iso3'])

        print('---Export site customer data to .csv')
        export_to_csv(country['iso3'])
