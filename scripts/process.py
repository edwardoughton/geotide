"""
Process sites data.

Ed Oughton

February 2022

"""
import sys
import os
import configparser
import pandas as pd
import geopandas as gpd
import pyproj
from shapely.ops import transform
from shapely.geometry import shape, Point, mapping, LineString, MultiPolygon
import rasterio
import random

from misc import get_countries #, get_scenarios, get_regions

CONFIG = configparser.ConfigParser()
CONFIG.read(os.path.join(os.path.dirname(__file__),'..', 'scripts', 'script_config.ini'))
BASE_PATH = CONFIG['file_locations']['base_path']

DATA_RAW = os.path.join(BASE_PATH, 'raw')
DATA_PROCESSED = os.path.join(BASE_PATH, 'processed')
RESULTS = os.path.join(BASE_PATH, '..', 'results')


def calc_direct_costs(iso3):
    """
    Carry out cost analysis. 
    
    """
    filename = 'acled_{}.csv'.format(iso3)
    path_in = os.path.join(RESULTS, filename)
    data = pd.read_csv(path_in)
    coords = data['space_time_id'].unique()

    output = []
    
    for coord in coords:
        radios = []
        for idx, item in data.iterrows():
            if coord == item['space_time_id']:
                GID_0 = item['GID_0']
                NAME_0 = item['NAME_0']
                GID_2 = item['GID_2']
                NAME_1 = item['NAME_1']
                year = item['year']
                month = item['month']
                day = item['day']
                country = item['country']
                admin1 = item['admin1']
                location = item['location']
                infra_cate = item['infra_cate']
                infra = item['infra']
                sub_event_ = item['sub_event_']
                damage_inf = item['damage_inf']
                severity_i = item['severity_i']
                notes = item['notes']
                fatalities = item['fatalities']
                actor1 = item['actor1']
                assoc_acto = item['assoc_acto']
                inter1 = item['inter1']
                actor2 = item['actor2']
                assoc_ac_1 = item['assoc_ac_1']
                inter2 = item['inter2']
                latitude = item['latitude']
                longitude = item['longitude']
                coordinate = item['coordinate']
                source = item['source']
                radios.append(item['radio'])
                mcc = item['mcc']
                net = item['net']

        if 'LTE' in radios:
            radio = 'LTE'
        elif 'UMTS' in radios:
            radio = 'UMTS'
        else:
            radio = 'GSM'

        output.append({
            'GID_0': GID_0,
            'NAME_0': NAME_0,
            'GID_2': GID_2,
            'NAME_1': NAME_1,
            'year': year,
            'month': month,
            'day': day,
            'country': country,
            'admin1': admin1,
            'location': location,
            'infra_cate': infra_cate,
            'infra': infra,
            'sub_event_': sub_event_,
            'damage_inf': damage_inf,
            'severity_i': severity_i,
            'notes': notes,
            'fatalities': fatalities,
            'actor1': actor1,
            'assoc_acto': assoc_acto,
            'inter1': inter1,
            'actor2': actor2,
            'assoc_ac_1': assoc_ac_1,
            'inter2': inter2,
            'latitude': latitude,
            'longitude': longitude,
            'coordinate': coordinate,
            'source': source,
            'radio': radio,
            'mcc': mcc,
            'net': net,
        })

    output = pd.DataFrame(output)
    
    filename = 'site_count_info.csv'
    path_out = os.path.join(RESULTS, filename)
    output.to_csv(path_out, index=False)

    return


def indirect_cost_analysis(iso3):
    """
    Estimate indirect costs. 

    """

    filename = 'customers_per_site.csv'
    folder = os.path.join(RESULTS, iso3)
    path = os.path.join(folder, filename)
    data = pd.read_csv(path)

    unique_sites = data['space_time_id'].unique()

    data = data.to_dict('records')

    lut = {
        "-.8_-.6": 300,
        "-.6_-.4": 400,
        "-.4_-.2": 500,
        "-.2_0": 600,
        "0_.2": 700,
        ".2_.4": 800,
        ".4_.6": 900,
        ".6_.8": 1000,
        ".8_1": 1200,
    }

    output = []

    for unique_site in unique_sites:
        for item in data:

            if item['pop'] < 0:
                continue

            if unique_site == item['space_time_id']:

                for key, value in lut.items(): #needs double checking
                    
                    lower, upper = key.split('_')
                    
                    if float(lower) < item['rwi'] < float(upper):

                        income = value

                        break

                wealth = item['pop'] * income

                output.append({
                    'space_time_id': item['space_time_id'],
                    'pop': item['pop'],
                    'rwi': item['rwi'],
                    'error': item['error'],
                    'wealth': wealth #this looks wrong in the .csv, double check
                })
    
    output = pd.DataFrame(output)
    filename = 'indirect_by_site.csv'
    folder = os.path.join(RESULTS, iso3)
    path = os.path.join(folder, filename)
    
    output.to_csv(path, index=False)

    return


if __name__ == "__main__":

    countries = get_countries()

    failures = []
    for idx, country in countries.iterrows():

        if not country['iso3'] in ['BFA']:#, 'MLI', 'NER']:
           continue

        print('Working on {}'.format(country['iso3']))

        calc_direct_costs(country['iso3']))

        # indirect_cost_analysis(country['iso3'])

