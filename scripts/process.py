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

from misc import get_countries, income_lut, parameters

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
    path_in = os.path.join(RESULTS, iso3, filename)
    data = pd.read_csv(path_in)
    space_time_ids = data['space_time_id'].unique()

    cost_lut = {
        'ran': 40000,
        'tower': 15000,
        'power': 5000,
        'backhaul': 30000,
        'labor': 1000,
    }

    output = []
    
    for space_time_id in space_time_ids:
        radios = []
        for idx, item in data.iterrows():
            if space_time_id == item['space_time_id']:
                GID_0 = item['GID_0']
                NAME_0 = item['NAME_0']
                GID_2 = item['GID_2']
                NAME_1 = item['NAME_1']
                year = item['year']
                month = item['month']
                day = item['day']
                country_iso3 = item['country']
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
                # latitude = item['latitude']
                # longitude = item['longitude']
                # coordinate = item['coordinate']
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
            'space_time_id': space_time_id,
            'GID_0': GID_0,
            'NAME_0': NAME_0,
            'GID_2': GID_2,
            'NAME_1': NAME_1,
            'year': year,
            'month': month,
            'day': day,
            'country': country_iso3,
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
            # 'latitude': latitude,
            # 'longitude': longitude,
            # 'coordinate': coordinate,
            'source': source,
            'radio': radio,
            'mcc': mcc,
            'net': net,
            'ran_usd': cost_lut['ran'],
            'tower_usd': cost_lut['tower'],
            'power_usd': cost_lut['power'],
            'backhaul_usd': cost_lut['backhaul'],
            'labor_usd': cost_lut['labor'],
            'total_usd': (
                cost_lut['ran'] + cost_lut['tower'] + 
                cost_lut['power'] + cost_lut['backhaul'] + 
                cost_lut['labor']  
            )
        })

    output = pd.DataFrame(output)
    
    filename = 'direct_site_costs_{}.csv'.format(country['iso3'])
    path_out = os.path.join(RESULTS, iso3, filename)
    output.to_csv(path_out, index=False)

    return


def indirect_cost_analysis(iso3, income_lut, parameters):
    """
    Estimate indirect costs. 

    """
    filename = 'customers_per_site.csv'
    folder = os.path.join(RESULTS, iso3)
    path = os.path.join(folder, filename)
    data = pd.read_csv(path)

    unique_sites = data['space_time_id'].unique()

    data = data.to_dict('records')

    country_income_lut = income_lut[iso3]

    output = []

    for unique_site in unique_sites:
        for item in data:

            if item['pop'] < 0:
                continue

            if unique_site == item['space_time_id']:

                for key, value in country_income_lut.items(): 
                    
                    lower, upper = key.split('_')
                    
                    if float(lower) < item['rwi'] < float(upper):

                        income = value

                        break

                wealth_usd = item['pop'] * income

                output.append({
                    'space_time_id': item['space_time_id'],
                    'pop': item['pop'],
                    'rwi': item['rwi'],
                    'error': item['error'],
                    'wealth_usd': int(wealth_usd),
                    'wealth_cfa': int(round(wealth_usd * parameters['usd_to_cfa'], 0)),#usd to west african franc 
                })
    
    output = pd.DataFrame(output)
    output['usd_to_cfa'] = parameters['usd_to_cfa']
    filename = 'indirect_by_site_long_{}.csv'.format(iso3)
    folder = os.path.join(RESULTS, iso3)
    path = os.path.join(folder, filename)
    output.to_csv(path, index=False)

    filename = 'indirect_by_site_{}.csv'.format(iso3)
    folder = os.path.join(RESULTS, iso3)
    path = os.path.join(folder, filename)
    output = output.groupby(['space_time_id'])[('pop', 'wealth_usd', 'wealth_cfa')].sum().reset_index()
    output['usd_to_cfa'] = parameters['usd_to_cfa']
    output.to_csv(path, index=False)

    return


def combined_data(iso3):
    """
    Combine data. 
    
    """
    filename = 'direct_site_costs_{}.csv'.format(iso3)
    path = os.path.join(RESULTS, iso3, filename)
    data = pd.read_csv(path)

    filename = 'indirect_by_site_{}.csv'.format(iso3)
    folder = os.path.join(RESULTS, iso3)
    path = os.path.join(folder, filename)
    indirect = pd.read_csv(path)

    data = data.merge(
        indirect, 
        left_on='space_time_id',
        right_on='space_time_id'
    )

    filename = 'combined_results_{}.csv'.format(iso3)
    folder = os.path.join(RESULTS, iso3)
    path = os.path.join(folder, filename)
    data.to_csv(path, index=False)

    return


def site_cost_benefit_metrics(iso3, parameters):
    """
    Generate cost benefit metrics.    
    
    """
    filename = 'combined_results_{}.csv'.format(iso3)
    folder = os.path.join(RESULTS, iso3)
    path = os.path.join(folder, filename)
    data = pd.read_csv(path)
    data = data.to_dict('records')

    output = []

    for item in data:

        direct_cost = item['total_usd']
        
        indirect_cost_usd = (
            (item['wealth_usd']/365) * 
            parameters['restoration_days'] * 
            (parameters['dependence_percent']/100)
        )

        value_at_risk = direct_cost + indirect_cost_usd

        protection_costs = (
            parameters['annual_protection_costs'] * 
            parameters['time_period'] *
            parameters['total_sites_{}'.format(iso3)] *
            (parameters['sites_to_protect'] / 100) /
            len(data)
        )

        output.append({
            'space_time_id': item['space_time_id'],
            'GID_0': item['GID_0'],
            'NAME_0': item['NAME_0'],
            'GID_2': item['GID_2'],
            'NAME_1': item['NAME_1'],
            'year': item['year'],
            'month': item['month'],
            'day': item['day'],
            'country': iso3,
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
            'source': item['source'],
            'radio': item['radio'],
            'mcc': item['mcc'],
            'net': item['net'],
            'ran_usd': item['ran_usd'],
            'tower_usd': item['tower_usd'],
            'power_usd': item['power_usd'],
            'backhaul_usd': item['backhaul_usd'],
            'labor_usd': item['labor_usd'],
            'total_usd': item['total_usd'],
            'pop': round(item['pop']),
            'wealth_usd': item['wealth_usd'],
            'wealth_cfa': item['wealth_cfa'],
            'usd_to_cfa': item['usd_to_cfa'], 
            'direct_cost': round(direct_cost),
            'indirect_cost_usd': round(indirect_cost_usd),
            'value_at_risk': round(value_at_risk),
            'protection_costs': round(protection_costs),
        })

    output = pd.DataFrame(output)

    output['benefit-cost-ratio'] = round(output['value_at_risk'] / output['protection_costs'],1)

    filename = 'cba_site_results_{}.csv'.format(iso3)
    folder = os.path.join(RESULTS, iso3)
    path = os.path.join(folder, filename)
    output.to_csv(path, index=False)

    return


def aggregate_country_results():
    """
    Aggregate cost benefit metrics by country. 

    """
    output = []

    filenames = [
        'cba_site_results_',
        # 'combined_results_'
    ]

    for short_filename in filenames:

        for iso3 in ['BFA','MLI','NER']:

            filename = short_filename + '{}.csv'.format(iso3)
            folder = os.path.join(RESULTS, iso3)
            path = os.path.join(folder, filename)
            data = pd.read_csv(path)
            data = data.to_dict('records')
            output = output + data

        output = pd.DataFrame(output)

        filename = '{}.csv'.format(short_filename[:-1])
        folder = os.path.join(RESULTS)
        path = os.path.join(folder, filename)
        output.to_csv(path, index=False)

    return


def country_cost_benefit_metrics():
    """
    Aggregate cost benefit metrics by country. 

    """
    output = []

    for iso3 in ['BFA','MLI','NER']:
            
        filename = 'cba_site_results_{}.csv'.format(iso3)
        folder = os.path.join(RESULTS, iso3)
        path = os.path.join(folder, filename)
        data = pd.read_csv(path)

        del data['space_time_id']
        data = data.sum()

        value_at_risk = data['value_at_risk']
        protection_costs = data['protection_costs']

        b_n_k_c = value_at_risk
        c_n_k_c = protection_costs
        cbr_n_k_c = value_at_risk / protection_costs

        output.append({
            'iso3': iso3,
            'b_n_k_c': b_n_k_c, 
            'c_n_k_c': c_n_k_c, 
            'cbr_n_k_c': cbr_n_k_c
        })

    output = pd.DataFrame(output)

    filename = 'cba_results.csv'
    folder = os.path.join(RESULTS)
    path = os.path.join(folder, filename)
    output.to_csv(path, index=False)

    return


if __name__ == "__main__":

    countries = get_countries()

    failures = []
    for idx, country in countries.iterrows():

        if not country['iso3'] in ['BFA', 'MLI', 'NER']:
           continue

        print('Working on {}'.format(country['iso3']))

        calc_direct_costs(country['iso3'])

        indirect_cost_analysis(country['iso3'], income_lut, parameters)

        combined_data(country['iso3'])

        site_cost_benefit_metrics(country['iso3'], parameters)

    aggregate_country_results()

    # country_cost_benefit_metrics()
