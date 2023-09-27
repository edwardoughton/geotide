"""
Visualize data.

Written by Ed Oughton

June 2022

"""
import os
import configparser
import numpy as np
import pandas as pd
import geopandas as gpd
import matplotlib.pyplot as plt
# import seaborn as sns
import contextily as ctx
from pylab import * #is this needed
import imageio

CONFIG = configparser.ConfigParser()
CONFIG.read(os.path.join(os.path.dirname(__file__), '..', 'scripts', 'script_config.ini'))
BASE_PATH = CONFIG['file_locations']['base_path']

DATA_PROCESSED = os.path.join(BASE_PATH, 'processed')
VIS = os.path.join(BASE_PATH, '..', 'vis', 'figures')


def plot_events():
    """
    Plot ACLED events. 

    """
    # #Loading regional data by pop density geotype
    # path = os.path.join(DATA_PROCESSED, 'all_regional_data.csv')
    # data = pd.read_csv(path)

    all_shapes = pd.DataFrame()
    for iso3 in ['BFA', 'MLI', 'NER']:
        path = os.path.join(DATA_PROCESSED, iso3, 'acled', 'acled_{}.shp'.format(iso3))
        data = gpd.read_file(path, crs='epsg:4326')
        all_shapes = pd.concat([all_shapes, data])

    #convert single digit days/months to have leading zeros
    all_shapes['day'] = all_shapes['day'].astype(str)
    all_shapes['day'] = all_shapes['day'].str.zfill(2)
    all_shapes['month'] = all_shapes['month'].astype(str)
    all_shapes['month'] = all_shapes['month'].str.zfill(2)

    all_shapes['date'] = (
        all_shapes['year'].astype(str) + 
        all_shapes['month'].astype(str) +
        all_shapes['day'].astype(str)
    )
    all_shapes['date'] = all_shapes['date'].astype(int)

    country_shapes = pd.DataFrame()
    for iso3 in ['BFA', 'MLI', 'NER']:
        path = os.path.join(DATA_PROCESSED, iso3, 'national_outline.shp')
        data = gpd.read_file(path, crs='epsg:4326')
        country_shapes = pd.concat([country_shapes, data])

    bins = [0,
            20200101,20200701,
            20210101,20210701,
            20220101,20220701,
            20230101,
            1e9]
    labels = ['Pre-2020',
              '2020 Q1-Q2','2020 Q3-Q4',
              '2021 Q1-Q2','2021 Q3-Q4',
              '2022 Q1-Q2','2022 Q3-Q4',
              '2023 Q1-Q2'#,'2023 Q3-Q4'
            ]

    all_shapes['bin'] = pd.cut(
        all_shapes['date'], 
        bins=bins,
        labels=labels
    )

    plt.rcParams["font.family"] = "Times New Roman"
    fig, ax = plt.subplots(1, 1, figsize=(12,7))

    minx, miny, maxx, maxy = country_shapes.total_bounds
    ax.set_xlim(minx, maxx-6)
    ax.set_ylim(miny-1, maxy-6)

    plt.figure()

    base = all_shapes.plot(
        column='bin', 
        ax=ax, 
        cmap='inferno_r', 
        linewidth=0.1,
        legend=True, 
        edgecolor='grey'
        )
    country_shapes.plot(ax=base, facecolor="none", edgecolor='black', linewidth=0.75)

    handles, labels = ax.get_legend_handles_labels()

    fig.legend(handles[::-1], labels[::-1])

    ctx.add_basemap(ax, crs=all_shapes.crs, source=ctx.providers.CartoDB.Voyager)

    fig.suptitle(
        str('ACLED Conflict Events Affecting Telecommunications\nInfrastructure in Burkina Faso, Mali and Niger'), 
        fontsize=18, fontname='Times New Roman')

    fig.tight_layout()
    filename = 'acled.png'
    fig.savefig(os.path.join(VIS, filename))

    plt.close(fig)

    return

def plot_gif_maps():
    """
    Plot maps and create gif
    
    """
        # #Loading regional data by pop density geotype
    # path = os.path.join(DATA_PROCESSED, 'all_regional_data.csv')
    # data = pd.read_csv(path)

    all_shapes = pd.DataFrame()
    for iso3 in ['BFA', 'MLI', 'NER']:
        path = os.path.join(DATA_PROCESSED, iso3, 'acled', 'acled_{}.shp'.format(iso3))
        data = gpd.read_file(path, crs='epsg:4326')
        all_shapes = pd.concat([all_shapes, data])

    #convert single digit days/months to have leading zeros
    all_shapes['day'] = all_shapes['day'].astype(str)
    all_shapes['day'] = all_shapes['day'].str.zfill(2)
    all_shapes['month'] = all_shapes['month'].astype(str)
    all_shapes['month'] = all_shapes['month'].str.zfill(2)

    all_shapes['date'] = (
        all_shapes['year'].astype(str) + 
        all_shapes['month'].astype(str) +
        all_shapes['day'].astype(str)
    )
    all_shapes['date'] = all_shapes['date'].astype(int)

    country_shapes = pd.DataFrame()
    for iso3 in ['BFA', 'MLI', 'NER']:
        path = os.path.join(DATA_PROCESSED, iso3, 'national_outline.shp')
        data = gpd.read_file(path, crs='epsg:4326')
        country_shapes = pd.concat([country_shapes, data])

    unique_dates = all_shapes['date']

    for year in range(2018,2024):
        
        for month in range(1,13):

            if year == 2023 and month > 5:
                return

            year_month_day = "{}{}{}".format(year, str(month).zfill(2), 31)

            shapes = all_shapes[all_shapes['date'] <= int(year_month_day)]

            plot_map(country_shapes, shapes, year, month)

    return


def plot_map(country_shapes, shapes, date_year, date_month):
    """
    Plot map. 

    """
    bins = [0,
            20200101,20200701,
            20210101,20210701,
            20220101,20220701,
            20230101,
            1e9]
    labels = ['Pre-2020',
              '2020 Q1-Q2','2020 Q3-Q4',
              '2021 Q1-Q2','2021 Q3-Q4',
              '2022 Q1-Q2','2022 Q3-Q4',
              '2023 Q1-Q2'#,'2023 Q3-Q4'
            ]

    shapes['bin'] = pd.cut(
        shapes['date'], 
        bins=bins,
        labels=labels
    )
    plt.rcParams["font.family"] = "Times New Roman"
    fig, ax = plt.subplots(1, 1, figsize=(12,7))

    minx, miny, maxx, maxy = country_shapes.total_bounds
    ax.set_xlim(minx+1, maxx-6)
    ax.set_ylim(miny, maxy-6)

    plt.figure()

    base = shapes.plot(
        column='bin', 
        ax=ax, 
        cmap='viridis_r', 
        linewidth=0.1,
        legend=True, 
        edgecolor='grey'
        )
    country_shapes.plot(ax=base, facecolor="none", edgecolor='black', linewidth=0.75)

    handles, labels = ax.get_legend_handles_labels()

    fig.legend(handles[::-1], labels[::-1])

    ctx.add_basemap(ax, crs=shapes.crs, source=ctx.providers.CartoDB.Voyager)

    date = "{} {}".format(find_month(date_month), date_year)

    fig.suptitle(
        str('ACLED Conflict Events Affecting Telecommunications Infrastructure\nin Burkina Faso, Mali and Niger\n{}'.format(date)), 
        fontsize=18, fontname='Times New Roman')

    fig.tight_layout()
    filename = '{}{}.png'.format(date_year, str(date_month).zfill(2))
    folder = os.path.join(VIS, 'gif_maps')
    if not os.path.exists(folder):
        os.mkdir(folder)
    fig.savefig(os.path.join(folder, filename))

    plt.close(fig)


def find_month(date_month):
    """
    Given a monthly number, return the correct monthly string. 

    """
    if str(date_month) == "1":
        return "January"
    elif str(date_month) == "2":
        return "February"
    elif str(date_month) == "3":
        return "March"
    elif str(date_month) == "4":
        return "April"
    elif str(date_month) == "5":
        return "May"
    elif str(date_month) == "6":
        return "June"
    elif str(date_month) == "7":
        return "July"
    elif str(date_month) == "8":
        return "August"
    elif str(date_month) == "9":
        return "September"
    elif str(date_month) == "10":
        return "October"
    elif str(date_month) == "11":
        return "November"
    elif str(date_month) == "12":
        return "December"
    else:
        print(date_month)
        return "Did not recognize month"

    return


def generate_gif():
    """
    Stitch all maps together into gif. 

    """
    images = []

    folder = os.path.join(VIS, 'gif_maps')

    filenames = os.listdir(folder)

    for filename in filenames:

        path = os.path.join(folder, filename)

        images.append(imageio.imread(path))

    kargs = {'duration': .25}
    imageio.mimsave(os.path.join(VIS, 'gif.gif'), images, **kargs)

    return print('Generated .gif')


if __name__ == '__main__':

    plot_events()

    plot_gif_maps()

    generate_gif()

    print('Complete')
