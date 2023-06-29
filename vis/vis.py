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

    base = all_shapes.plot(column='bin', ax=ax, cmap='inferno_r', linewidth=0.1,
        legend=True, edgecolor='grey')
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


if __name__ == '__main__':


    plot_events()


    print('Complete')
