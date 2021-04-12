# -*- coding: utf-8 -*-

"""

Created on Wed Feb  5 14:09:05 2020
@author: cristinarusso

"""

#%% Importing neccessary libraries and modules

import math
import numpy as np
import pylab as plt
import xarray as xr
import pandas as pd
from scipy.ndimage import gaussian_filter1d
from shapely.geometry import Point
from shapely.geometry import Polygon as SPolygon, Point, LineString
import os
from scipy import interpolate as interp
from scipy.io.netcdf import netcdf_file
from matplotlib import pyplot
from datetime import datetime
import pandas as pd 
import sys
from geojson import FeatureCollection
if not sys.warnoptions: # switching off warnings
    import warnings
    warnings.simplefilter("ignore")

#%% Setting Paths
    
    

pname ='/home/ocims_platform/current_auto/LACCE_split_current_data/AC_LACCE_daily.nc' 
pname_dist ='/home/ocims_platform/current_auto/distance_data/'

save_path = '/home/ocims_platform/ocims_tsite/static/main/'

save_fname = 'nearest_distance_from_coast.json'
save_fname_dist = 'city_distance.json'
#%% Distance function

def distance(a, b):
   """
   Calculates the distance between two GPS points (decimal)
   @param a: 2-tuple of point A
   @param b: 2-tuple of point B
   @return: distance in m
   """

   from math import radians, atan2, sqrt
  
   r = 6371000             # average earth radius in m
   dLat = radians(a[1]-b[1])
   dLon = radians(a[0]-b[0])
   x = np.sin(dLat/2) ** 2 + \
       np.cos(radians(a[0])) * np.cos(radians(b[0])) *\
       np.sin(dLon/2) ** 2
   y = 2 * atan2(sqrt(x), sqrt(1-x))
   d = r * y

   return d*0.000539957

def findNearest(ilon,ilat,lon,lat):
        """

        loni,lati,di = findNearest(ilon,ilat,lon,lat):
        find closest point within lat and lon data
        Input:
        ilon, ilat (input longitude and latitude to match)
        lon, lat = longitude and latitude in which to look
        lon and lat must have same dimension
        return:
        loni: longitude of closest point within lon
        lati: latitude of closest point within lat
        di: distance between ilon,ilat and loni, lati
        """

        d=[]
        for i, j in zip(lon,lat):
                d.append(distance([ilon,ilat],[i,j]))
        ind, = np.where(np.array(d)==min(d))
#       print(lon[ind],lat[ind])
        ind=ind[-1]
        return lon[ind], lat[ind],np.array(d)[ind]

#%% City Names and coordinates
city_names = ['RB','DB','EL','PE']

city_start_coordinates = [[32.0383,-28.7807],
                          [31.0218,-29.8587],
                          [27.8546,-33.0292],
                          [25.6022,-33.9608]]

#%% Importing data

ds          = xr.open_dataset(pname)
time        = ds.time
ac_lat    = ds.ac_lat.values
ac_lon    = ds.ac_lon.values

#%% Find Nearest
for i in range(len(city_names)):
    exec(city_names[i]+'_nearest_lon,'+city_names[i]+'_nearest_lat,'+city_names[i]+'_nearest_distance'+'= findNearest('+str(city_start_coordinates[i][0])+','+str(city_start_coordinates[i][1])+',ac_lon,ac_lat)')

# %%
for i in range(len(city_names)):
    exec(city_names[i]+'distance = np.expand_dims(np.array('+city_names[i]+'_nearest_distance),0)')
    exec('ac_lat_'+city_names[i]+ ' = np.expand_dims(np.array('+city_names[i]+'_nearest_lat),0)')
    exec('ac_lon_'+city_names[i]+ ' = np.expand_dims(np.array('+city_names[i]+'_nearest_lon),0)')

# %% creating arrays from json format

dist_ds=xr.open_dataset(pname_dist+'city_distances.nc')

RB_distance = dist_ds.distance_RB.values.tolist();RB_distance.insert(0,'distance')
DB_distance = dist_ds.distance_DB.values.tolist();DB_distance.insert(0,'distance')
EL_distance = dist_ds.distance_EL.values.tolist();EL_distance.insert(0,'distance')
PE_distance = dist_ds.distance_PE.values.tolist();PE_distance.insert(0,'distance')

dates = dist_ds.time.values
time=[]

for x in np.arange(len(dates)):
    ts = pd.to_datetime(str(dates[x]))
    d = ts.strftime('%Y-%m-%d')
    time.append(d)

time.insert(0,'x')

ac_section = [[ac_lon_RB[0],ac_lat_RB[0]],[ac_lon_DB[0],ac_lat_DB[0]],[ac_lon_EL[0],ac_lat_EL[0]],[ac_lon_PE[0],ac_lat_PE[0]]]
distance_array = [RB_distance,DB_distance,EL_distance,PE_distance]

# %% Looping through values to create diction in geojson format

x={}

for i in np.arange(len(city_start_coordinates)):
    city_cord = [city_start_coordinates[i][0],city_start_coordinates[i][1]];
    ac_cord = [np.float64(ac_section[i][0]),np.float64(ac_section[i][1])];
    tmp_cord = str(city_cord)+','+ str(ac_cord)
    x[i]={
        "type": "Feature",
            "geometry": {
               "type": "LineString",
                  "coordinates":[city_cord,ac_cord]
                    },
            "properties": {
                  "distance": distance_array[i],
                  "time": time
                  }
    }


xy = []
for ip in np.arange(len(city_start_coordinates)):
            xy.append(x[ip])

geojson_data = FeatureCollection(xy)
 
# %% Saving geojson file
with open(save_path+save_fname, 'w', encoding='utf-8') as f:
        f.write(str(geojson_data))




