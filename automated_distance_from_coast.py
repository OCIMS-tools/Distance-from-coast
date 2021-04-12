# -*- coding: utf-8 -*-

"""

Created on Wed Feb  5 14:09:05 2020



@author: cristinarusso

"""

#%% Importing neccessary libraries and modules

import shutil
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
import sys
from geojson import FeatureCollection
if not sys.warnoptions: # switching off warnings
    import warnings
    warnings.simplefilter("ignore")

#%% Setting Paths
pname ='/home/ocims_platform/current_auto/LACCE_split_current_data/AC_LACCE_daily.nc'
pname_dist='/home/ocims_platform/current_auto/distance_data/'

save_path = '/home/ocims_platform/ocims_tsite/static/main/'
save_fname ='distance_from_coast.json'
path_out ='/home/ocims_platform/current_auto/distance_data/' 

#%% Distance function

def distance(a, b):

   """
   Calculates the distance between two GPS points (decimal)
   @param a: 2-tuple of point A
   @param b: 2-tuple of point B
   @return: distance in nautical miles
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

def terminal_coordinates(starting_longitude,starting_latitude,distance,bearing):
    
    ''' This function returns the terminal geographical coordinates, 
        along a specific distance (in km) and at a specified bearing (in deg).

        N.B. This function is an adaptation of code found on stackoverflow:
             https://stackoverflow.com/questions/7222382/get-lat-long-given-current-point-distance-and-bearing

        INPUT:
            starting_latitude = latitude in degress of the starting coordinate pair
            starting_longitude = longitude in degress of the starting coordinate pair
            distance = desired distance between start and end location
            bearing = bearing between start and end location

        OUTPUT:
            lat_end = terminal latitude in degrees
            lon_end = terminal longitude in degrees

        '''
    import math
    
    R = 6378.1 #Radius (km) of the Earth
    brng = bearing
    brng = math.radians(brng) #Bearing is 90 degrees converted to radians.
    d = distance #Distance in km
    
    lat_start = math.radians(starting_latitude) #Current lat point converted to radians
    lon_start = math.radians(starting_longitude) #Current long point converted to radians
    
    lat_end = math.asin( math.sin(lat_start)*math.cos(d/R) +
         math.cos(lat_start)*math.sin(d/R)*math.cos(brng))
    
    lon_end = lon_start + math.atan2(math.sin(brng)*math.sin(d/R)*math.cos(lat_start),
                 math.cos(d/R)-math.sin(lat_start)*math.sin(lat_end))
    
    lat_end = math.degrees(lat_end)
    lon_end = math.degrees(lon_end)
    
    return lon_end,lat_end

#%% City Names

city_names = ['RB','DB','EL','PE']

city_start_coordinates = [[32.0383,-28.7807],
                          [31.0218,-29.8587],
                          [27.8546,-33.0292],
                          [25.6022,-33.9608]]

distance_from_cities = 300
bearing_from_cities = 135 #perpendicular to current

city_end_coordinates = []
for i in range(len(city_names)):
    city_end_coordinates.append(terminal_coordinates(city_start_coordinates[i][0],city_start_coordinates[i][1],distance_from_cities,bearing_from_cities))


#%% Importing data

ds        = xr.open_dataset(pname)
time      = ds.time
ac_lat    = ds.ac_lat.values
ac_lon    = ds.ac_lon.values
    
#%% RICHARDS BAY

for i in range(len(city_start_coordinates)):
    
    exec('ac_lat_'+city_names[i]+'= []')
    exec('ac_lon_'+city_names[i]+'= []')
    
    ## Defining line start and end coordinates

    p1 = Point(city_start_coordinates[i])
    p2 = Point(city_end_coordinates[i])
    
    ## Making line segment

    l = LineString([p1, p2])
    
    ## Make list with core coordinates

    h = []

    for j in range(len(ac_lat)):
        h.append((ac_lon[j],ac_lat[j]))
    
    ## Make line segments with core coordinates

    a = []
    for j in range(len(h)-1):
        a.append(LineString([h[j],h[j+1]]))
    
    
    ## Find where cities line intersects with core coordinates   
    p = []   
    for j, line in enumerate((a)):
        res = l.intersection(line)
        if np.size(res)> 0:
            p.append(res)
    
    ## Save intersection coordinates

    intersection_lat = []
    intersection_lon = []

    for j in range(len(p)):
        intersection_lon.append(p[j].coords.xy[0][0]) #longitude
        intersection_lat.append(p[j].coords.xy[1][0]) #latitude
    
    ## If there is no lat coordinate, means line doesn't intersect, append nan

    if np.size(intersection_lat)==0:
        exec('ac_lat_' +city_names[i] +'= np.nan')
    else:
        c = np.nanmin(intersection_lat)
        exec('ac_lat_' +city_names[i] +'= c')
    
    ## If there is no lon coordinate, means line doesn't intersect, append nan

    if np.size(intersection_lon)==0:
        exec('ac_lon_' +city_names[i] +'= np.nan')
    else:
        c = np.nanmin(intersection_lon)
        exec('ac_lon_' +city_names[i] +'= c')
    
#%% Calculate Distances

for i in range(len(city_names)):
    x = eval('ac_lon_'+city_names[i])
    y = eval('ac_lat_'+city_names[i])
    exec(city_names[i]+'_distance'+' = np.expand_dims(np.array(distance(city_start_coordinates[i],[x,y])),0)')
    print(eval(city_names[i]+'_distance'))
#%%

for i in range(len(city_names)):
    exec('ac_lon_'+city_names[i]+ ' = np.expand_dims('+'ac_lon_'+city_names[i]+',0)')
    exec('ac_lat_'+city_names[i]+ ' = np.expand_dims('+'ac_lat_'+city_names[i]+',0)')


# %% Creating array for JSON format

dist_ds=xr.open_dataset(pname_dist+'city_distances.nc')

RB_dist_array = dist_ds.distance_RB.values.tolist();RB_dist_array.insert(0,'distance')
DB_dist_array = dist_ds.distance_DB.values.tolist();DB_dist_array.insert(0,'distance')
EL_dist_array = dist_ds.distance_EL.values.tolist();EL_dist_array.insert(0,'distance')
PE_dist_array = dist_ds.distance_PE.values.tolist();PE_dist_array.insert(0,'distance')

dates = dist_ds.time.values
time=[]

for x in np.arange(len(dates)):
    ts = pd.to_datetime(str(dates[x]))
    d = ts.strftime('%Y-%m-%d')
    time.append(d)

time.insert(0,'x')

ac_section = [[ac_lon_RB[0],ac_lat_RB[0]],[ac_lon_DB[0],ac_lat_DB[0]],[ac_lon_EL[0],ac_lat_EL[0]],[ac_lon_PE[0],ac_lat_PE[0]]]
distance_array = [RB_dist_array,DB_dist_array,EL_dist_array,PE_dist_array]

# %% Looping through arrays to create diction in geojson format

x={}

for i in np.arange(len(city_start_coordinates)):
    city_cord = [city_start_coordinates[i][0],city_start_coordinates[i][1]];
    ac_cord = [ac_section[i][0],ac_section[i][1]];
    tmp_cord = str(city_cord)+','+ str(ac_cord)
    x[i]={
        "type": "Feature",
            "geometry": {
               "type": "LineString",
                  "coordinates":[city_cord,ac_cord]
                    },
            "properties": {
                  "distance": distance_array[i],
                  "time":time
                  }
    }


xy = []

for ip in np.arange(len(city_start_coordinates)):
            xy.append(x[ip])

geojson_data = FeatureCollection(xy)

# %% Saving geojson file
 
with open(save_path+save_fname, 'w', encoding='utf-8') as f:
        f.write(str(geojson_data))

# Import originalnetcdf file
fname_old = 'city_distances2.nc'
df=xr.open_dataset(path_out+fname_old)

# Sorting out the time dimension
# adding new date to variable with old dates
time_new =ds.time
time_old = df.time
time_together=np.concatenate((time_old,time_new),axis=0)
time_together = np.squeeze(time_together)

# Importing old transports
distance_RB_old = df.distance_RB
distance_DB_old = df.distance_DB
distance_EL_old = df.distance_EL
distance_PE_old = df.distance_PE


# Expanding dimension of old variables for concatenation
#tot_20 =np.expand_dims(tot_20,axis=0)
#tot_96 = np.expand_dims(tot_96,axis=0)
#tot_172 = np.expand_dims(tot_172,axis=0)
#tot_198 = np.expand_dims(tot_198,axis=0)
#tot_248 = np.expand_dims(tot_248,axis=0)
#tot_beng = np.expand_dims(tot_beng,axis=0)

# Expanding dimension of new variables for concatenation
RB_distance = np.concatenate((distance_RB_old, distance_RB_new),axis=0)
DB_distance = np.concatenate((distance_DB_old, distance_DB_new),axis=0)
EL_distance = np.concatenate((distance_EL_old, distance_EL_new),axis=0)
PE_distance = np.concatenate((distance_PE_old, distance_PE_new),axis=0)

# Creating NetCDF for timeseries
dist_dict = {'distance_RB':(['time'],RB_distance),'distance_DB':(['time'],DB_distance),'distance_EL':(['time'],EL_distance),'distance_PE':(['time'],PE_distance),
                  'time':time_together}
sec_data = xr.Dataset(dist_dict)
fname = 'city_distances.nc'
sec_data.to_netcdf(path_out+fname,  unlimited_dims={'time':True}, format='NETCDF4_CLASSIC',mode='w')

# Making copy of netcdf and saving into same directory
shutil.copy(path_out+fname, path_out+'city_distances2.nc')

# Making copy of netcdf and saving into same directory
shutil.copy(path_out+fname, path_out+'city_distances3.nc')

