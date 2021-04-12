# -*- coding: utf-8 -*-
"""
Created on Wed Feb  5 09:56:40 2020

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
from shapely.geometry.polygon import Polygon
import os
from scipy import interpolate as interp
from scipy.io.netcdf import netcdf_file
from matplotlib import pyplot
from datetime import datetime
import sys
if not sys.warnoptions: # switching off warnings
    import warnings
    warnings.simplefilter("ignore")

#%% Setting the paths
path = '/home/ocims_platform/current_auto/LACCE_data/'
filename = 'LACCE_daily.nc'
path_out = '/home/ocims_platform/current_auto/LACCE_split_current_data/'
filename_out='AC_'+filename

#%% Importing data
ds          = xr.open_dataset(path+filename)
core_lat    = ds.core_lat.values
core_lon    = ds.core_lon.values
date        = ds.time.values
retro_lon   = np.nanmin(core_lon)
retro_loc   = [ n for n,i in enumerate(core_lon) if i==retro_lon ][0]
ac_lat_c    = core_lat[0:retro_loc]
ac_lon_c    = core_lon[0:retro_loc]

#%% Creating daily ACCE netcdf file 
# Define some dummy data
lon         = np.ones(len(ac_lon_c))
lat         = np.ones(len(ac_lon_c))
ac_lon_c    = np.expand_dims(ac_lon_c,1)
ac_lat_c    = np.expand_dims(ac_lat_c,1)
time = np.expand_dims(date,0)


#%% Creating daily ACCE netcdf file
# Write out data to a new netCDF file with some attributes
#filename = netcdf_file(path_out+filename_out, 'w')
split_current_dict = {'ac_lat':(['lat','time'],ac_lat_c),'ac_lon':(['lon','time'],ac_lon_c),'time':time,'lat':lat,'lon':lon}
sec_data = xr.Dataset(split_current_dict)
sec_data.to_netcdf(path_out+filename_out,  unlimited_dims={'time':True}, format='NETCDF4_CLASSIC',mode='w')


