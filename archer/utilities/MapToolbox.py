# -*- coding: utf-8 -*-
"""
Created on Thu Feb 25 17:10:01 2016

@author: wimmers
"""

from __future__ import print_function
import numpy as np
import netCDF4
import os


def distance_deg(lat1, lon1, lat2, lon2):
    lat_dist_deg = lat1 - lat2
    avg_lat = (lat1 + lat2) / 2
    lon_dist_deg = (lon1 - lon2) * np.cos(np.pi/180 * avg_lat)
    return np.sqrt(lat_dist_deg**2 + lon_dist_deg**2)


def nan_around_coasts_89(lon_nopc_mx, lat_nopc_mx, bth_mx):

    # Load coast map
    this_dir = os.path.dirname(os.path.realpath(__file__))
    ncfile = netCDF4.Dataset(os.path.join(this_dir, '../etc/new_world.nc'))
    coasts_lon = ncfile.variables['coasts_lon'][:]
    coasts_lat = ncfile.variables['coasts_lat'][:]
    coasts_lon = coasts_lon[~np.isnan(coasts_lon)] # Remove nans
    coasts_lat = coasts_lat[~np.isnan(coasts_lat)]
    ncfile.close

    # Estimate effective field of view for swath points if not defined
    num_rows, num_cols = np.shape(bth_mx)
    center_col_num = int(np.round(num_cols / 2))
    nadir_row_dist = distance_deg(lat_nopc_mx[0, center_col_num], lon_nopc_mx[0, center_col_num], 
        lat_nopc_mx[1, center_col_num], lon_nopc_mx[1, center_col_num])

    fov_eff = nadir_row_dist/2 * (1.4+0.05) # Distance on the diagonal + fudge for parallelogram-ness (.05 to 0.1)

    # Identify the coastal lon/lat points in the swath domain
    min_lon = np.min(lon_nopc_mx[:])
    max_lon = np.max(lon_nopc_mx[:])
    min_lat = np.min(lat_nopc_mx[:])
    max_lat = np.max(lat_nopc_mx[:])

    is_coastal_domain1 = np.logical_and(coasts_lon > min_lon, coasts_lon < max_lon)
    is_coastal_domain2 = np.logical_and(coasts_lat > min_lat, coasts_lat < max_lat)
    is_coastal_domain = np.logical_and(is_coastal_domain1, is_coastal_domain1)

    coasts_lon = coasts_lon[is_coastal_domain]
    coasts_lat = coasts_lat[is_coastal_domain]

    # Calculate the swath mask as the points affected by the coasts
    swath_mask = np.full(np.shape(bth_mx), False, dtype='bool') # all False values

    for ri in range(0, num_rows):
        for ci in range(0, num_cols):

            if bth_mx[ri, ci] < 260:
                continue

            in_range = np.logical_and(np.abs(coasts_lat - lat_nopc_mx[ri,ci]) < fov_eff,
                np.abs(coasts_lon - lon_nopc_mx[ri,ci]) < 1.6*fov_eff)

            if np.sum(in_range) == 0:
                continue

            if np.min(distance_deg(lat_nopc_mx[ri, ci], lon_nopc_mx[ri, ci], coasts_lat, coasts_lon)) < fov_eff:
                swath_mask[ri, ci] = True

    cleared_swath = bth_mx + 0
    cleared_swath[swath_mask] = np.nan

    print('Total points masked for coastal boundaries:', np.sum(swath_mask, (0,1)))

    return cleared_swath







