
from __future__ import print_function
import numpy as np


def lookupNearest(x0, y0):
    xi = numpy.abs(x-x0).argmin()
    yi = numpy.abs(y-y0).argmin()
    return data[yi,xi]

def lookupNearest_arr(x0_arr, y0_arr, x_arr, y_arr, data):

    # x0_arr, y0_arr: Set of points for nearest neighbor interpolation
    # x_arr, y_arr: Coordinates of data
    # data: Data grid

    xi_arr = [numpy.abs(x_arr-x0).argmin() for x0 in x0_arr]
    yi_arr = [numpy.abs(y_arr-y0).argmin() for y0 in y0_arr]

    return data[yi_arr, xi_arr]

def calc_inten_85(out_dict, score_dict):

    rpd = np.pi/180
    deg_arr = range(0, 360.1, 5)

    lon_arr = np.flatten(score_dict['lon_grid1'][0,:])
    lat_arr = np.flatten(score_dict['lat_grid1'][:,0])
    lon_grid = score_dict['lon_grid1']
    lat_grid = score_dict['lat_grid1']
    data_grid = score_dict['data_grid1']

    # Normalized for 89-GHz sensors
    #?
    #?

    # Use either weak or normal center point
    if out_dict['center_lon'] is None:
        archer_lon = out_dict['weak_center_lon']
        archer_lat = out_dict['weak_center_lat']
    else:
        archer_lon = out_dict['center_lon']
        archer_lat = out_dict['center_lat']

    # Find the ring with the strongest eyewall
    max_circle_bt = 1e12
    best_eyewall_bt_arr = np.zeros(np.shape([72]))

    ring_radius_deg = out_dict['ring_radius_deg']
    for rad_circle in range(ring_radius_deg, ring_radius_deg+0.31, 0.05):

        lon_circle_arr = archer_lon + np.cos(rpd*deg_arr) * rad_circle / np.cos(rpd*archer_lat)
        lat_circle_arr = archer_lat + np.sin(rpd*deg_arr) * rad_circle

        circle_bt_arr = lookupNearest_arr(lon_circle_arr, lat_circle_arr, 
            lon_arr, lat_arr, data_grid)

        if np.nanmax(circle_bt_arr) < max_circle_bt:
            max_circle_bt = np.nanmax(circle_bt_arr)
            best_eye_rad_deg = rad_circle
            best_eyewall_bt_arr = circle_bt_arr


    # Find the highest BT *inside* that ring
    is_inside_ring = ((lon_grid - archer_lon)/np.cos(rpd*archer_lat))**2 + \
        (lat_grid - archer_lat)**2 <= best_eye_rad_deg**2
    max_eye_bt = np.nanmax(data_grid[is_inside_ring])

    # Calculate the scores
    eyewall_arr = best_eyewall_bt_arr
    eyewall_arr_nn = eyewall_arr[~np.isnan(eyewall_arr)]
    fraction_eyewall_nn = len(eyewall_arr_nn) / len(eyewall_arr)
    if fraction_eyewall_nn < 0.85:
        eye_score = None
        diff_score = None
    else:
        fraction_below_ceiling = np.sum(eyewall_arr_nn<232 and max_eye_bt-eyewall_arr_nn>10)
        all_below_255 = np.all(eyewall_arr_nn<255)
        if all_below_255:
            fraction_below_hotspot = np.sum(max_eye_bt-eyewall_arr_nn > 20) / len(eyewall_arr_nn)
        else:
            fraction_below_hotspot = 0
        eye_score = np.max([fraction_below_ceiling, fraction_below_hotspot])
        diff_score = max_eye_bt - np.max(eyewall_arr_nn)

    final_score = (eye_score > 0.85) * 15 + diff_score

    # This prevents ridiculous eyewall diameters
    if eye_score < 0.50 or np.isnan(eye_score):
        best_eye_rad_deg = None


    return ?


