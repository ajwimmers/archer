
import matplotlib
matplotlib.use('agg') # AJW 16Aug2018
import matplotlib.pyplot as plt
import matplotlib.patches as patches
import matplotlib.gridspec as gridspec
from matplotlib import cm
from mpl_toolkits.axes_grid1 import make_axes_locatable

plt.close('all')

matplotlib.rcParams['axes.titlesize'] = 10
matplotlib.rcParams['axes.labelsize'] = 10
matplotlib.rcParams['legend.fontsize'] = 10
matplotlib.rcParams['font.size'] = 10

import numpy as np
import os, sys
import time
import netCDF4
from importlib import reload

import archer.utilities.DisplayToolbox as ditbx
import archer.utilities.MapToolbox as mptbx
import archer.utilities.ScoreFuncs as sftbx
import archer.utilities.Conversions as cotbx
#import archer.utilities.NetcdfToolbox as nctbx
import archer.utilities.NavToolbox as nvtbx


def archer3_mw(image, attrib, first_guess, para_fix=True, display_filename=None):

    # Function output dictionaries
    in_dict = {}
    out_dict = {}
    score_dict = {}


    # Channel-specific settings
    if attrib['archer_channel_type'] == '89GHz':

        if first_guess['vmax'] < 45:
            ring_weight = 0.005 # very very small
        elif first_guess['vmax'] >= 45 and first_guess['vmax'] < 84:
            ring_weight = 0.0694
        elif first_guess['vmax'] >= 84:
            ring_weight = 0.0263
        else:
            print('Fatal error: op_vmax must be a scalar')

        structure_height_km = 10
        mask_val = 245 # 245 avoids ice
        penalty_weight = 1.0

        # Archer acts on BT. This sends it coast-masked BT.
        image['bt_grid'] = mptbx.nan_around_coasts_89(
            image['lon_grid'], image['lat_grid'], image['data_grid']) 


    elif attrib['archer_channel_type'] == '37GHz': # Not as well validated

        if first_guess['vmax'] < 45:
            ring_weight = 0.005
        elif first_guess['vmax'] >= 45 and first_guess['vmax'] < 84:
            ring_weight = 0.0694
        elif first_guess['vmax'] >= 84:
            ring_weight = 0.0263
        else:
            print('Fatal error: op_vmax must be a scalar')

        structure_height_km = 3
        mask_val = 50 # Basically, no mask
        penalty_weight = 1.0

        # Archer acts on BT. This sends it original BT.
        image['bt_grid'] = image['data_grid']


    # Parallax fix

    if para_fix:

        if attrib['scan_type'] == 'Conical':
            image['lon_pc_grid'], image['lat_pc_grid'] = nvtbx.parallax_fix_conical(
                image['lon_grid'], image['lat_grid'], 
                attrib['sensor'], structure_height_km)

        elif attrib['scan_type'] == 'Crosstrack':
            image['lon_pc_grid'], image['lat_pc_grid'] = nvtbx.parallax_fix_crosstrack(
                image['lon_grid'], image['lat_grid'], 
                attrib['sensor'], attrib['archer_channel_type'], 
                structure_height_km)

    else:

        image['lon_pc_grid'], image['lat_pc_grid'] = \
            image['lon_grid'], image['lat_grid']



    print('Computing center-fix on whole image...')

    # Write all the input data to a dictionary

    # Note here that the naming of "default" lat/lon changes between in_dict and image. 
    # image lat/lon grid is *unaltered* lat/lon. However, in_dict lat/lon is the lat/lon
    # *to be used in Archer*.

    in_dict['sensor'] = attrib['archer_channel_type']
    in_dict['lon_mx'] = image['lon_pc_grid']
    in_dict['lat_mx'] = image['lat_pc_grid']
    in_dict['bt_mx'] = image['bt_grid']
    in_dict['time'] = first_guess['time']
    in_dict['op_lon'] = first_guess['lon']
    in_dict['op_lat'] = first_guess['lat']
    in_dict['op_vmax'] = first_guess['vmax']
    in_dict['ring_weight'] = ring_weight # Just for display purposes


    # Calculate gridded score components
    score_dict = sftbx.combo_parts_calc_3_0(in_dict, penalty_weight=penalty_weight)


    # Clean up the ring score to allow a combo center in cloud-masked areas

    score_dict['ring_score_grid'][np.isnan(score_dict['ring_score_grid'])] = 0


    # Calculate combo score, target point

    combo_grid = (score_dict['spiral_score_grid'] - score_dict['penalty_grid']) + \
        ring_weight * score_dict['ring_score_grid']
    score_dict['combo_score_grid'] = combo_grid

    combo_score = np.max(combo_grid)
    max_idx = np.argmax(combo_grid)
    i_combo_max, j_combo_max = sftbx.ind2sub(np.shape(combo_grid), max_idx)

    lon_combo_max = score_dict['lon_grid1'][i_combo_max, j_combo_max]
    lat_combo_max = score_dict['lat_grid1'][i_combo_max, j_combo_max]

    ring_radius_deg = score_dict['ring_radius_grid'][i_combo_max, j_combo_max]
    ring_score = score_dict['ring_score_grid'][i_combo_max, j_combo_max]
    score_by_radius_arr = np.squeeze(score_dict['ring_score_grid_full'][i_combo_max, j_combo_max, :])
    gradient_grid = np.squeeze(score_dict['radial_gradient_4d'][i_combo_max, j_combo_max, :])


    # Calculate confidence score using combo grid w/o the distance penalty.
    # The confidence score is the maximum value minus the highest score CONFIDENCE_DIST_DEG away

    confidence_grid = (score_dict['spiral_score_grid'] - 0*score_dict['penalty_grid']) + \
        ring_weight * score_dict['ring_score_grid']
    confidence_max = np.nanmax(confidence_grid)
    confidence_max_idx = np.argmax(confidence_grid)
    i_conf_max, j_conf_max = sftbx.ind2sub(np.shape(confidence_grid), confidence_max_idx)

    lon_conf_max = score_dict['lon_grid1'][i_conf_max, j_conf_max]
    lat_conf_max = score_dict['lat_grid1'][i_conf_max, j_conf_max]

    dist_squared_grid = (score_dict['lat_grid1'] - lat_conf_max) ** 2 + \
        (np.cos(np.pi/180 * lat_conf_max) * (score_dict['lon_grid1'] - lon_conf_max)) ** 2 

    CONFIDENCE_DIST_DEG = 0.75
    confidence_score = confidence_max - \
        np.nanmax(confidence_grid[dist_squared_grid > CONFIDENCE_DIST_DEG ** 2])

    alpha = cotbx.confidence_to_alpha(confidence_score, attrib['archer_channel_type'], 0, in_dict['op_vmax'])


    # Represent the center fix uncertainty in terms of radius of 50%
    # confidence and radius of 95% confidence

    x_arr = np.arange(0, 10, 0.01)
    cdf_arr = 1 - (alpha * x_arr +1) * np.exp(-alpha * x_arr)
    nearest_idx = np.argmin(np.abs(cdf_arr - 0.50))
    rad_50 = x_arr[nearest_idx]
    nearest_idx = np.argmin(np.abs(cdf_arr - 0.95))
    rad_95 = x_arr[nearest_idx]


    # Calculate the probability of having detected an eye

    if attrib['archer_channel_type'] == '89GHz':
       eye_prob_stat = confidence_score * ring_score;
       calib_stat_arr = [0, 5, 10, 15, 20,  30,  40,  50,  60,  70,  75,  80]
       calib_perc_arr = [0, 9, 27, 44, 56,  72,  82,  90,  94,  99, 100, 100]
       eye_prob = np.interp(eye_prob_stat, calib_stat_arr, calib_perc_arr)
    else:
       eye_prob = np.nan


    # Pack up the output values

    uses_target = sftbx.quality_check(score_dict)
    out_dict['uses_target'] = uses_target

    if uses_target:
        # This is an official center-fix
        out_dict['center_lon'] = lon_combo_max
        out_dict['center_lat'] = lat_combo_max
        out_dict['weak_center_lon'] = None
        out_dict['weak_center_lat'] = None

    else:
        # This is a center-fix if you must, but it's not official because it's probably corrupted
        out_dict['center_lon'] = None
        out_dict['center_lat'] = None
        out_dict['weak_center_lon'] = lon_combo_max
        out_dict['weak_center_lat'] = lat_combo_max

    out_dict['archer_channel_type'] = attrib['archer_channel_type']
    out_dict['confidence_score'] = confidence_score
    out_dict['ring_radius_deg'] = ring_radius_deg
    out_dict['alpha_parameter'] = alpha
    out_dict['radius50percCertDeg'] = rad_50
    out_dict['radius95percCertDeg'] = rad_95


    # Plot test fig
    if display_filename is not None:
        ditbx.plot_diag(image, attrib, in_dict, out_dict, score_dict, display_filename=display_filename)

 
    return in_dict, out_dict, score_dict

