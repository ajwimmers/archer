
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


def archer4_visir(image, attrib, first_guess, para_fix=True, display_filename=None):

    # Function output dictionaries
    in_dict = {}
    out_dict = {}
    score_dict = {}


    # Channel-specific settings
    if attrib['archer_channel_type'] == 'IR':

        structure_height_km = 16
        mask_val = 265 # 265: Keep low cloud only
        ring_weight = 0.0167
        penalty_weight = 0.33
        image['bt_grid'] = image['data_grid'] # Archer acts on BT

    elif attrib['archer_channel_type'] == 'SWIR':

        structure_height_km = 16
        mask_val = -1e6 # 265: Keep low cloud only
        ring_weight = 0.0167
        penalty_weight = 0.33
        image['bt_grid'] = image['data_grid'] # Archer acts on BT

    elif attrib['archer_channel_type'] in ['Vis', 'DNB']:

        structure_height_km = 2
        mask_val = -1e6
        ring_weight = 0.0020
        penalty_weight = 0.33

        if attrib['scan_type'] == 'Geo' and attrib['archer_channel_type'] == 'Vis':

            # Normalize for solar zenith angle. This matters because a difference
            # in average brightness value from image to image does indeed throw off
            # the results.

            cos_solar_zenith_grid = nvtbx.cos_solar_zenith(
                image['lon_grid'], image['lat_grid'], first_guess['time'])

            # Sun must be >1.5 degrees above horizon for that conversion to make sense
            cos_solar_zenith_grid[cos_solar_zenith_grid < 0.0262] = np.NaN

            # Convert to pseudo brightness temp for Archer to act on
            bv_norm = image['data_grid'] / cos_solar_zenith_grid**0.5
            image['bt_grid'] = 350 - 0.75*bv_norm

            # Filter out nighttime cases
            num_pix = np.product(np.shape(image['bt_grid']))
            num_nan = np.sum(np.isnan(image['bt_grid']), axis=(0,1))
            if num_nan / num_pix > 0.3:
                print('Too dark for ARCHER. Exiting.')
                return in_dict, out_dict, score_dict

        else:
            # Assume polar Vis/DNB imagery is already normalized:
            image['bt_grid'] = 350 - 0.75*image['data_grid']


    # Parallax fix

    if para_fix:

        if 'zen_grid' in image.keys() and 'azm_grid' in image.keys():

            image['lon_pc_grid'], image['lat_pc_grid'] = nvtbx.parallax_fix_allnav(
                image['lon_grid'], image['lat_grid'], 
                image['zen_grid'], image['azm_grid'], structure_height_km)

        else:

            if attrib['scan_type'] is 'Geo':
                image['lon_pc_grid'], image['lat_pc_grid'] = nvtbx.parallax_fix_geo(
                    image['lon_grid'], image['lat_grid'], 
                    attrib['nadir_lon'], attrib['sensor'], structure_height_km)

            elif attrib['scan_type'] is 'Crosstrack':
                image['lon_pc_grid'], image['lat_pc_grid'] = nvtbx.parallax_fix_crosstrack(
                    image['lon_grid'], image['lat_grid'], 
                    attrib['sensor'], attrib['archer_channel_type'], 
                    structure_height_km)

            else: 
                print('attrib[''scan_type''] is not valid and cannot parallax fix. Exiting.')
                return in_dict, out_dict, score_dict


    else:

        image['lon_pc_grid'], image['lat_pc_grid'] = \
            image['lon_grid'], image['lat_grid']


    # Reduce resolution of the image data to prevent lags and artifacts
    image = nvtbx.reduce_res(image, step_km=4)


    # Write all the input data to a dictionary
    in_dict['sensor'] = attrib['archer_channel_type']
    in_dict['op_lon'] = first_guess['lon']
    in_dict['op_lat'] = first_guess['lat']
    in_dict['op_vmax'] = first_guess['vmax']
    in_dict['ring_weight'] = ring_weight # Just for display purposes


    # Iterate through feature-level and then surface-level center-fixes, keep the best:
    # BTW, I am not proud of this code construction. It stems from something that was 
    # not very straightforward in the original Matlab syntax.
    for level in ['feature', 'surface']:

        # Note here that the naming of "default" lat/lon changes between in_dict and image. 
        # image lat/lon grid is *unaltered* lat/lon. However, in_dict lat/lon is the lat/lon
        # *to be used in Archer*.

        if level is 'feature':

            print('Computing center-fix on full image...')
            in_dict['lon_mx'] = image['lon_pc_grid']
            in_dict['lat_mx'] = image['lat_pc_grid']
            in_dict['bt_mx'] = image['bt_grid']

        elif level is 'surface':

            print('Computing center-fix on cloud-masked image (w/o parallax fix) ...')
            in_dict['lon_mx'] = image['lon_grid']
            in_dict['lat_mx'] = image['lat_grid']
            in_dict['bt_mx'][in_dict['bt_mx'] < mask_val] = np.NaN


            # With IR imagery, we have to consider pixels with partial clouds/convection.
            # Dilate these features by 2 (NOTE: for ~10 km resolution):

            if archer_channel_type is 'IR':

                dilate_val = 2
                out_mx = in_dict['bt_mx']

                n_rows, n_cols = np.shape(in_dict['bt_mx'])

                for row_idx in range(n_rows):
                    for col_idx in range(n_cols):

                        box_min_row = np.max([0, row_idx - dilate_val])
                        box_min_col = np.max([0, col_idx - dilate_val])
                        box_max_row = np.min([n_rows-1, row_idx + dilate_val])
                        box_max_col = np.min([n_cols-1, col_idx + dilate_val])
                        box = in_dict['bt_mx'][box_min_row:box_max_row, box_min_col:box_max_col]

                        if np.any(np.isnan(box)):
                            out_mx[row_idx, col_idx] = np.NaN

                in_dict['bt_mx'] = out_mx


        # Calculate gridded score components
        score_dict = sftbx.combo_parts_calc_3_0(in_dict, penalty_weight=penalty_weight)
        uses_target = sftbx.quality_check(score_dict)
        out_dict['uses_target'] = uses_target


        if level is 'feature':
            print('TBD: Save weak version here')


        if uses_target:

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

            alpha = cotbx.confidence_to_alpha(
                confidence_score, attrib['archer_channel_type'], 0, in_dict['op_vmax'])


            # Represent the center fix uncertainty in terms of radius of 50%
            # confidence and radius of 95% confidence

            x_arr = np.arange(0, 10, 0.01)
            cdf_arr = 1 - (alpha * x_arr +1) * np.exp(-alpha * x_arr)
            nearest_idx = np.argmin(np.abs(cdf_arr - 0.50))
            rad_50 = x_arr[nearest_idx]
            nearest_idx = np.argmin(np.abs(cdf_arr - 0.95))
            rad_95 = x_arr[nearest_idx]


            # Calculate the probability of having detected an eye (n/a for geo)
            if attrib['archer_channel_type'] == 'IR':
               eye_prob_stat = confidence_score * ring_score;
               calib_stat_arr = [0, .1, .5,  1, 1.5,  2,  3,  4,  5,  6,  7,  8,  9, 10, 11, 12,  13,  14,  15]
               calib_perc_arr = [0,  1,  7, 15,  40, 47, 50, 60, 65, 70, 75, 80, 85, 95, 98, 99, 100, 100, 100]
               eye_prob = np.interp(eye_prob_stat, calib_stat_arr, calib_perc_arr)
            else:
               eye_prob = np.nan


            # Pack up the output values

            #out_dict['TC_atcf_dtg'] = TC_atcf_dtg
            out_dict['archer_channel_type'] = attrib['archer_channel_type']
            out_dict['tc_time'] = first_guess['time'] 
            out_dict['target_lon'] = lon_combo_max
            out_dict['target_lat'] = lat_combo_max
            out_dict['weak_target_lon'] = np.NaN
            out_dict['weak_target_lat'] = np.NaN
            out_dict['confidence_score'] = confidence_score
            out_dict['ring_radius_deg'] = ring_radius_deg
            out_dict['alpha_parameter'] = alpha
            out_dict['radius50percCertDeg'] = rad_50
            out_dict['radius95percCertDeg'] = rad_95


            # If it made it this far successfully, then no need to try the next option ('surface')
            break 

        else:

            alpha = -1


        """
        Find some way to convert this:
        if isempty(out)
            fprintf(1, 'No center-fix found with either method. \n\n');
        end
        """

    # Plot test fig
    if display_filename is not None:
        ditbx.plot_diag(image, attrib, in_dict, out_dict, score_dict, display_filename=display_filename)


    return in_dict, out_dict, score_dict

