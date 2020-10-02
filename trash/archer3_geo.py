
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


def archer3_geo(image, attrib, first_guess, display_filename=None):

    # Temporarily convert to the following variable names:
    bt_mx = image['bt_grid']
    lon_nopc_mx = image['lon_grid']
    lat_nopc_mx = image['lat_grid']
    #time_arr = image['time_arr'] # Currently unused

    sensor = attrib['sensor']
    archer_channel_type = attrib['archer_channel_type']

    op_lon = first_guess['lon']
    op_lat = first_guess['lat']
    op_vmax = first_guess['vmax']
    op_time = first_guess['time']

    # Outputs
    in_dict = {}
    out_dict = {}
    score_dict = {}

    # Channel-specific settings
    if archer_channel_type == 'IR':

        structure_height_km = 16
        mask_val = 265 # 265: Keep low cloud only
        ring_weight = 0.0167
        penalty_weight = 0.33
        color_lo = [180, 180, 180]
        color_hi = [300, 300, 300]
        color_label = 'Brightness temp, K'

    elif archer_channel_type == 'SWIR':

        structure_height_km = 16
        mask_val = -1e6 # 265: Keep low cloud only
        ring_weight = 0.0167
        penalty_weight = 0.33
        color_lo = [180, 180, 180]
        color_hi = [300, 300, 300]
        color_label = 'Brightness temp, K'

    elif archer_channel_type == 'Vis':

        structure_height_km = 2
        mask_val = -1e6
        ring_weight = 0.0020
        penalty_weight = 0.33
        color_lo = [0,   150, 150]
        color_hi = [255, 320, 320]
        color_label = 'Pseudo brightness temp, K'

        # Normalize for solar zenith angle. This matters because a difference
        # in average brightness value from image to image does indeed throw off
        # the results.

        cos_solar_zenith_grid = nvtbx.cos_solar_zenith(
            lon_nopc_mx, lat_nopc_mx, op_time)

        # Sun must be >1.5 degrees above horizon
        cos_solar_zenith_grid[cos_solar_zenith_grid < 0.0262] = np.NaN

        bv_norm = bt_mx / cos_solar_zenith_grid**0.5
        bt_mx = 350 - 0.75*bv_norm

        # Filter out nighttime cases
        num_pix = np.product(np.shape(bt_mx))
        num_nan = np.sum(np.isnan(bt_mx), axis=(0,1))
        if num_nan / num_pix > 0.3:
            print('Looks too dark. Exiting.')
            return in_dict, out_dict, score_dict


    # Parallax fix
    lon_mx, lat_mx = nvtbx.parallax_fix_geo(
        lon_nopc_mx, lat_nopc_mx, attrib['nadir_lon'], sensor, structure_height_km)

    # Write all the input data to a dictionary
    in_dict['sensor'] = archer_channel_type
    in_dict['op_lon'] = op_lon
    in_dict['op_lat'] = op_lat
    in_dict['op_vmax'] = op_vmax


    # Iterate through feature-level and then surface-level center-fixes, keep the best:
    # BTW, I am not proud of this code construction. It stems from something that was 
    # not very straightforward in the original Matlab syntax.
    for level in ['feature', 'surface']:

        if level is 'feature':

            print('Computing center-fix on parallax-corrected image...')
            in_dict['lon_mx'] = lon_mx
            in_dict['lat_mx'] = lat_mx
            in_dict['bt_mx'] = bt_mx

        elif level is 'surface':

            print('Computing  center-fix on cloud-masked image...')
            in_dict['lon_mx'] = lon_nopc_mx
            in_dict['lat_mx'] = lat_nopc_mx
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

            alpha = cotbx.confidence_to_alpha(confidence_score, archer_channel_type, 0, in_dict['op_vmax'])


            # Represent the center fix uncertainty in terms of radius of 50%
            # confidence and radius of 95% confidence

            x_arr = np.arange(0, 10, 0.01)
            cdf_arr = 1 - (alpha * x_arr +1) * np.exp(-alpha * x_arr)
            nearest_idx = np.argmin(np.abs(cdf_arr - 0.50))
            rad_50 = x_arr[nearest_idx]
            nearest_idx = np.argmin(np.abs(cdf_arr - 0.95))
            rad_95 = x_arr[nearest_idx]


            # Calculate the probability of having detected an eye (n/a for geo)
            if archer_channel_type == 'IR':
               eye_prob_stat = confidence_score * ring_score;
               calib_stat_arr = [0, .1, .5,  1, 1.5,  2,  3,  4,  5,  6,  7,  8,  9, 10, 11, 12,  13,  14,  15]
               calib_perc_arr = [0,  1,  7, 15,  40, 47, 50, 60, 65, 70, 75, 80, 85, 95, 98, 99, 100, 100, 100]
               eye_prob = np.interp(eye_prob_stat, calib_stat_arr, calib_perc_arr)
            else:
               eye_prob = np.nan


            # Pack up the output values

            #out_dict['TC_atcf_dtg'] = TC_atcf_dtg
            out_dict['archer_channel_type'] = archer_channel_type
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

        if archer_channel_type is 'Vis':
            disp_grid = image['bt_grid'] # Varies from 0 to 255
            #sub0_lo = 0
            #sub0_hi = 255
        else:
            disp_grid = bt_mx # Varies from 180 to 300 usually
            #sub0_lo = color_lo
            #sub0_hi = color_hi


        # Turn interactive plotting off
        plt.ioff()

        fig = plt.figure()

        # Fix for pcolor offsets
        lon_mx_disp, lat_mx_disp, bt_mx_disp = ditbx.pcolorCenterShift(lon_mx, lat_mx, disp_grid)

        gs = gridspec.GridSpec(1, 3, width_ratios=(1, 1, 1))

        ax0 = plt.subplot(gs[0,0])
        # This explains the nan-related runtime error for the following line: https://github.com/matplotlib/matplotlib/issues/9892
        img0 = plt.pcolormesh(lon_mx_disp, lat_mx_disp, bt_mx_disp, vmin=color_lo[0], vmax=color_hi[0], cmap='jet_r')
        #img0 = plt.scatter(lon_mx, lat_mx, c=bt_noco_mx, vmin=160, vmax=280, cmap='jet_r', edgecolors='none', marker=',')
        spiral_score_grid_disp = score_dict['spiral_score_grid'] + 0
        spiral_score_grid_disp[spiral_score_grid_disp < -1e8] = np.nan
        plt.contour(score_dict['lon_grid1'], score_dict['lat_grid1'], spiral_score_grid_disp, 
            np.arange(0,50), colors='gray')
        plt.plot(op_lon, op_lat, 'w+')
        plt.title('Spiral score contours, alpha = %4.1f' % (alpha))

        aspect_lon = np.cos(np.pi/180*op_lat) # Aspect ratio
        plt.xlim(op_lon-2.5/aspect_lon, op_lon+2.5/aspect_lon)
        plt.ylim(op_lat-2.5, op_lat+2.5)

        # Colorbar, h/t https://joseph-long.com/writing/colorbars/
        divider = make_axes_locatable(ax0)
        cax0 = divider.append_axes("bottom", size="5%", pad=0.30)
        fig.colorbar(img0, cax=cax0, orientation='horizontal')

        # The next two are plotted on the product grid:
        # Fix for pcolor offsets
        lon_grid_disp, lat_grid_disp, bt_grid_disp = ditbx.pcolorCenterShift(
            score_dict['lon_grid1'], score_dict['lat_grid1'], score_dict['data_grid1'] )

        ax1 = plt.subplot(gs[0,1])
        # This explains the nan-related runtime error for the following line: https://github.com/matplotlib/matplotlib/issues/9892
        img1 = plt.pcolormesh(lon_grid_disp, lat_grid_disp, bt_grid_disp, vmin=color_lo[1], vmax=color_hi[1], cmap='jet_r')

        ring_score_grid_disp = score_dict['ring_score_grid'] + 0
        ring_score_grid_disp[ring_score_grid_disp == 0] = np.nan
        plt.contour(score_dict['lon_grid1'], score_dict['lat_grid1'], ring_weight * ring_score_grid_disp, 
            np.arange(0,50), colors='gray')
        #plt.plot(op_lon, op_lat, 'w+')
        plt.title(('%s  Vmax = %d kt\nRing score contours') % 
            (archer_channel_type, op_vmax) )

        ring_score_max_idx = np.argmax(score_dict['ring_score_grid'])
        i_rs_max, j_rs_max = sftbx.ind2sub(np.shape(score_dict['ring_score_grid']), ring_score_max_idx)
        rs_max_lon = score_dict['lon_grid1'][i_rs_max, j_rs_max]
        rs_max_lat = score_dict['lat_grid1'][i_rs_max, j_rs_max]
        plt.plot(rs_max_lon, rs_max_lat, 'ws', markerfacecolor='none')


        # Colorbar, h/t https://joseph-long.com/writing/colorbars/
        divider = make_axes_locatable(ax1)
        cax1 = divider.append_axes("bottom", size="5%", pad=0.30)
        cbar = fig.colorbar(img1, cax=cax1, orientation='horizontal')
        cbar.set_label(color_label)


        ax2 = plt.subplot(gs[0,2])
        # This explains the nan-related runtime error for the following line: https://github.com/matplotlib/matplotlib/issues/9892
        img2 = plt.pcolormesh(lon_grid_disp, lat_grid_disp, bt_grid_disp, vmin=color_lo[2], vmax=color_hi[2], cmap='jet_r')

        if uses_target:
            combo_score_grid_disp = combo_grid + 0
            combo_score_grid_disp[combo_score_grid_disp < -1e8] = np.nan
            plt.contour(score_dict['lon_grid1'], score_dict['lat_grid1'], combo_score_grid_disp, 
                np.arange(0,50), colors='gray')
            plt.plot(op_lon, op_lat, 'w+')
            plt.plot(lon_combo_max, lat_combo_max, 'ws', markerfacecolor='none')

            if not np.isnan(out_dict['target_lon']):
                ring_ang_arr = np.arange(0, 361, 5)
                lon_scale = np.cos(out_dict['target_lat'] * np.pi/180)
                ring_lon_arr = out_dict['target_lon'] + np.cos(ring_ang_arr * np.pi/180) * ring_radius_deg / lon_scale
                ring_lat_arr = out_dict['target_lat'] + np.sin(ring_ang_arr * np.pi/180) * ring_radius_deg 
                plt.plot(ring_lon_arr, ring_lat_arr, 'm-')

        plt.title('Combined score contours')

        # Colorbar, h/t https://joseph-long.com/writing/colorbars/
        divider = make_axes_locatable(ax2)
        cax2 = divider.append_axes("bottom", size="5%", pad=0.30)
        fig.colorbar(img2, cax=cax2, orientation='horizontal')

        #print(np.round(ring_weight * score_dict['ring_score_grid'] * 100)[::4, ::4].astype('int'))

        # Colorbar, h/t https://joseph-long.com/writing/colorbars/
        divider = make_axes_locatable(ax2)
        cax2 = divider.append_axes("bottom", size="5%", pad=0.30)
        fig.colorbar(img2, cax=cax2, orientation='horizontal')


        fig.set_size_inches(15, 5)
        #plt.show()
        plt.savefig(display_filename, dpi=100)


    return in_dict, out_dict, score_dict

