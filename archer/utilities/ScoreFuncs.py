from __future__ import print_function
import numpy as np
import archer.utilities.InterpToolbox as intbx
#from importlib import reload
#reload(intbx)


def ind2sub(array_shape, ind):
    rows = (ind.astype('int') // array_shape[1])
    cols = (ind.astype('int') % array_shape[1]) # or numpy.mod(ind.astype('int'), array_shape[1])
    return (rows, cols)


def distance_deg(lon1, lat1, lon2, lat2):
    lat_dist_deg = lat1 - lat2
    avg_lat = np.mean(lat2)
    lon_dist_deg = (lon1 - lon2) * np.cos(np.pi/180 * avg_lat)
    return np.sqrt(lat_dist_deg**2 + lon_dist_deg**2)


def combo_parts_calc_3_0(mi_dict, penalty_weight=1.0):

    #Functions:
    #spiralCenterLowRes.m - Find the spiral center (large range)
    #spiralCenterLowRes2nd.m - Hone in on the center and get a good score (short range)
    #ringFitScoresOverhaul.m - Assign scores according to the best ring fit
    #distanceDeg.m - Find the distance between two points (in g.c.d.)

    # Resampling, spiral and ring parameters
    if mi_dict['sensor'] == '89GHz' or mi_dict['sensor'] == '37GHz' or mi_dict['sensor'] == '183GHz':
        lon_inc = 0.05
        lat_inc = 0.05
        perim_deg = 1.6
        filter_radius_deg = 2.0
        spiral_search_radius_deg = 2.0
        max_radius_deg = 2.50
    else:
        lon_inc = 0.025
        lat_inc = 0.025
        perim_deg = 2.5
        filter_radius_deg = 2.5
        spiral_search_radius_deg = 2.0
        max_radius_deg = 0.50

    # Universal spiral parameters
    spiral_weight = 15
    spiral_offset = 20
    spiral_spacing_deg = 0.05

    # Universal ring parameters
    ring_weight = 250.0
    min_radius_deg = 0.05

    # Remove any false data
    if mi_dict['sensor'] in ['89GHz', '37GHz',  '183GHz', 'IR']:
        with np.errstate(invalid='ignore'):
            mi_dict['bt_mx'][(mi_dict['bt_mx'] < 80)] = np.nan

    # Resample the swath to a regular grid, centred on the fx point, and with 
    # the correct aspect ratio at the fx point
    print('Regridding data to equal lon/lat aspect ratio ... ')
    #lon_mx_offset = (mi_dict['lon_mx'] - mi_dict['op_lon']) * np.cos(np.pi/180 * mi_dict['op_lat'])
    #lat_mx_offset = mi_dict['lat_mx'] - mi_dict['op_lat']
    #good_points = np.logical_and(~np.isnan(lon_mx_offset), ~np.isnan(lat_mx_offset))
    #lon_mx_offset_nn = lon_mx_offset[good_points]
    #lat_mx_offset_nn = lat_mx_offset[good_points]
    #bth_mx_nn = mi_dict['bth_mx'][good_points]
    #print('np.shape(lon_mx_offset_nn) = ', np.shape(lon_mx_offset_nn))

    # Represent data as a regular grid centered at (0,0). This is not the prettiest approach, but
    # it stays faithful to the legacy code.
    x_arr_offset_gcd = np.arange(-perim_deg, perim_deg+1e-6, lon_inc)
    y_arr_offset_gcd = np.arange(perim_deg, -perim_deg-1e-6, -lat_inc)
    x_grid_offset_gcd, y_grid_offset_gcd = np.meshgrid(x_arr_offset_gcd, y_arr_offset_gcd)
    lon_arr1 = x_arr_offset_gcd / np.cos(np.pi/180 * mi_dict['op_lat']) + mi_dict['op_lon']
    lat_arr1 = y_arr_offset_gcd + mi_dict['op_lat']
    lon_grid1, lat_grid1 = np.meshgrid(lon_arr1, lat_arr1)
    """
    data_grid1 = intbx.InterpSectionToGlobalRectGrid(
        lon_mx_offset, lat_mx_offset, mi_dict['bt_mx'], x_arr_offset_gcd, y_arr_offset_gcd,
        interp_type='linear') # Doesn't work somehow
    """
    from scipy.interpolate import griddata
    #gridCoordinates = list(zip(lon_mx_offset.ravel(), lat_mx_offset.ravel()))
    #data_grid1 = griddata(gridCoordinates, mi_dict['bt_mx'].ravel(), 
    #    (x_grid_offset_gcd, y_grid_offset_gcd), 'linear')

    lon_mx_buff, lat_mx_buff, bt_mx_buff = intbx.addEdgeBuffer(
        mi_dict['lon_mx'], mi_dict['lat_mx'], mi_dict['bt_mx'], 'linear')
    #lon_mx_buff, lat_mx_buff, bt_mx_buff = (mi_dict['lon_mx'], mi_dict['lat_mx'], mi_dict['bt_mx'])
    gridCoordinates = list(zip(lon_mx_buff.ravel(), lat_mx_buff.ravel()))
    data_grid1 = griddata(gridCoordinates, bt_mx_buff.ravel(), (lon_grid1, lat_grid1), 'linear')
    #print('np.shape(data_grid1) = ', np.shape(data_grid1))
    #print('data_grid1 = ', data_grid1)

    # Spiral center
    print('Calculating spiral center ...')
    sp_grid1, fraction_input = spiral_center_calc(x_grid_offset_gcd, y_grid_offset_gcd, data_grid1, mi_dict['sensor'],
        mi_dict['op_lon'], mi_dict['op_lat'], 
        filter_radius_deg, spiral_search_radius_deg, spiral_spacing_deg)
    spiral_score_grid = spiral_weight * sp_grid1 - spiral_offset

    # Ring center
    print('Calculating ring center ... ')

    # Add in the penalty for distance from first guess
    penalty_grid = penalty_weight * distance_deg(mi_dict['op_lon'], mi_dict['op_lat'], lon_grid1, lat_grid1)
    spiral_score_grid_with_penalty = spiral_score_grid - penalty_grid

    # Find the swarm of points to test for ring fitting
    spiral_fit_buffer = 1.5 # Include all points within this range of the max value ...
    swarm_reach = 0.25 # ... plus this distance out

    spiral_score_grid_with_penalty_nonan = spiral_score_grid_with_penalty
    spiral_score_grid_with_penalty_nonan[np.isnan(spiral_score_grid_with_penalty)] = -1e9
    is_inside_buffer = spiral_score_grid_with_penalty_nonan > \
        np.nanmax(spiral_score_grid_with_penalty) - spiral_fit_buffer

    #print('np.sum(is_inside_buffer) = ', np.sum(is_inside_buffer))
    #print(is_inside_buffer[::8, ::8])
   

    # Expand the buffer by swarm_reach:
    # I had to change this on from the matlab version, which was too matlaby
    is_in_bounds = np.full(np.shape(lat_grid1), False, dtype=bool)
    num_rows, num_cols = np.shape(lat_grid1)
    for row_idx in range(num_rows):
        for col_idx in range(num_cols):
            lon_val = lon_grid1[row_idx, col_idx]
            lat_val = lat_grid1[row_idx, col_idx]
            dist_buff_arr = distance_deg(lon_val, lat_val, \
                lon_grid1[is_inside_buffer], lat_grid1[is_inside_buffer])
            #if np.nanmin(dist_buff_arr) < swarm_reach:
            if np.min(dist_buff_arr) < swarm_reach:
                is_in_bounds[row_idx, col_idx] = True

    ring_score_dict = ring_score_calc(x_grid_offset_gcd, y_grid_offset_gcd, data_grid1, \
        mi_dict['sensor'], is_in_bounds, min_radius_deg, max_radius_deg)


    combo_score_dict = {}
    combo_score_dict['lon_grid1'] = lon_grid1
    combo_score_dict['lat_grid1'] = lat_grid1
    combo_score_dict['data_grid1'] = data_grid1
    combo_score_dict['spiral_score_grid'] = spiral_score_grid
    combo_score_dict['penalty_grid'] = penalty_grid
    combo_score_dict['ring_score_grid'] = ring_weight * ring_score_dict['ring_score_grid']
    combo_score_dict['ring_radius_grid'] = ring_score_dict['ring_radius_grid']
    combo_score_dict['ring_score_grid_full'] = ring_score_dict['ring_score_grid_full']
    combo_score_dict['radial_gradient_4d'] = ring_score_dict['radial_gradient_4d']
    combo_score_dict['fraction_input'] = fraction_input

    return combo_score_dict



def spiral_center_calc(x_grid_offset_gcd, y_grid_offset_gcd, data_grid1, sensor_type,
        op_lon, op_lat, filter_radius_deg, spiral_search_radius_deg, spiral_spacing_deg):

    alpha = 5 * np.pi/180

    # Sensor-specfic settings
    if sensor_type == '37GHz':
        outside_factor = 0.62
        data_grid1 = -data_grid1
    elif sensor_type == 'IR' or sensor_type == 'Vis':
        outside_factor = 0.50
    else:
        outside_factor = 0.62

    # Cut down to a usable disk, surrounded by nans
    in_filter_disk = x_grid_offset_gcd**2 + y_grid_offset_gcd**2 <= filter_radius_deg**2
    disk_img = np.nan * data_grid1
    disk_img[in_filter_disk] = data_grid1[in_filter_disk]

    # Make 1D arrays of just the clean points. "clean" means no nans
    is_clean = ~np.isnan(disk_img)
    disk_img_clean = disk_img[is_clean]
    x_grid_offset_gcd_clean = x_grid_offset_gcd[is_clean]
    y_grid_offset_gcd_clean = y_grid_offset_gcd[is_clean]
    lon_inc = x_grid_offset_gcd[0,1] - x_grid_offset_gcd[0,0]
    lat_inc = y_grid_offset_gcd[0,0] - y_grid_offset_gcd[1,0]

    grad_n, grad_e = np.gradient(disk_img, lon_inc, lat_inc)
    grad_n = -grad_n
    grad_n_clean = grad_n[is_clean]
    grad_e_clean = grad_e[is_clean]
    grad_orig_mag_clean = np.sqrt(grad_n_clean**2 + grad_e_clean**2)
    grad_log_mag_clean = np.log(1 + grad_orig_mag_clean)
    with np.errstate(divide='ignore', invalid='ignore'):
        grad_log_reduction_clean = grad_log_mag_clean / grad_orig_mag_clean
    grad_n_log_clean = grad_log_reduction_clean * grad_n_clean
    grad_e_log_clean = grad_log_reduction_clean * grad_e_clean


    # 1. Iterate the cross product score on a coarse grid
    off_arr = np.arange(-spiral_search_radius_deg, spiral_search_radius_deg+1e-6, spiral_spacing_deg)

    all_center_mean_cross = np.zeros((len(off_arr), len(off_arr))) * np.nan
    all_center_xs = np.zeros((len(off_arr), len(off_arr))) * np.nan
    all_center_ys = np.zeros((len(off_arr), len(off_arr))) * np.nan

    # Search out (search_radius_deg) degrees from the center point
    for row_idx, x_off in enumerate(off_arr):

        for col_idx, y_off in enumerate(off_arr):

            all_center_xs[row_idx, col_idx] = x_off
            all_center_ys[row_idx, col_idx] = y_off

            if x_off**2 + y_off**2 > (spiral_search_radius_deg + 2*spiral_spacing_deg/3)**2:
                pass # This takes out corners to save time

            else:

                proxy_x_clean = x_grid_offset_gcd_clean - x_off
                proxy_y_clean = y_grid_offset_gcd_clean - y_off

                spiral_x_clean = (alpha * proxy_x_clean + np.sign(op_lat) * proxy_y_clean) \
                    / np.sqrt((1+alpha**2) * (proxy_x_clean**2 + proxy_y_clean**2))
                spiral_y_clean = (alpha * proxy_y_clean - np.sign(op_lat) * proxy_x_clean) \
                    / np.sqrt((1+alpha**2) * (proxy_x_clean**2 + proxy_y_clean**2))

                raw_cross_score = spiral_x_clean * grad_n_log_clean - \
                    spiral_y_clean * grad_e_log_clean
                cross_score_clean = np.maximum(0, -raw_cross_score) + \
                    outside_factor * np.maximum(0, raw_cross_score)
                is_nan_cross = np.isnan(raw_cross_score)
                cross_score_clean[is_nan_cross] = np.nan
                #norm_mean_cross_score_clean = np.mean(cross_score_clean[~np.isnan(cross_score_clean)])
                norm_mean_cross_score_clean = np.nanmean(cross_score_clean)

                all_center_mean_cross[row_idx, col_idx] = norm_mean_cross_score_clean

    # 2. Search for the best full-resolution grid cell by (cubic?) interpolation
    sp_grid = intbx.InterpSectionToGlobalRectGrid(all_center_xs, all_center_ys, all_center_mean_cross, 
        x_grid_offset_gcd[0,:], y_grid_offset_gcd[:,0], 'linear')

    # Clean out the dodgy edge values
    sp_grid[x_grid_offset_gcd**2 + y_grid_offset_gcd**2 >= spiral_search_radius_deg**2] = np.nan

    fraction_input = len(disk_img_clean)/np.sum(in_filter_disk)

    return sp_grid, fraction_input


def ring_score_calc(x_grid_offset_gcd, y_grid_offset_gcd, data_grid1, \
        sensor_type, is_in_bounds, min_radius_deg, max_radius_deg):

    # Parameters
    inc = 5
    ang_deg_arr = np.arange(0, 360, inc)
    n_ang = len(ang_deg_arr)
    ring_point_thresh = 0.425 * n_ang

    # Add reversal step for 37GHz b/c eyes are *colder*:
    print('sensor_type = ', sensor_type)
    if '37' in sensor_type:
        data_grid1 = 450 - data_grid1

    # Calculate the modified gradient field
    # Somehow, the step of *1.14* is better than 1 for the *final* result. Couldn't figure out why.
    grad_n, grad_e = np.gradient((data_grid1**.333), -1.14, 1.14) # (This has to reverse the Matlab formula)

    # Translate degrees to pixels
    deg_per_pix = np.abs(y_grid_offset_gcd[0,0] - y_grid_offset_gcd[1,0])

    # Initialize variables related to the score grids
    n_rows, n_cols = np.shape(data_grid1)
    ring_score_grid = np.nan * np.zeros(np.shape(data_grid1))
    ring_radius_grid = np.zeros(np.shape(data_grid1))
    max_eye_bt_grid = np.nan * np.zeros(np.shape(data_grid1))

    # Unit vectors pointed radially inward
    ring_unit_vector_x = -np.cos(np.pi/180 * ang_deg_arr)
    ring_unit_vector_y = -np.sin(np.pi/180 * ang_deg_arr)

    # Create offset grids 
    mid_row = np.round(n_rows/2) # =32
    mid_col = np.round(n_cols/2)
    off_col_grid, off_row_grid = np.meshgrid(range(0, n_cols) - mid_col, range(0, n_rows) - mid_row) # -32...32

    # Convert offset grids to offset arrays
    is_in_radius_range = \
        x_grid_offset_gcd**2 + y_grid_offset_gcd**2 < (max_radius_deg + deg_per_pix)**2
    off_col_pts = off_col_grid[is_in_radius_range]
    off_row_pts = off_row_grid[is_in_radius_range]
    off_x_pts = x_grid_offset_gcd[is_in_radius_range]
    off_y_pts = y_grid_offset_gcd[is_in_radius_range]

    # Build the score grids. Iterate by radius, and within that, iterate by location
    print('Radius (deg) = ', end='')
    rad_arr = np.arange(min_radius_deg, max_radius_deg+1e-6, 0.05)

    # Extra part for ERC calcs
    n_rad = len(rad_arr)
    ring_score_grid_full = np.zeros((n_rows, n_cols, n_rad)) * np.nan
    radial_gradient_4d = np.zeros((n_rows, n_cols, n_rad, n_ang)) * np.nan

    for rad_idx in reversed(range(n_rad)):

        radius_deg = rad_arr[rad_idx]
        print('%4.2f' % (radius_deg), '', end='')

        # Calculate the row/col offsets for any center point at this radius
        row_off_arr = np.zeros(np.shape(ang_deg_arr))
        col_off_arr = np.zeros(np.shape(ang_deg_arr))
        for ang_idx, ang_deg in enumerate(ang_deg_arr):
            ring_x_pt = radius_deg * np.cos(np.pi/180 * ang_deg)
            ring_y_pt = radius_deg * np.sin(np.pi/180 * ang_deg)
            distance_arr = (off_x_pts - ring_x_pt)**2 + (off_y_pts - ring_y_pt)**2
            nearest_idx = np.argmin(distance_arr)
            row_off_arr[ang_idx] = off_row_pts[nearest_idx]
            col_off_arr[ang_idx] = off_col_pts[nearest_idx]

        radius_factor = radius_deg ** 0.1

        # Iterate by location
        for i, line in enumerate(is_in_bounds):
            for j, val in enumerate(line):
                
                if val == 0:
                    continue

                # Resolve the row/cols of the candidate ring at this point
                ring_row_arr = (i + row_off_arr).astype(int)
                ring_col_arr = (j + col_off_arr).astype(int)

                # Assign the gradients on the ring, or NaN if off the image
                is_in_image_box = (ring_row_arr >= 0) & (ring_row_arr < n_rows) & \
                    (ring_col_arr >= 0) & (ring_col_arr < n_cols)
                ring_gradient_x_arr = np.nan * np.zeros(np.shape(ring_col_arr))
                ring_gradient_y_arr = np.nan * np.zeros(np.shape(ring_row_arr))

                ring_gradient_x_arr[is_in_image_box] = \
                    grad_e[ring_row_arr[is_in_image_box], ring_col_arr[is_in_image_box]]
                ring_gradient_y_arr[is_in_image_box] = \
                    grad_n[ring_row_arr[is_in_image_box], ring_col_arr[is_in_image_box]]

                # Calculate the score at this location
                dot_product_arr = ring_unit_vector_x * ring_gradient_x_arr + \
                    ring_unit_vector_y * ring_gradient_y_arr
                dot_score_arr = dot_product_arr[~np.isnan(dot_product_arr)]

                # Put the score in the score grid if it comes from enough valid points
                n_real = len(dot_score_arr)
                if n_real <= ring_point_thresh:
                    # If it doesn't qualify for the smallest, then none of
                    # the other diameters get to count either (will be NaN in the
                    # max or inequality function)
                    # ring_score = np.nan

                    # Try to be nicer about it instead:
                    ring_score = 0

                else:
                    ring_score = radius_factor * np.mean(dot_score_arr)

                # If the high score was beaten, then set this value in the score grid
                if rad_idx == len(rad_arr)-1 or ring_score > ring_score_grid[i, j]:
                    ring_score_grid[i, j] = ring_score
                    ring_radius_grid[i, j] = radius_deg

                #print('i, j, ring score = ', i, j, ring_score)

                # Save info for "score by radius" array
                ring_score_grid_full[i, j, rad_idx] = ring_score

                # Save all info in massive radial_gradient_4d tensor
                radial_gradient_4d[i, j, rad_idx, :] = dot_product_arr

    print('')

    # Assinging the warmest pixel corresponding to each ring
    for i, line in enumerate(is_in_bounds):
        for j, val in enumerate(line):
            
            if val == 0:
                continue

            # Ring of points
            radius_deg = ring_radius_grid[i, j]

            # Assign warmest pixel value to output grid
            is_in_eye_grid = (x_grid_offset_gcd - x_grid_offset_gcd[i,j])**2 + \
                (y_grid_offset_gcd - y_grid_offset_gcd[i,j])**2 <= radius_deg**2

            if np.sum(is_in_eye_grid) > 0:
                max_eye_bt_grid[i, j] = np.nanmax(data_grid1[is_in_eye_grid])
            else:
                max_eye_bt_grid[i, j] = np.nan
   

    # Put output into dictionary
    ring_score_dict = {}
    ring_score_dict['ring_score_grid'] = ring_score_grid
    ring_score_dict['ring_radius_grid'] = ring_radius_grid
    ring_score_dict['max_eye_bt_grid'] = max_eye_bt_grid
    ring_score_dict['ring_score_grid_full'] = ring_score_grid_full
    ring_score_dict['radial_gradient_4d'] = radial_gradient_4d

    return ring_score_dict


def quality_check(score_dict):

    n_rows, n_cols = np.shape(score_dict['spiral_score_grid'])
    spiral_score_grid = score_dict['spiral_score_grid'] + 0
    spiral_score_grid[np.isnan(spiral_score_grid)] = -1e9
    max_idx = np.argmax(spiral_score_grid)
    i_max_score, j_max_score = ind2sub((n_rows, n_cols), max_idx)
    #i_max_score, j_max_score = np.unravel_index(max_idx, (n_rows, n_cols)) # This works too

    if score_dict['fraction_input'] < 0.5:
        uses_target = False

    elif np.isnan(spiral_score_grid[i_max_score, j_max_score]):
        uses_target = False

    elif i_max_score <= 1 or j_max_score <= 1:
        uses_target = False

    elif i_max_score >= n_rows-2 or j_max_score >= n_cols-2:
        uses_target = False

    elif np.isnan(spiral_score_grid[i_max_score-2, j_max_score]):
        uses_target = False

    elif np.isnan(spiral_score_grid[i_max_score+2, j_max_score]):
        uses_target = False

    elif np.isnan(spiral_score_grid[i_max_score, j_max_score-2]):
        uses_target = False

    elif np.isnan(spiral_score_grid[i_max_score, j_max_score+2]):
        uses_target = False

    else:
        uses_target = True

    return uses_target



