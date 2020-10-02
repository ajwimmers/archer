from __future__ import print_function
import time, calendar
import numpy as np
import netCDF4
from os import remove
from os.path import isfile
import utilities.ScoreFuncs as sftbx


def write_archer_out_nc(filename, mi_dict, score_dict, stats_dict, processed_dict):

    # This is necessary to fix a bug in the Netcdf4 module with file permission settings
    if isfile(filename):
        remove(filename)

    rootgrp = netCDF4.Dataset(filename, 'w' ,format='NETCDF4')
    numRows, numCols, numRad = np.shape(score_dict['ring_score_grid_full'])
    rootgrp.createDimension('lat', numRows)
    rootgrp.createDimension('lon', numCols)
    rootgrp.createDimension('time', None)
    rootgrp.createDimension('radius', numRad)

    nc_tc_time = rootgrp.createVariable('tc_time', 'i8', ('time',))
    #nc_tc_time[:] = stats_dict['tc_time']
    nc_tc_time = stats_dict['tc_time']

    nc_sensor_type = rootgrp.createVariable('sensor_type', 'c')
    nc_sensor_type = stats_dict['sensor_type']
    #nc_sensor_type[:] = stats_dict['sensor_type']

    nc_fx_lon = rootgrp.createVariable('fx_lon', 'f8')
    nc_fx_lon[:] = mi_dict['op_lon']
    nc_fx_lat = rootgrp.createVariable('op_lat', 'f8')
    nc_fx_lat[:] = mi_dict['op_lat']
    nc_est_vmax = rootgrp.createVariable('op_vmax', 'f8')
    nc_est_vmax[:] = mi_dict['op_vmax']

    nc_fraction_input = rootgrp.createVariable('fraction_input', 'f8')
    nc_fraction_input[:] = score_dict['fraction_input']

    nc_target_lon = rootgrp.createVariable('target_lon', 'f8')
    nc_target_lon[:] = stats_dict['target_lon']
    nc_target_lat = rootgrp.createVariable('target_lat', 'f8')
    nc_target_lat[:] = stats_dict['target_lat'] 
    nc_weak_target_lon = rootgrp.createVariable('weak_target_lon', 'f8')
    nc_weak_target_lon[:] = stats_dict['weak_target_lon']
    nc_weak_target_lat = rootgrp.createVariable('weak_target_lat', 'f8')
    nc_weak_target_lat[:] = stats_dict['weak_target_lat']
    nc_confidence_score = rootgrp.createVariable('confidence_score', 'f8')
    nc_confidence_score[:] = stats_dict['confidence_score']
    nc_alpha_parameter = rootgrp.createVariable('alpha_parameter', 'f8')
    nc_alpha_parameter[:] = stats_dict['alpha_parameter']
    nc_radius50percCertDeg = rootgrp.createVariable('radius50percCertDeg', 'f8')
    nc_radius50percCertDeg[:] = stats_dict['radius50percCertDeg']
    nc_radius95percCertDeg = rootgrp.createVariable('radius95percCertDeg', 'f8')
    nc_radius95percCertDeg[:] = stats_dict['radius95percCertDeg']

    nc_lon_arr = rootgrp.createVariable('lon_arr', 'f8', ('lon',))
    nc_lon_arr[:] = score_dict['lon_grid1'][0,:]
    nc_lat_arr = rootgrp.createVariable('lat_arr', 'f8', ('lat',))
    nc_lat_arr[:] = score_dict['lat_grid1'][:,0] # OR REVERSE THIS MAYBE?

    nc_data_grid = rootgrp.createVariable('data_grid', 'f8', ('lat','lon',))
    nc_data_grid[:] = score_dict['data_grid1']
    nc_spiral_score_grid = rootgrp.createVariable('spiral_score_grid', 'f8', ('lat','lon',), least_significant_digit=3)
    nc_spiral_score_grid[:] = score_dict['spiral_score_grid']
    nc_penalty_grid = rootgrp.createVariable('penalty_grid', 'f8', ('lat','lon',), least_significant_digit=3)
    nc_penalty_grid[:] = score_dict['penalty_grid']
    nc_ring_score_grid = rootgrp.createVariable('ring_score_grid', 'f8', ('lat','lon',), least_significant_digit=3)
    nc_ring_score_grid[:] = score_dict['ring_score_grid']

    # Get combo_max vars, locs
    max_idx = np.argmax(processed_dict['combo_score_grid'])
    i_combo_max, j_combo_max = sftbx.ind2sub(np.shape(processed_dict['combo_score_grid']), max_idx)

    nc_combo_score = rootgrp.createVariable('combo_score', 'f8')
    nc_combo_score[:] = np.nanmax(processed_dict['combo_score_grid'])
    nc_ring_score = rootgrp.createVariable('ring_score', 'f8')
    nc_ring_score[:] = np.nanmax(score_dict['ring_score_grid'])
    nc_ring_radius_deg = rootgrp.createVariable('nc_ring_radius_deg', 'f8')
    nc_ring_radius_deg[:] = score_dict['ring_radius_grid'][i_combo_max, j_combo_max]
    nc_score_by_radius_arr = rootgrp.createVariable('score_by_radius_arr', 'f8', ('radius',))
    nc_score_by_radius_arr[:] = np.squeeze(score_dict['ring_score_grid_full'][i_combo_max, j_combo_max, :])
    #nc_gradient_grid = rootgrp.createVariable('gradient_grid', 'f8', ('lat','lon',), least_significant_digit=3)
    #nc_gradient_grid[:] = # I forget what this is, and whether the dimensions (above line) are right.
    nc_ring_radius_grid = rootgrp.createVariable('ring_radius_grid', 'f8', ('lat','lon',), least_significant_digit=3)
    nc_ring_radius_grid[:] = score_dict['ring_radius_grid']

    rootgrp.close()
    return 'completed'


    """
    # Pack up the output values

    processed_dict = {}
    processed_dict['combo_score_grid'] = combo_grid

    out_dir = './'
    out_filename = out_dir + os.path.splitext(os.path.basename(filename))[0] + '_archer.nc'

    write_status = nctbx.write_archer_out_nc(out_filename, mi_dict, score_dict, stats_dict, processed_dict)

    output = {}
    output['tc_time'] = tc_time
    output['lon_arr'] = score_dict['lon_grid1'][0,:]
    output['lat_arr'] = score_dict['lat_grid1'][:,0]
    output['data_grid'] = score_dict['data_grid1']
    output['spiral_score_grid'] = score_dict['spiral_score_grid']
    output['penalty_grid'] = score_dict['penalty_grid']
    output['ring_score_grid'] = score_dict['ring_score_grid']
    output['combo_score_grid'] = combo_grid
    output['combo_score'] = combo_score
    output['ring_score'] = ring_score
    output['score_by_radius_arr'] = score_by_radius_arr
    output['gradient_grid'] = gradient_grid
    output['fraction_input'] = score_dict['fraction_input']
    output['sensor_type'] = sensor_type
    output['sat'] = sat ???

    fx_lon, fx_lat, est_vmax ::: mi_dict

    target_lon/lat ::: lon_combo_max/lat_combo_max
    weaktargetLon/Lat ::: skip for now

    ring_radius_grid :::  combo_score_dict['ring_radius_grid']
    From earlier: ring_radius_deg = score_dict['ring_radius_deg'][i_combo_max, j_combo_max]

    confidence_score: confidence_score
    poissonAlpha (should be alpha parameter): alpha

    radius50percCertDeg: rad50
    radius95percCertDeg: rad95

    

    """
