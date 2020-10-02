
from scipy.io import loadmat
import time, calendar, datetime
import numpy as np
import netCDF4, h5py # Note that you might get funny behavior if you import these in reverse order
from archer.archer4  import archer4

# This instructs the script to calculate all the examples. Otherwise, edit the 
# script on a case by case basis:
calc_all = False

# Refer to the header of archer4.py for an explanation of all the algorithm inputs/outputs.



def convert_matlab_time(matlab_datenum):

    days_0000_to_1970 = 719528
    secs_in_day = 86400

    # Seconds since 1970 for time module object:
    time_secs = [calendar.timegm(time.gmtime((i - days_0000_to_1970) * secs_in_day)) for i in matlab_datenum]

    # Convert list to numpy array
    if len(time_secs) > 1:
        time_secs = np.asarray(time_secs)
    else:
        time_secs = time_secs[0]

    return time_secs


if False or calc_all:

    filename = './data/AMSU-B_N16_L1B_BRIGHTNESS_20050921_2038.nc'
    with netCDF4.Dataset(filename, 'r') as nc:
        bt_mx  = nc.variables['amsub_ch1'][:]
        lat_mx = nc.variables['latitude'][:]
        lon_mx = nc.variables['longitude'][:]

    # Set approximate first guesses
    op_lon = -87.0
    op_lat = 24.4
    op_vmax = 125
    image_time = calendar.timegm(datetime.datetime(2005, 9, 21, 20, 38, 0).timetuple())


    # Collect for ARCHER

    image = {}
    image['data_grid'] = bt_mx
    image['lat_grid'] = lat_mx
    image['lon_grid'] = lon_mx

    attrib = {}
    attrib['sat'] = 'NOAA-16'
    attrib['sensor'] = 'AMSU-B'
    attrib['scan_type'] = 'Crosstrack'
    attrib['archer_channel_type'] = '89GHz'

    first_guess = {}
    first_guess['source'] = 'fx'
    first_guess['time'] = image_time
    first_guess['vmax'] = op_vmax
    first_guess['lat'] = op_lat
    first_guess['lon'] = op_lon


    in_dict, out_dict, score_dict = \
        archer4(image, attrib, first_guess, para_fix=True, display_filename='imgs/demo_amsub_89.png')


if False or calc_all:

    filename = './data/MHS_N18_L1B_BRIGHTNESS_20050921_1910.nc'
    with netCDF4.Dataset(filename, 'r') as nc:
        bt_mx  = nc.variables['mhs_ch1'][:]
        lat_mx = nc.variables['latitude'][:]
        lon_mx = nc.variables['longitude'][:]

    # Set approximate first guesses
    op_lon = -87.0
    op_lat = 24.4
    op_vmax = 125
    image_time = calendar.timegm(datetime.datetime(2005, 9, 21, 19, 10, 0).timetuple())


    # Collect for ARCHER

    image = {}
    image['data_grid'] = bt_mx
    image['lat_grid'] = lat_mx
    image['lon_grid'] = lon_mx

    attrib = {}
    attrib['sat'] = 'NOAA-18'
    attrib['sensor'] = 'MHS'
    attrib['scan_type'] = 'Crosstrack'
    attrib['archer_channel_type'] = '89GHz'

    first_guess = {}
    first_guess['source'] = 'fx'
    first_guess['time'] = image_time
    first_guess['vmax'] = op_vmax
    first_guess['lat'] = op_lat
    first_guess['lon'] = op_lon


    in_dict, out_dict, score_dict = \
        archer4(image, attrib, first_guess, para_fix=True, display_filename='imgs/demo_mhs_89.png')


if False or calc_all:

    data = loadmat('./data/20200208T094041_ssmis18_85.mat')

    image = {}
    image['data_grid'] = data['mi'][0][0]['bthMx']
    image['lat_grid'] = data['mi'][0][0]['rawLatMx']
    image['lon_grid'] = data['mi'][0][0]['rawLonMx']

    attrib = {}
    attrib['sat'] = 'F18'
    attrib['sensor'] = 'SSMIS'
    attrib['scan_type'] = 'Conical'
    attrib['archer_channel_type'] = '89GHz'

    first_guess = {}
    first_guess['source'] = 'fx'
    first_guess['time'] = convert_matlab_time(data['mi'][0][0]['tcTime']) 
    first_guess['vmax'] = data['mi'][0][0]['estVmax']
    first_guess['lat'] = data['mi'][0][0]['estCenterLat']
    first_guess['lon'] = data['mi'][0][0]['estCenterLon']


    in_dict, out_dict, score_dict = \
        archer4(image, attrib, first_guess, para_fix=True, display_filename='imgs/demo_conic_89h.png')


if False or calc_all:

    data = loadmat('./data/20200208T094041_ssmis18_85.mat')

    image = {}
    image['data_grid'] = data['mi'][0][0]['bthMx']
    image['lat_grid'] = data['mi'][0][0]['rawLatMx']
    image['lon_grid'] = data['mi'][0][0]['rawLonMx']

    # Convert to partial scan just to test something weirder
    image['data_grid'][:, 73:] = np.nan


    attrib = {}
    attrib['sat'] = 'F18'
    attrib['sensor'] = 'SSMIS'
    attrib['scan_type'] = 'Conical'
    attrib['archer_channel_type'] = '89GHz'

    first_guess = {}
    first_guess['source'] = 'fx'
    first_guess['time'] = convert_matlab_time(data['mi'][0][0]['tcTime']) 
    first_guess['vmax'] = data['mi'][0][0]['estVmax']
    first_guess['lat'] = data['mi'][0][0]['estCenterLat']
    first_guess['lon'] = data['mi'][0][0]['estCenterLon']


    in_dict, out_dict, score_dict = \
        archer4(image, attrib, first_guess, para_fix=True, display_filename='imgs/demo_conic_partial_89h.png')


if False or calc_all:

    data = loadmat('./data/20200208T061548_amsr2_37.mat')

    image = {}
    image['data_grid'] = data['mi'][0][0]['bthMx']
    image['lat_grid'] = data['mi'][0][0]['rawLatMx']
    image['lon_grid'] = data['mi'][0][0]['rawLonMx']

    attrib = {}
    attrib['sat'] = 'Aqua'
    attrib['sensor'] = 'AMSR2'
    attrib['scan_type'] = 'Conical'
    attrib['archer_channel_type'] = '37GHz'

    first_guess = {}
    first_guess['source'] = 'fx'
    first_guess['time'] = convert_matlab_time(data['mi'][0][0]['tcTime']) 
    first_guess['vmax'] = data['mi'][0][0]['estVmax']
    first_guess['lat'] = data['mi'][0][0]['estCenterLat']
    first_guess['lon'] = data['mi'][0][0]['estCenterLon']


    in_dict, out_dict, score_dict = \
        archer4(image, attrib, first_guess, para_fix=True, display_filename='imgs/demo_conic_37h.png')


if False or calc_all:

    data = loadmat('./data/20200208T083000_ir.mat')

    image = {}
    image['data_grid'] = data['aSetSaved'][0][0]['tb']
    image['lat_grid'] = data['aSetSaved'][0][0]['latMatrix']
    image['lon_grid'] = data['aSetSaved'][0][0]['lonMatrix']

    attrib = {}
    attrib['sat'] = 'Himawari-8'
    attrib['sensor'] = 'Imager'
    attrib['scan_type'] = 'Geo'
    attrib['archer_channel_type'] = 'IR'
    attrib['nadir_lon'] = 140.7

    first_guess = {}
    first_guess['source'] = 'fx'
    first_guess['time'] = convert_matlab_time(data['aSetSaved'][0][0]['matDate']) 
    first_guess['vmax'] = data['aSetSaved'][0][0]['estVmax'][0][0]
    first_guess['lat'] = data['aSetSaved'][0][0]['estCenterLat'][0][0]
    first_guess['lon'] = data['aSetSaved'][0][0]['estCenterLon'][0][0]


    in_dict, out_dict, score_dict = \
        archer4(image, attrib, first_guess, para_fix=True, display_filename='imgs/demo_geo_ir.png')


if False or calc_all:

    data = loadmat('./data/20200208T083000_sw.mat')

    image = {}
    image['data_grid'] = data['aSetSaved'][0][0]['tb']
    image['lat_grid'] = data['aSetSaved'][0][0]['latMatrix']
    image['lon_grid'] = data['aSetSaved'][0][0]['lonMatrix']

    attrib = {}
    attrib['sat'] = 'Himawari-8'
    attrib['sensor'] = 'Imager'
    attrib['scan_type'] = 'Geo'
    attrib['archer_channel_type'] = 'SWIR'
    attrib['nadir_lon'] = 140.7 #41.5

    first_guess = {}
    first_guess['source'] = 'fx'
    first_guess['time'] = convert_matlab_time(data['aSetSaved'][0][0]['matDate']) 
    first_guess['vmax'] = data['aSetSaved'][0][0]['estVmax'][0][0]
    first_guess['lat'] = data['aSetSaved'][0][0]['estCenterLat'][0][0]
    first_guess['lon'] = data['aSetSaved'][0][0]['estCenterLon'][0][0]


    in_dict, out_dict, score_dict = \
        archer4(image, attrib, first_guess, para_fix=True, display_filename='imgs/demo_geo_sw.png')


if False or calc_all:

    data = loadmat('./data/20200208T083000_vis.mat')

    image = {}
    image['data_grid'] = data['aSetSaved'][0][0]['bv']
    image['lat_grid'] = data['aSetSaved'][0][0]['latMatrix']
    image['lon_grid'] = data['aSetSaved'][0][0]['lonMatrix']

    attrib = {}
    attrib['sat'] = 'Himawari-8'
    attrib['sensor'] = 'Imager'
    attrib['scan_type'] = 'Geo'
    attrib['archer_channel_type'] = 'Vis'
    attrib['nadir_lon'] = 140.7 #41.5

    first_guess = {}
    first_guess['source'] = 'fx'
    first_guess['time'] = convert_matlab_time(data['aSetSaved'][0][0]['matDate']) 
    first_guess['vmax'] = data['aSetSaved'][0][0]['estVmax'][0][0]
    first_guess['lat'] = data['aSetSaved'][0][0]['estCenterLat'][0][0]
    first_guess['lon'] = data['aSetSaved'][0][0]['estCenterLon'][0][0]


    in_dict, out_dict, score_dict = \
        archer4(image, attrib, first_guess, para_fix=True, display_filename='imgs/demo_geo_vis.png')


if False or calc_all:

    # Apply data grabbed from a CLASS file and converted to netcdf with satpy:

    filename = './data/viirs_sdr_i01_damien_test.nc'
    with netCDF4.Dataset(filename, 'r') as nc:
        bv_mx  = nc.variables['I01'][:]
        lat_mx = nc.variables['latitude'][:]
        lon_mx = nc.variables['longitude'][:]

    # Set approximate first guesses
    op_lon = 116.8
    op_lat = -20.8
    op_vmax = 90
    image_time = calendar.timegm(datetime.datetime(2020, 2, 8, 5, 43, 0).timetuple())


    # Collect for ARCHER

    # Adjust brightness value to match geo: ####################
    bv_mx = bv_mx * 255 / 100 * 0.8 # Assume it varies from 0-100, and use 0.8 as an adjustment.

    image = {}
    image['data_grid'] = bv_mx
    image['lat_grid'] = lat_mx
    image['lon_grid'] = lon_mx

    attrib = {}
    attrib['sat'] = 'JPSS-1'
    attrib['sensor'] = 'VIIRS'
    attrib['scan_type'] = 'Crosstrack'
    attrib['archer_channel_type'] = 'Vis'

    first_guess = {}
    first_guess['source'] = 'fx'
    first_guess['time'] = image_time
    first_guess['vmax'] = op_vmax
    first_guess['lat'] = op_lat
    first_guess['lon'] = op_lon


    in_dict, out_dict, score_dict = \
        archer4(image, attrib, first_guess, para_fix=True, display_filename='imgs/demo_viisr_vis.png')


if False or calc_all:

    # Apply data grabbed from a CLASS file and converted to netcdf with satpy:
    # HOWEVER, IT LOOKS TO ME LIKE THIS SWIR IMAGE IS INCORRECTLY CALIBRATED FROM SATPY.
    # THE RANGE SEEMS TOO NARROW.

    filename = './data/viirs_sdr_i04_damien_test.nc'
    with netCDF4.Dataset(filename, 'r') as nc:
        bt_mx  = nc.variables['I04'][:]
        lat_mx = nc.variables['latitude'][:]
        lon_mx = nc.variables['longitude'][:]

    # Set approximate first guesses
    op_lon = 116.8
    op_lat = -20.8
    op_vmax = 90
    image_time = calendar.timegm(datetime.datetime(2020, 2, 8, 5, 43, 0).timetuple())


    # Collect for ARCHER

    image = {}
    image['data_grid'] = bt_mx
    image['lat_grid'] = lat_mx
    image['lon_grid'] = lon_mx

    attrib = {}
    attrib['sat'] = 'NOAA-20'
    attrib['sensor'] = 'VIIRS'
    attrib['scan_type'] = 'Crosstrack'
    attrib['archer_channel_type'] = 'SWIR'

    first_guess = {}
    first_guess['source'] = 'fx'
    first_guess['time'] = image_time
    first_guess['vmax'] = op_vmax
    first_guess['lat'] = op_lat
    first_guess['lon'] = op_lon


    in_dict, out_dict, score_dict = \
        archer4(image, attrib, first_guess, para_fix=True, display_filename='imgs/demo_viisr_sw.png')


if False or calc_all:

    # Apply data grabbed from a CLASS file and converted to netcdf with satpy:

    filename = './data/viirs_sdr_i05_damien_test.nc'
    with netCDF4.Dataset(filename, 'r') as nc:
        bt_mx  = nc.variables['I05'][:]
        lat_mx = nc.variables['latitude'][:]
        lon_mx = nc.variables['longitude'][:]

    # Set approximate first guesses
    op_lon = 116.8
    op_lat = -20.8
    op_vmax = 90
    image_time = calendar.timegm(datetime.datetime(2020, 2, 8, 5, 43, 0).timetuple())


    # Collect for ARCHER

    image = {}
    image['data_grid'] = bt_mx
    image['lat_grid'] = lat_mx
    image['lon_grid'] = lon_mx

    attrib = {}
    attrib['sat'] = 'NOAA-20'
    attrib['sensor'] = 'VIIRS'
    attrib['scan_type'] = 'Crosstrack'
    attrib['archer_channel_type'] = 'IR'

    first_guess = {}
    first_guess['source'] = 'fx'
    first_guess['time'] = image_time
    first_guess['vmax'] = op_vmax
    first_guess['lat'] = op_lat
    first_guess['lon'] = op_lon


    in_dict, out_dict, score_dict = \
        archer4(image, attrib, first_guess, para_fix=True, display_filename='imgs/demo_viisr_ir.png')


if False or calc_all:

    # Apply data grabbed from a CLASS file and converted to netcdf with satpy:

    filename = './data/viirs_sdr_DNB_damien_test.nc'
    with netCDF4.Dataset(filename, 'r') as nc:
        bv_mx  = nc.variables['DNB'][:]
        lat_mx = nc.variables['latitude'][:]
        lon_mx = nc.variables['longitude'][:]

    # Set approximate first guesses
    op_lon = 116.8
    op_lat = -20.8
    op_vmax = 90
    image_time = calendar.timegm(datetime.datetime(2020, 2, 8, 5, 43, 0).timetuple())


    # Adjust brightness value to match geo: ####################
    bv_mx = bv_mx / 0.018 * 255 # This is ad hoc. Assume that it has been normalized to 0-255 like geo

    # Collect for ARCHER

    image = {}
    image['data_grid'] = bv_mx
    image['lat_grid'] = lat_mx
    image['lon_grid'] = lon_mx

    attrib = {}
    attrib['sat'] = 'NOAA-20'
    attrib['sensor'] = 'VIIRS'
    attrib['scan_type'] = 'Crosstrack'
    attrib['archer_channel_type'] = 'DNB'

    first_guess = {}
    first_guess['source'] = 'fx'
    first_guess['time'] = image_time
    first_guess['vmax'] = op_vmax
    first_guess['lat'] = op_lat
    first_guess['lon'] = op_lon


    in_dict, out_dict, score_dict = \
        archer4(image, attrib, first_guess, para_fix=True, display_filename='imgs/demo_viisr_dnb.png')


if False or calc_all:

    # Apply 89GHz data grabbed from a CLASS file:

    filename = './data/GATMO-SATMS_j01_d20200208_t0536586_e0544583_b11518_c20200415211227285695_noac_ops.h5'

    with h5py.File(filename, 'r') as f:
        bt_mx =  f['All_Data']['ATMS-SDR_All']['BrightnessTemperature'][()][:,:,15]
        lat_mx = f['All_Data']['ATMS-SDR-GEO_All']['Latitude'][()]
        lon_mx = f['All_Data']['ATMS-SDR-GEO_All']['Longitude'][()]
        zen_mx = f['All_Data']['ATMS-SDR-GEO_All']['SatelliteZenithAngle'][()]
        azm_mx = f['All_Data']['ATMS-SDR-GEO_All']['SatelliteAzimuthAngle'][()]

    bt_mx  =  bt_mx / 100. - 270 # Just a rough guess to get this going

    # Vmax
    op_lon = 116.8
    op_lat = -20.8
    op_vmax = 90
    image_time = calendar.timegm(datetime.datetime(2020, 2, 8, 5, 43, 0).timetuple())


    # Collect for ARCHER

    image = {}
    image['data_grid'] = bt_mx
    image['lat_grid'] = lat_mx
    image['lon_grid'] = lon_mx
    image['zen_grid'] = zen_mx
    image['azm_grid'] = azm_mx

    attrib = {}
    attrib['sat'] = 'NOAA-20'
    attrib['sensor'] = 'ATMS'
    attrib['scan_type'] = 'Crosstrack'
    attrib['archer_channel_type'] = '89GHz'

    first_guess = {}
    first_guess['source'] = 'fx'
    first_guess['time'] = image_time
    first_guess['vmax'] = op_vmax
    first_guess['lat'] = op_lat
    first_guess['lon'] = op_lon


    in_dict, out_dict, score_dict = \
        archer4(image, attrib, first_guess, para_fix=True, display_filename='imgs/demo_atms_89.png')


if True or calc_all:

    # Apply 183 GHz data grabbed from a CLASS file:

    filename = './data/GATMO-SATMS_j01_d20200208_t0536586_e0544583_b11518_c20200415211227285695_noac_ops.h5'

    with h5py.File(filename, 'r') as f:
        bt_mx =  f['All_Data']['ATMS-SDR_All']['BrightnessTemperature'][()][:,:,17]
        lat_mx = f['All_Data']['ATMS-SDR-GEO_All']['Latitude'][()]
        lon_mx = f['All_Data']['ATMS-SDR-GEO_All']['Longitude'][()]
        zen_mx = f['All_Data']['ATMS-SDR-GEO_All']['SatelliteZenithAngle'][()]
        azm_mx = f['All_Data']['ATMS-SDR-GEO_All']['SatelliteAzimuthAngle'][()]

    bt_mx  =  bt_mx / 100. - 270 # Just a rough guess to get this going

    # Vmax
    op_lon = 116.8
    op_lat = -20.8
    op_vmax = 90
    image_time = calendar.timegm(datetime.datetime(2020, 2, 8, 5, 43, 0).timetuple())


    # Collect for ARCHER

    image = {}
    image['data_grid'] = bt_mx
    image['lat_grid'] = lat_mx
    image['lon_grid'] = lon_mx
    image['zen_grid'] = zen_mx
    image['azm_grid'] = azm_mx

    attrib = {}
    attrib['sat'] = 'NOAA-20'
    attrib['sensor'] = 'ATMS'
    attrib['scan_type'] = 'Crosstrack'
    attrib['archer_channel_type'] = '183GHz'

    first_guess = {}
    first_guess['source'] = 'fx'
    first_guess['time'] = image_time
    first_guess['vmax'] = op_vmax
    first_guess['lat'] = op_lat
    first_guess['lon'] = op_lon


    in_dict, out_dict, score_dict = \
        archer4(image, attrib, first_guess, para_fix=True, display_filename='imgs/demo_atms_183.png')


if False or calc_all:

    data = loadmat('./data/20150826T163957_ssmis16.mat')

    image = {}
    image['data_grid'] = data['mi'][0][0]['bthMx']
    image['lat_grid'] = data['mi'][0][0]['rawLatMx']
    image['lon_grid'] = data['mi'][0][0]['rawLonMx']
    # State the lon grid properly as -180 to 180:
    image['lon_grid'] = np.mod(image['lon_grid'] + 180, 360) - 180

    attrib = {}
    attrib['sat'] = 'F16'
    attrib['sensor'] = 'SSMIS'
    attrib['scan_type'] = 'Conical'
    attrib['archer_channel_type'] = '89GHz'

    first_guess = {}
    first_guess['source'] = 'fx'
    first_guess['time'] = convert_matlab_time(data['mi'][0][0]['tcTime']) 
    first_guess['vmax'] = data['mi'][0][0]['estVmax'][0][0]
    first_guess['lat'] = data['mi'][0][0]['estCenterLat'][0][0]
    first_guess['lon'] = data['mi'][0][0]['estCenterLon'][0][0]


    in_dict, out_dict, score_dict = \
        archer4(image, attrib, first_guess, para_fix=True, display_filename='imgs/demo_conic_89h_antem.png')

    print(out_dict)