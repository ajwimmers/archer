
from scipy.io import loadmat
from archer.archer3_mw  import archer3_mw
from archer.archer3_visir import archer3_visir
import time, calendar
import numpy as np
import h5py


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


if False:

    data = loadmat('./data/20200208T094041_ssmis18_85.mat')

    image = {}
    image['bt_grid'] = data['mi'][0][0]['bthMx']
    image['lat_grid'] = data['mi'][0][0]['rawLatMx']
    image['lon_grid'] = data['mi'][0][0]['rawLonMx']
    image['time_arr'] = convert_matlab_time(data['mi'][0][0]['timeArray']) 

    attrib = {}
    attrib['sat'] = 'F18'
    attrib['sensor'] = 'SSMIS'
    attrib['scan_type'] = 'Conical'
    attrib['archer_channel_type'] = '89(H)GHz'

    first_guess = {}
    first_guess['source'] = 'fx' # 'fx', 'bt', 'archer'
    first_guess['time'] = convert_matlab_time(data['mi'][0][0]['tcTime']) 
    first_guess['vmax'] = data['mi'][0][0]['estVmax']
    first_guess['lat'] = data['mi'][0][0]['estCenterLat']
    first_guess['lon'] = data['mi'][0][0]['estCenterLon']


    in_dict, out_dict, score_dict = \
        archer3_mw(image, attrib, first_guess, display_filename='imgs/test1_89h.png')

if False:

    data = loadmat('./data/20200208T061548_amsr2_37.mat')

    image = {}
    image['bt_grid'] = data['mi'][0][0]['bthMx']
    image['lat_grid'] = data['mi'][0][0]['rawLatMx']
    image['lon_grid'] = data['mi'][0][0]['rawLonMx']
    image['time_arr'] = convert_matlab_time(data['mi'][0][0]['timeArray']) 

    attrib = {}
    attrib['sat'] = 'Aqua'
    attrib['sensor'] = 'AMSR2'
    attrib['scan_type'] = 'Conical'
    attrib['archer_channel_type'] = '37(H)GHz'

    first_guess = {}
    first_guess['source'] = 'fx' # 'fx', 'bt', 'archer'
    first_guess['time'] = convert_matlab_time(data['mi'][0][0]['tcTime']) 
    first_guess['vmax'] = data['mi'][0][0]['estVmax']
    first_guess['lat'] = data['mi'][0][0]['estCenterLat']
    first_guess['lon'] = data['mi'][0][0]['estCenterLon']


    in_dict, out_dict, score_dict = \
        archer3_mw(image, attrib, first_guess, display_filename='imgs/test1_37h.png')


if False:

    data = loadmat('./data/20200208T083000_ir.mat')

    image = {}
    image['bt_grid'] = data['aSetSaved'][0][0]['tb']
    image['lat_grid'] = data['aSetSaved'][0][0]['latMatrix']
    image['lon_grid'] = data['aSetSaved'][0][0]['lonMatrix']
    #image['time_arr'] = Does not exist

    attrib = {}
    attrib['sat'] = 'Himawari-8'
    attrib['sensor'] = 'Imager'
    attrib['scan_type'] = 'Geo'
    attrib['archer_channel_type'] = 'IR'
    attrib['nadir_lon'] = 140.7 #41.5

    first_guess = {}
    first_guess['source'] = 'fx' # 'fx', 'bt', 'archer'
    first_guess['time'] = convert_matlab_time(data['aSetSaved'][0][0]['matDate']) 
    first_guess['vmax'] = data['aSetSaved'][0][0]['estVmax']
    first_guess['lat'] = data['aSetSaved'][0][0]['estCenterLat']
    first_guess['lon'] = data['aSetSaved'][0][0]['estCenterLon']


    in_dict, out_dict, score_dict = \
        archer3_visir(image, attrib, first_guess, display_filename='imgs/test1_ir.png')


if False:

    data = loadmat('./data/20200208T083000_sw.mat')

    image = {}
    image['bt_grid'] = data['aSetSaved'][0][0]['tb']
    image['lat_grid'] = data['aSetSaved'][0][0]['latMatrix']
    image['lon_grid'] = data['aSetSaved'][0][0]['lonMatrix']
    #image['time_arr'] = Does not exist

    attrib = {}
    attrib['sat'] = 'Himawari-8'
    attrib['sensor'] = 'Imager'
    attrib['scan_type'] = 'Geo'
    attrib['archer_channel_type'] = 'SWIR'
    attrib['nadir_lon'] = 140.7 #41.5

    first_guess = {}
    first_guess['source'] = 'fx' # 'fx', 'bt', 'archer'
    first_guess['time'] = convert_matlab_time(data['aSetSaved'][0][0]['matDate']) 
    first_guess['vmax'] = data['aSetSaved'][0][0]['estVmax']
    first_guess['lat'] = data['aSetSaved'][0][0]['estCenterLat']
    first_guess['lon'] = data['aSetSaved'][0][0]['estCenterLon']


    in_dict, out_dict, score_dict = \
        archer3_visir(image, attrib, first_guess, display_filename='imgs/test1_sw.png')


if False:

    data = loadmat('./data/20200208T083000_vis.mat')

    image = {}
    image['bt_grid'] = data['aSetSaved'][0][0]['bv']
    image['lat_grid'] = data['aSetSaved'][0][0]['latMatrix']
    image['lon_grid'] = data['aSetSaved'][0][0]['lonMatrix']
    #image['time_arr'] = Does not exist

    attrib = {}
    attrib['sat'] = 'Himawari-8'
    attrib['sensor'] = 'Imager'
    attrib['scan_type'] = 'Geo'
    attrib['archer_channel_type'] = 'Vis'
    attrib['nadir_lon'] = 140.7 #41.5

    first_guess = {}
    first_guess['source'] = 'fx' # 'fx', 'bt', 'archer', 'NHCanfx'
    first_guess['time'] = convert_matlab_time(data['aSetSaved'][0][0]['matDate']) 
    first_guess['vmax'] = data['aSetSaved'][0][0]['estVmax']
    first_guess['lat'] = data['aSetSaved'][0][0]['estCenterLat']
    first_guess['lon'] = data['aSetSaved'][0][0]['estCenterLon']


    in_dict, out_dict, score_dict = \
        archer3_visir(image, attrib, first_guess, display_filename='imgs/test1_vis.png')


if True:

    # Apply data grabbed from a CLASS file:

    filename = './data/GATMO-SATMS_j01_d20200208_t0536586_e0544583_b11518_c20200415211227285695_noac_ops.h5'

    with h5py.File(filename, 'r') as f:
         bt_mx = f['All_Data']['ATMS-SDR_All']['BrightnessTemperature'][()]
        lat_mx = f['All_Data']['ATMS-SDR-GEO_All']['Latitude'][()]
        lon_mx = f['All_Data']['ATMS-SDR-GEO_All']['Longitude'][()]
        zen_mx = f['All_Data']['ATMS-SDR-GEO_All']['SatelliteZenithAngle'][()]
        azm_mx = f['All_Data']['ATMS-SDR-GEO_All']['SatelliteAzimuthAngle'][()]

    bt_mx  =  bt_mx / 190. - 50 # Just a guess to get this going

    # Vmax
    op_lon = 116.8
    op_lat = -20.8
    op_vmax = 90


    # Collect for ARCHER

    image = {}
    image['bt_grid'] = bt_mx
    image['lat_grid'] = lat_mx
    image['lon_grid'] = lon_mx
    #image['time_arr'] = ?
    image['zen_grid'] = zen_mx
    image['azm_grid'] = azm_mx

    attrib = {}
    attrib['sat'] = 'JPSS-1'
    attrib['sensor'] = 'ATMS'
    attrib['scan_type'] = 'Crosstrack'
    attrib['archer_channel_type'] = '89GHz'

    first_guess = {}
    first_guess['source'] = 'fx' # 'fx', 'bt', 'archer', 'NHCanfx'
    first_guess['time'] = np.NaN # I won't bother with this yet...
    first_guess['vmax'] = op_vmax
    first_guess['lat'] = op_lat
    first_guess['lon'] = op_lon


    in_dict, out_dict, score_dict = \
        archer3_mw(image, attrib, first_guess, display_filename='imgs/demo_atms_89.png')
    """