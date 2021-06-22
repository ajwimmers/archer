"""
This script demos and tests out the ad hoc algorithm for 85-92 GHz imagery intensity scoring
"""


from scipy.io import loadmat
import time, calendar, datetime
import numpy as np
import netCDF4, h5py # Note that you might get funny behavior if you import these in reverse order
from archer.archer4  import archer4

import archer.utilities.CalcInten85 as calc85


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
        archer4(image, attrib, first_guess, para_fix=True) #, display_filename='imgs/demo_conic_partial_89h.png')

    # Test intensity scoring
    bestEyeRadDeg, finalScore, scoreStats = calc85.CalcInten85(attrib['sensor'], out_dict, score_dict)
    finalScore, sfcDiamKM = calc85.PostProcessCalcInten85(in_dict, out_dict, finalScore, scoreStats)



if True or calc_all:

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
        archer4(image, attrib, first_guess, para_fix=True) #, display_filename='imgs/demo_conic_89h.png')

    # Test intensity scoring
    bestEyeRadDeg, finalScore, scoreStats = calc85.CalcInten85(attrib['sensor'], out_dict, score_dict)
    finalScore, sfcDiamKM = calc85.PostProcessCalcInten85(in_dict, out_dict, finalScore, scoreStats)


