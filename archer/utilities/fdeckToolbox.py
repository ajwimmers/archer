from __future__ import print_function
import numpy as np
import datetime


def generate_string(attrib, in_dict, out_dict, sector_info=None):

    if sector_info is not None:
        bn = sector_info['storm_basin']
        cy = sector_info['storm_num']
    else:
        # Placeholders
        bn = 'BB' 
        cy = 'NN' 

    image_dt = datetime.datetime.utcfromtimestamp(in_dict['time'])
    YYYYMMDDHHMM = image_dt.strftime('%y%m%d%H%M')
    
    fix_format = ' 70'

    # Convert strings for fix_type
    if out_dict['archer_channel_type'] == 'Vis':
        fix_type = 'VISI'
    elif out_dict['archer_channel_type'] == 'IR':
        fix_type = 'INFR'
    elif out_dict['archer_channel_type'] == 'SWIR':
        fix_type = 'SWIR'
    elif attrib['sensor'] == 'SSMI':
        fix_type = 'SSMI'
    elif attrib['sensor'] == 'SSMIS':
        fix_type = 'SSMS'
    elif attrib['sensor'] == 'GMI':
        fix_type = 'GPMI'
    elif attrib['sensor'] == 'AMSR2':
        fix_type = 'AMS2'
    else:
        fix_type = 'UNKN'

    center_int = '         C'

    # Use either the ARCHER center, or None if there is none
    if out_dict['center_lat'] is None:

        flag_indic = 'C'
        lat_str = '     '
        lon_str = '      '

    else:

        flag_indic = ' '
        if out_dict['center_lat'] > 0:
            nshem = 'N'
        else:
            nshem = 'S'

        if out_dict['center_lon'] > 0:
            ewhem = 'E'
        else:
            ewhem = 'W'

        lat_str = '{:4d}'.format(int(np.round(np.abs(out_dict['center_lat']*100)))) + nshem
        lon_str = '{:5d}'.format(int(np.round(np.abs(out_dict['center_lon']*100)))) + ewhem

    height_ob = '     ' # Keep empty?

    # Position confidence
    if out_dict['radius50percCertDeg'] <= 0.15: # Maybe bring down as far as 0.125?
        posit_conf = '1'
    elif out_dict['radius50percCertDeg'] <= 0.30: # Pretty good threshold
        posit_conf = '2'
    else:
        posit_conf = '3'

    # Observation sources
    spaces_23 = '                       ';
    
    if attrib['sensor'][:3] == 'SSM':
        obs_sources = spaces_23 + 'm'
    elif attrib['sensor'] == 'TMI':
        obs_sources = spaces_23 + 't'
    elif attrib['sensor'] == 'ASC':
        obs_sources = spaces_23 + 'c'
    elif out_dict['archer_channel_type'] == 'Vis':
        obs_sources = spaces_23 + 'v'
    elif out_dict['archer_channel_type'] == 'IR':
        obs_sources = spaces_23 + 'i'
    else:
        obs_sources = spaces_23 + 'x'


    # Final comment, NaN conditionals
    if out_dict['ring_radius_deg'] is None or np.isnan(out_dict['ring_radius_deg']):
        ringRadComm = 9.99
    else:
        ringRadComm = out_dict['ring_radius_deg']

    if out_dict['radius50percCertDeg'] is None or np.isnan(out_dict['radius50percCertDeg']):
        ringRad50Comm = 9.99
    else:
        ringRad50Comm = out_dict['radius50percCertDeg']
    
    if out_dict['radius95percCertDeg'] is None or np.isnan(out_dict['radius95percCertDeg']):
        ringRad95Comm = 9.99
    else:
        ringRad95Comm = out_dict['radius95percCertDeg']
   
    if out_dict['eye_prob'] is None or np.isnan(out_dict['eye_prob']):
        ringEPComm = -99
    else:
        ringEPComm = int(np.round(out_dict['eye_prob']))

    sensorComm = 'ARCH';
                    

    comments = 'irad={:4.2f} | r50={:4.2f} | r95={:4.2f} | ep={:3d} | src={}'.format(
        ringRadComm, ringRad50Comm, ringRad95Comm, ringEPComm, sensorComm)


    # Write it all out
    fdeck_str = '{}, {}, {}, {}, {}, {}, {}, {}, {}, {}, {}, {}, {}, {}, {}'.format(
        bn, cy, YYYYMMDDHHMM, fix_format, fix_type, center_int, flag_indic, lat_str, lon_str, height_ob, posit_conf,  # First 11 vars
        '   ,  ,     ,  ,     ,    ,     ,     ,     ,     ,     ,  ,  ,  ,  ,  ,    ,    ,  ,  CIMS, AUT', # Dummy line
        '   ,             ,             ,    ,     ', obs_sources, comments) # Dummy line plus two more strings


    return fdeck_str
    