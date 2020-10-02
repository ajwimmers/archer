from __future__ import print_function
import numpy as np
import time, pandas, os
from pyproj import Geod ### REQUIRES SPECIAL INSTALL: conda install pyproj


def antemeridian_decross(lonGrid, first_guess_lon):

    # If the TC is close to the antemeridian, choose a side (positive lon)
    if np.any(lonGrid > 170) and np.any(lonGrid < -170):
        
        print('Revamping lon grid for an antemeridian crossing')
        if first_guess_lon < 0:
            lonGrid[lonGrid > 100] = lonGrid[lonGrid > 100] - 360

        else:
            lonGrid[lonGrid < -100] = lonGrid[lonGrid < -100] + 360

    return lonGrid


def antemeridian_restore(lon):

    # Insure that the longitude value is between -180 and 180
    lon = np.mod(lon + 180, 360) - 180

    return lon


def reduce_step(lonGrid, latGrid, step_km=4):

    # Calculate the row and col step necessary to bring the image to step_km resolution
    # For IR, step_km is 4 (km)

    rpd = np.pi/180;
    km_per_gcd = 6370*rpd

    num_rows, num_cols = np.shape(lonGrid)

    mid_row = int(np.round(num_rows / 2))
    mid_col = int(np.round(num_cols / 2))

    # Use the middle of the grid as the reference
    row_dist_km = km_per_gcd * np.sqrt( 
        ((lonGrid[mid_row, mid_col] - lonGrid[mid_row+1, mid_col]) / np.cos(latGrid[mid_row, mid_col]))**2 +
        (latGrid[mid_row, mid_col] - latGrid[mid_row+1, mid_col])**2 )

    col_dist_km = km_per_gcd * np.sqrt( 
        ((lonGrid[mid_row, mid_col] - lonGrid[mid_row, mid_col+1]) / np.cos(latGrid[mid_row, mid_col]))**2 +
        (latGrid[mid_row, mid_col] - latGrid[mid_row, mid_col+1])**2 )

    # Use "floor" to keep it conservative
    new_row_step = int(max(1, np.floor(step_km / row_dist_km)))
    new_col_step = int(max(1, np.floor(step_km / col_dist_km)))

    return new_row_step, new_col_step


def reduce_res(image, step_km=4):

    new_row_step, new_col_step = reduce_step(image['lon_grid'], image['lat_grid'], step_km=step_km)

    if new_row_step > 1 or new_col_step > 1:

        print('Reducing resolution by factors of', new_row_step, ',', new_col_step)

        for key in image.keys():
            image[key] = image[key][::new_row_step, ::new_col_step]

    return image


def extrap_row1(xMx, yMx):

    xNm2=xMx[:,2]
    xNm1=xMx[:,1]
    xN=xMx[:,0]
    yNm2=yMx[:,2]
    yNm1=yMx[:,1]
    yN=yMx[:,0]

    xNp=2*xNm1-xNm2
    yNp=2*yNm1-yNm2

    xNp1p=2*xN-xNm1
    yNp1p=2*yN-yNm1

    alpha = np.arctan((yN-yNm1)/(xN-xNm1)) - np.arctan((yNp-yNm1)/(xNp-xNm1))

    xNp1 = xN + (xNp1p-xN)*np.cos(alpha) - (yNp1p-yN)*np.sin(alpha)
    yNp1 = yN + (xNp1p-xN)*np.sin(alpha) + (yNp1p-yN)*np.cos(alpha)

    return xNp1, yNp1


def extrap_rowN(xMx, yMx):

    xNm2=xMx[:,-3]
    xNm1=xMx[:,-2]
    xN=xMx[:,-1]
    yNm2=yMx[:,-3]
    yNm1=yMx[:,-2]
    yN=yMx[:,-1]

    xNp=2*xNm1-xNm2
    yNp=2*yNm1-yNm2

    xNp1p=2*xN-xNm1
    yNp1p=2*yN-yNm1

    alpha = np.arctan((yN-yNm1)/(xN-xNm1)) - np.arctan((yNp-yNm1)/(xNp-xNm1))

    xNp1 = xN + (xNp1p-xN)*np.cos(alpha) - (yNp1p-yN)*np.sin(alpha)
    yNp1 = yN + (xNp1p-xN)*np.sin(alpha) + (yNp1p-yN)*np.cos(alpha)

    return xNp1, yNp1


def parallax_fix_conical(lonGrid, latGrid, sensor, structure_height_km):

    # Constants
    rpd = np.pi/180;
    km_per_gcd = 6370*rpd

    if sensor == 'SSMI':
        view_angle_deg = 53.1
        flip_elements = False

    elif sensor == 'SSMIS':
        view_angle_deg = 53.1
        flip_elements = True

    elif sensor == 'TMI':
        view_angle_deg = 53.1
        flip_elements = False

    elif sensor == 'AMSRE' or sensor == 'AMSR2':
        view_angle_deg = 55.2
        flip_elements = True

    elif sensor == 'GMI':
        view_angle_deg = 52.8
        flip_elements = True


    # Calculate nudges and nudge the image to make a first-approximation correction for
    # parallax *according to the most important features in the image*

    if flip_elements:
        lonGrid = lonGrid[:,::-1]
        latGrid = latGrid[:,::-1]

    lonCol0, latCol0 = extrap_row1(lonGrid, latGrid)
    lonGrid = np.asarray(np.append(np.matrix(lonCol0).T, lonGrid, axis=1)) # np.matrix allows transpose
    latGrid = np.asarray(np.append(np.matrix(latCol0).T, latGrid, axis=1))

    lonColN1, latColN1 = extrap_rowN(lonGrid, latGrid)
    lonGrid = np.asarray(np.append(lonGrid, np.matrix(lonColN1).T, axis=1))
    latGrid = np.asarray(np.append(latGrid, np.matrix(latColN1).T, axis=1))

    # Calculate parallax offset based on scanline orientations and view angle
    dist2p = structure_height_km * np.tan(rpd * view_angle_deg) / km_per_gcd
    cos13 = np.cos(rpd * 0.5 * (latGrid[:, :-2] + latGrid[:, 2:]))
    delX13 = cos13 * (lonGrid[:, :-2] - lonGrid[:, 2:])
    delY13 = latGrid[:, :-2] - latGrid[:, 2:]
    dist13 = np.sqrt( delX13**2 + delY13**2 )

    # Apply offset
    lonGrid = lonGrid[:, 1:-1] - dist2p * delY13 / dist13 / cos13
    latGrid = latGrid[:, 1:-1] + dist2p * delX13 / dist13

    if flip_elements:
        lonGrid = lonGrid[:,::-1]
        latGrid = latGrid[:,::-1]
       
    return lonGrid, latGrid


def parallax_fix_geo(lonGrid, latGrid, nadir_lon, sensor, structure_height_km):

    # Constants
    rpd = np.pi/180;
    km_per_gcd = 6370*rpd

    xre = 6370 * np.cos(np.pi/180*latGrid) * np.cos(np.pi/180*lonGrid)
    yre = 6370 * np.cos(np.pi/180*latGrid) * np.sin(np.pi/180*lonGrid)
    zre = 6370 * np.sin(np.pi/180*latGrid)
    xrd = (35788+6370) * np.cos(nadir_lon*np.pi/180) - xre
    yrd = (35788+6370) * np.sin(nadir_lon*np.pi/180) - yre
    zrd = -zre
    remag = np.sqrt((xre**2) + (yre**2) + (zre**2))
    rdmag = np.sqrt((xrd**2) + (yrd**2) + (zrd**2))
    cosz = (xre*xrd + yre*yrd + zre*zrd) / (remag*rdmag)
    delX_m = 1000 * structure_height_km * np.tan(np.arccos(cosz))

    # Calculate the direction of nudging to correct for parallax
    g = Geod(ellps='sphere')
    _, azimuth21, _ = g.inv(
        nadir_lon * np.ones(np.shape(lonGrid)), 0 * np.ones(np.shape(lonGrid)), 
        lonGrid, latGrid)

    # Apply the nudge
    new_lonGrid, new_latGrid, _ = g.fwd(lonGrid, latGrid, azimuth21, delX_m)

    #return lonGrid, latGrid # Control
    return new_lonGrid, new_latGrid


def parallax_fix_crosstrack(
    lonGrid, latGrid, sensor, archer_channel_type, structure_height_km):

    # Constants
    rpd = np.pi/180;
    km_per_gcd = 6370*rpd

    this_dir = os.path.dirname(os.path.realpath(__file__))
    etc_dir = os.path.join(this_dir, '../etc/')

    if sensor == 'ATMS':
        scan_angle_arr = np.abs(np.linspace(-52.77, 52.77, 96))
        nadir_lon_arr = np.mean(lonGrid[:, 47:48], axis=1) # Correct row/col?
        nadir_lat_arr = np.mean(latGrid[:, 47:48], axis=1)
        scan_angle_grid, nadir_lon_grid = np.meshgrid(scan_angle_arr, nadir_lon_arr)
        _              , nadir_lat_grid = np.meshgrid(scan_angle_arr, nadir_lat_arr)

    elif sensor == 'AMSU-B' or sensor == 'MHS':
        scan_angle_arr = pandas.read_csv(
            os.path.join(etc_dir, 'amsub90scanangles.csv'), header=None)
        nadir_lon_arr = np.mean(lonGrid[:, 44:45], axis=1)
        nadir_lat_arr = np.mean(latGrid[:, 44:45], axis=1)
        scan_angle_grid, nadir_lon_grid = np.meshgrid(scan_angle_arr, nadir_lon_arr)
        _              , nadir_lat_grid = np.meshgrid(scan_angle_arr, nadir_lat_arr)

    elif sensor == 'VIIRS' and archer_channel_type == 'DNB':
        scan_angle_arr = pandas.read_csv(
            os.path.join(etc_dir, 'viisr_4064scanangles.csv'), header=None)
        nadir_lon_arr = np.mean(lonGrid[:, 2031:2032], axis=1)
        nadir_lat_arr = np.mean(latGrid[:, 2031:2032], axis=1)
        scan_angle_grid, nadir_lon_grid = np.meshgrid(scan_angle_arr, nadir_lon_arr)
        _              , nadir_lat_grid = np.meshgrid(scan_angle_arr, nadir_lat_arr)

    elif sensor == 'VIIRS':
        scan_angle_arr = pandas.read_csv(
            os.path.join(etc_dir, 'viisr_6400scanangles.csv'), header=None)
        nadir_lon_arr = np.mean(lonGrid[:, 3199:3200], axis=1)
        nadir_lat_arr = np.mean(latGrid[:, 3199:3200], axis=1)
        scan_angle_grid, nadir_lon_grid = np.meshgrid(scan_angle_arr, nadir_lon_arr)
        _              , nadir_lat_grid = np.meshgrid(scan_angle_arr, nadir_lat_arr)

    # Calculate the distance of nudging to correct for parallax
    delX_m = 1000 * structure_height_km * np.tan(rpd * scan_angle_grid)

    # Calculate the direction of nudging to correct for parallax
    g = Geod(ellps='sphere')
    _, azimuth21, _ = g.inv(nadir_lon_grid, nadir_lat_grid, lonGrid, latGrid)

    # Apply the nudge
    new_lonGrid, new_latGrid, _ = g.fwd(lonGrid, latGrid, azimuth21, delX_m)

    #return lonGrid, latGrid # Control
    return new_lonGrid, new_latGrid


def parallax_fix_allnav(lonGrid, latGrid, zenGrid, azmGrid, structure_height_km):

    # Constants
    rpd = np.pi/180;
    km_per_gcd = 6370*rpd

    # Calculate the distance of nudging to correct for parallax
    delX_m = 1000 * structure_height_km * np.tan(rpd * zenGrid)

    # Apply the nudge in the opposite of the azimuthal direction
    g = Geod(ellps='sphere')
    new_lonGrid, new_latGrid, _ = g.fwd(lonGrid, latGrid, azmGrid, -delX_m)

    #return lonGrid, latGrid # Control
    return new_lonGrid, new_latGrid


def cos_solar_zenith(lonGrid, latGrid, time_epoch_secs):

    # Calculates cosine of solar zenith angle over a grid of points
    #
    # Reference: http://en.wikipedia.org/wiki/Insolation
    # Note: The declination calculation on the wiki page looks like it's in
    # error. It uses the angle of orbit, which looks out of place in the
    # calculations. I used my own calculation for "del", and it checks out on
    # the equinox/solstice days.
    #
    # AJW (2011)

    rpd = np.pi/180

    # Relevant times
    img_time = time.gmtime(time_epoch_secs)
    hr_float = img_time.tm_hour + img_time.tm_min/60. + img_time.tm_sec/3600.

    # Terms from the wiki Insolation page:

    # Obliquity of the earth (degrees)
    eta = 23.4398

    # Declination of the earth (degrees)
    delta = eta * np.sin( (284. + img_time.tm_yday) / 365. * 360. * rpd)

    # Angle relative to peak insolation longitude (hours*deg/hr)
    h = lonGrid + (hr_float - 12) * 15 

    cos_solar_zenith_grid = np.sin(latGrid * rpd) * np.sin(delta * rpd) + \
        np.cos(latGrid * rpd) * np.cos(delta * rpd) * np.cos(h * rpd)

    return cos_solar_zenith_grid




