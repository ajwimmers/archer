
from __future__ import print_function
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
from scipy.interpolate import interp1d
from scipy.interpolate import RegularGridInterpolator
#from importlib import reload



def amsre2tmiConversion(amsre_bts):
    """
    Convert amsre-like BTs to tmi-like BTs. From art/ver3/utilities.
    AJW, CIMSS, June 2021
    """

    amsre_training_points = np.arange(130, 286, 5)
    tmi_training_points = np.concatenate([np.arange(130,151,5), [155+1], np.arange(160,186,5)+8, [190+6.5], 
        [195+5.7], np.arange(200,221,5)+5, [225+4], [230+3], [235+2], [240+1], [245-1], [250-2], [255-3], np.arange(260,286,5)-4])

    print(len(amsre_training_points))
    print(len(tmi_training_points))

    f = interp1d(amsre_training_points, tmi_training_points, 'cubic')
    tmi_bts = f(amsre_bts)

    return tmi_bts



def CalcInten85(sensor, out_dict, score_dict):
    """
    Adapted from art/ver3/scripts/utilities/calcInten85.m
    This is a mess. I welcome someone cleaning this up.
    Note that the finalScore will be about 1-2K lower, because archer uses linear interpolation, not cubic.
    AJW, CIMSS, June 2021
    """

    # Normalize for 89-GHz sensors
    if sensor in ['AMSR2', 'AMSRE', 'GMI']:
        normalizedDataGrid = amsre2tmiConversion(score_dict['data_grid1'])
    else:
        normalizedDataGrid = score_dict['data_grid1']

    # Abbreviate some names
    lonGrid = score_dict['lon_grid1']
    latGrid = score_dict['lat_grid1']
    lonArr = lonGrid[0,:]

    # This makes it work with RegularGridInterpolator
    latGrid = np.flip(latGrid, axis=0)
    latArr = latGrid[:,0]
    normalizedDataGrid = np.flip(normalizedDataGrid, axis=0)
    f = RegularGridInterpolator((lonArr, latArr), normalizedDataGrid.T, 'nearest') # RGI expects 'ij' indexing


    # Use either weak or normal target point
    if out_dict['center_lon'] is None:
        archerLon = out_dict['weak_center_lon']
        archerLat = out_dict['weak_center_lat']
    else:
        archerLon = out_dict['center_lon']
        archerLat = out_dict['center_lat']

    # Find the ring with the strongest eyewall
    maxCircleBT = 1e12 # Really big number
    bestEyewallBTarr = np.array([])
    for radCircle in np.arange(out_dict['ring_radius_deg'], out_dict['ring_radius_deg']+0.301, 0.05):
        lonCircle = archerLon + np.cos(np.arange(0,360,5)*np.pi/180) * radCircle / np.cos(np.pi/180*archerLat)
        latCircle = archerLat + np.sin(np.arange(0,360,5)*np.pi/180) * radCircle
        circleBTarr = f(list(zip(lonCircle, latCircle)))
        circleBTarrNN = circleBTarr[~np.isnan(circleBTarr)]
        if np.max(circleBTarrNN) < maxCircleBT:
            maxCircleBT = np.max(circleBTarrNN)
            bestEyeRadDeg = radCircle
            bestEyewallBTarr = circleBTarr

    # Find the highest BT inside that ring
    ptsInsideRing = ((lonGrid-archerLon)/np.cos(np.pi/180*archerLat))**2 + (latGrid-archerLat)**2 <= bestEyeRadDeg**2
    maxEyeBT = np.max(np.concatenate([normalizedDataGrid[ptsInsideRing], bestEyewallBTarr]))

    # Calculate the scores as in comboCenterInten
    eyeWallArr = bestEyewallBTarr
    eyeWallArrNN = eyeWallArr[~np.isnan(eyeWallArr)]
    fractionEyeWallNN = len(eyeWallArrNN)/len(eyeWallArr)
    if fractionEyeWallNN < 0.85:
        eyeScore = np.NaN
        diffScore = np.NaN
        finalScore = np.NaN
    else:
        fractionBelowCeiling = np.sum( np.logical_and(eyeWallArrNN<232, maxEyeBT-eyeWallArrNN>10) ) / len(eyeWallArrNN)
        allBelow255 = np.sum(eyeWallArrNN<255) == len(eyeWallArrNN)
        if allBelow255:
            fractionBelowHotSpot = np.sum(maxEyeBT-eyeWallArrNN > 20) / len(eyeWallArrNN)
        else:
            fractionBelowHotSpot = 0

        eyeScore = np.max([fractionBelowCeiling, fractionBelowHotSpot])
        diffScore = maxEyeBT - np.max(eyeWallArrNN)
        finalScore = (eyeScore>0.85) * 15 + diffScore

    # This was added to prevent ridiculous eyewall diameters (10/2009)
    # And it was changed from .80 to .50 after other error checks were added (8/2010)
    if np.isnan(eyeScore) or eyeScore < .50:
        bestEyeRadDeg = np.NaN

    # Put together the scoreStats dict
    scoreStats = {}
    scoreStats['hotSpot'] = maxEyeBT
    scoreStats['weakLink'] = np.max(eyeWallArrNN)
    scoreStats['diffScore'] = diffScore
    scoreStats['completeness'] = eyeScore


    return bestEyeRadDeg, finalScore, scoreStats


def PostProcessCalcInten85(in_dict, out_dict, finalScore, scoreStats):
    """
    Adapted from art/ver3/scripts/utilities/compileIntenScore.m
    AJW, CIMSS, June 2021
    """

    # Mark dicy cases as negative so that it is not used in the ADT
    if in_dict['op_vmax'] < 55:
        finalScore = -finalScore

    # Prevent a false positive in some random structure of convective cells
    ringRadius = out_dict['ring_radius_deg']
    if scoreStats['hotSpot'] < 200:
        finalScore = -920
        ringRadius = np.NaN

    if np.isnan(ringRadius):
        sfcDiamKM = np.NaN
    elif ringRadius > 0.051 or ringRadius < 0:
        sfcDiamKM = ringRadius * 2 * 111.18 - 10 # Conversion to km
    else:
        sfcDiamKM = 5 # When ringRadius==0.05, the calculation must be different

    return finalScore, sfcDiamKM



