# -*- coding: utf-8 -*-
"""
Created on Wed Feb 24 10:56:54 2016

@author: wimmers
"""

from __future__ import print_function
import numpy as np
from scipy.interpolate import griddata
#from scipy import interpolate


def InterpSwathToGlobalRectGridSet(lonSwath, latSwath, dataSwath, footprintSwath, timeSwath, lonGridArr, latGridArr):

    """
    Special case of nearest-neighbor interpolation for three of the same grids
    """
    
    # Make the "index swath"
    numScans, numElems = np.shape(latSwath)
    indexSwath = np.arange(0, numScans*numElems).reshape(numScans, numElems)
    
    # Interpolate to the "index grid"
    indexInterpGrid = InterpSwathToGlobalRectGrid(lonSwath, latSwath, indexSwath, lonGridArr, latGridArr, 'nearest')

    # These are the grid indeces to work from
    IsInIndexGrid = np.logical_not( np.isnan(indexInterpGrid) )
    
    # These are the swath indeces to pull into the output grids
    indexSwathList = indexInterpGrid[IsInIndexGrid].astype(int)
 
    # Create the output grids in each the same way
    dataInterpGrid = np.nan * np.zeros( np.shape(indexInterpGrid) )
    dataInterpGrid[IsInIndexGrid] = dataSwath.flatten()[indexSwathList]

    footprintInterpGrid = np.nan * np.zeros( np.shape(indexInterpGrid) )
    footprintInterpGrid[IsInIndexGrid] = footprintSwath.flatten()[indexSwathList]

    timeInterpGrid = np.nan * np.zeros( np.shape(indexInterpGrid) )
    timeInterpGrid[IsInIndexGrid] = timeSwath.flatten()[indexSwathList]

    return dataInterpGrid, footprintInterpGrid, timeInterpGrid


def InterpSwathToGlobalRectGrid(lonSwath, latSwath, dataSwath, lonGridArr, latGridArr, interpType):
    """
    Computes an interpoated rectangular grid from a single swath, parsing the swath by zone (north pole,
    midlatitude, south pole) in order to get a fast and accurate interpolation in all zones.  
    
    Inputs:
        lonSwath, latSwath: Navigation of the swath
        dataSwath: Data, such as TPW
        lonGridArr, latGridArr: 1D arrays of the new grid navigation. NOTE! lonGridArr
            must not have redundant lons on the left and right. Use -180 and 179.75, for example.
        interpType: Interpolation type for griddata - 'linear', 'cubic', 'nearest', etc.
            NOTE! Currently 'cubic' does not work b/c it conflicts with the nan mask. Use
            'linear' or 'nearest' instead.
    
    Returns:
        dataInterpGrid: Interpolated grid of data in the grid navigation
        
    Feb 2016 (AJW)
    """
    
    # Useful constants
    numScans, numElems = np.shape(latSwath)
    centerElem = np.floor(numElems/2)
    
    # Initialize output array
    dataInterpGrid = np.NaN * np.zeros((np.size(latGridArr), np.size(lonGridArr)))  
    
    # Divide the swath into overlapping zone sections
    zoneCodeArr = np.zeros((numScans, 1)) # Midlatitude (rectangular)
    zoneCodeArr[latSwath[:,centerElem] >  60] =  1 # North (polar)
    zoneCodeArr[latSwath[:,centerElem] < -60] = -1 # South (polar)

    # Iterate by scan, and operate on a new section when we get to the end of it.
    # Allow the first and last scans to just be part of their neighboring scans.
    sectionInitIdx = 0
    for scanIdx in range(1, numScans-1):
        
        # This finds the end of the ongoing section (North polar, midlatitude or south polar)
        if zoneCodeArr[scanIdx] != zoneCodeArr[scanIdx-1] or scanIdx == numScans-3:
            
            sectionEndIdx = scanIdx + 2
            zoneCode = zoneCodeArr[scanIdx-1]
            #print 'Start Idx: ', sectionInitIdx, ', Final Idx: ', sectionEndIdx
            #print 'Zone Code: ', zoneCode
            
            # Interpolate this section in the proper way to the rectangular grid
            dataInterpSection = InterpSectionToGlobalRectGrid(
                lonSwath[sectionInitIdx:sectionEndIdx, :],
                latSwath[sectionInitIdx:sectionEndIdx, :],
                dataSwath[sectionInitIdx:sectionEndIdx, :],
                lonGridArr, latGridArr, interpType, zoneCode)           
            
            # Add the new interpolated section to the grid
            dataInterpGrid[np.logical_not(np.isnan(dataInterpSection))] = \
                dataInterpSection[np.logical_not(np.isnan(dataInterpSection))]
            
            # Set the new start index for the next section
            sectionInitIdx = scanIdx
            
    return dataInterpGrid
    
            
def InterpSectionToGlobalRectGrid(lonSection, latSection, dataSection, lonGridArr, latGridArr, interp_type='cubic', zoneCode=0):
    """
    Computes an interpoated rectangular grid from a single swath section, according to its zone (north pole,
    midlatitude, south pole).
    
    Inputs:
        lonSection, latSection: Navigation of the swath section
        dataSection: Data, such as TPW
        lonGridArr, latGridArr: 1D arrays of the new grid navigation
        interpType: Interpolation type for griddata - 'linear', 'cubic', 'nearest', etc.
        zoneCode: 1 (North pole zone), 0 (Midlatiude zone), or -1 (South pole zone)
    
    Returns:
        dataInterpSection: Interpolated grid of section's data in the grid navigation
        
    Feb 2016 (AJW)
    """

    # Useful constants
    RPD = np.deg2rad(1) # Radians per degree
    numScans, numElems = np.shape(latSection)
    centerElem = int(np.floor(numElems/2))
    centerScan = int(np.floor(numScans/2))
    lonGridRes = np.abs(lonGridArr[1] - lonGridArr[0])
    BUFFER_XY = 1.5 # Interpolate around the borders of the swath by this much 
    
    # Initialize grid and section variables
    lonGrid, latGrid = np.meshgrid(lonGridArr, latGridArr)
    dataInterpSection = np.nan * lonGrid
    xSectionBuff = np.nan * np.zeros((numScans + 4, numElems + 4))
    ySectionBuff = np.nan * np.zeros((numScans + 4, numElems + 4))
    dataSectionBuff = np.nan * np.zeros((numScans + 4, numElems + 4))
    
    # Convert according to zone:
    
    if zoneCode == 0:
        
        # Convert to a swath centered around zero
        centerLon = np.round(lonSection[centerScan, centerElem]/lonGridRes) * lonGridRes
        xSection = np.mod(lonSection - centerLon + 180, 360) - 180
        ySection = latSection
        
        # Work out the reordering for grids as well
        xGridOrigArr = np.mod(lonGridArr - centerLon + 180, 360) - 180
        centeredLonOrder = np.argsort(xGridOrigArr)
        centeredLonReorder = np.argsort(centeredLonOrder)
        #print centeredLonOrder
        #print centeredLonReorder

        # Keep *grid* in original coordinates, but then rearrange with 
        # centeredLonReorder at the end
        xGrid = lonGrid
        yGrid = latGrid
        
    if zoneCode == 1:
        
        # Convert from polar to rectangular coordinates
        xSection = (90 - latSection) * np.cos(RPD*lonSection)
        ySection = (90 - latSection) * np.sin(RPD*lonSection)
        xGrid = (90 - latGrid) * np.cos(RPD*lonGrid)
        yGrid = (90 - latGrid) * np.sin(RPD*lonGrid)
        
    if zoneCode == -1:
        
        # Convert from polar to rectangular coordinates
        xSection = (90 + latSection) * np.cos(RPD*lonSection)
        ySection = (90 + latSection) * np.sin(RPD*lonSection)
        xGrid = (90 + latGrid) * np.cos(RPD*lonGrid)
        yGrid = (90 + latGrid) * np.sin(RPD*lonGrid)
    
    # Use x, y, data sections with buffers of data and nans at the perimeter to 
    # set it up for clean interpolation at edges
    xSectionBuff, ySectionBuff, dataSectionBuff = addEdgeBuffer(xSection, ySection, dataSection, interp_type)        
    
    # Do the interpolation with griddata on the subsection defined by inSubGrid
    gridCoordinates = list(zip(xSectionBuff.ravel(), ySectionBuff.ravel()))
    
    inSubGrid = np.logical_and( 
        np.logical_and(xGrid > np.min(xSectionBuff.flatten() - BUFFER_XY), 
                       xGrid < np.max(xSectionBuff.flatten() + BUFFER_XY) ), 
        np.logical_and(yGrid > np.min(ySectionBuff.flatten() - BUFFER_XY), 
                       yGrid < np.max(ySectionBuff.flatten() + BUFFER_XY) ) )
    

    try:
        dataInterpSection[inSubGrid] = griddata(gridCoordinates, dataSectionBuff.ravel(), 
            (xGrid[inSubGrid], yGrid[inSubGrid]), method=interp_type) 
    except:
        print('Error: Could not interpolate with griddata on this section')
        dataInterpSection[inSubGrid] = np.nan
    
    # Shift back from recentering if remapped in zoneCode 0 (midlatitudes)
    if zoneCode == 0:
        dataInterpSection = dataInterpSection[:, centeredLonReorder]
    
    return dataInterpSection
    
    
def addEdgeBuffer(xSection, ySection, dataSection, interpType):
    """
    Put a buffer of data at the outer edge of the swath to allow the data to be interpolated
    as far as can be justified, and then put a buffer of nans around that to prevent interpolation
    artifacts.
    
    Inputs:
        xSection, ySection: Navigation of swath section
        dataSection: Corresponding data
        interpType: 'linear', 'nearest', 'cubic'
        
    Returns:
        xSectionBuff, ySectionBuff, dataSectionBuff: Same as inputs, but with 
            buffered edges
            
    Feb 2016 (AJW)
    """
 
    # Set up the buffered arrays
    numScans, numElems = np.shape(xSection)
    xSectionBuff = np.nan * np.zeros((numScans + 4, numElems + 4))
    ySectionBuff = np.nan * np.zeros((numScans + 4, numElems + 4))
    dataSectionBuff = np.nan * np.zeros((numScans + 4, numElems + 4))

    # Add a buffer of perimeter values to allow interpolation to true size of obs

    # Spacing of 0.5 puts a value on the edge of the swath
    if interpType == 'linear' or interpType == 'cubic':
        xSectionBuff[2:-2, 2:-2] = xSection
        ySectionBuff[2:-2, 2:-2] = ySection
        dataSectionBuff[2:-2, 2:-2] = dataSection
        xSectionBuff[1,:] = 1.5*xSectionBuff[2,:] - 0.5*xSectionBuff[3,:]
        ySectionBuff[1,:] = 1.5*ySectionBuff[2,:] - 0.5*ySectionBuff[3,:]
        dataSectionBuff[1,:] = dataSectionBuff[2,:]
        xSectionBuff[-2,:] = 1.5*xSectionBuff[-3,:] - 0.5*xSectionBuff[-4,:]
        ySectionBuff[-2,:] = 1.5*ySectionBuff[-3,:] - 0.5*ySectionBuff[-4,:]
        dataSectionBuff[-2,:] = dataSectionBuff[-3,:]
        xSectionBuff[:,1] = 1.5*xSectionBuff[:,2] - 0.5*xSectionBuff[:,3]
        ySectionBuff[:,1] = 1.5*ySectionBuff[:,2] - 0.5*ySectionBuff[:,3]
        dataSectionBuff[:,1] = dataSectionBuff[:,2]
        xSectionBuff[:,-2] = 1.5*xSectionBuff[:,-3] - 0.5*xSectionBuff[:,-4]
        ySectionBuff[:,-2] = 1.5*ySectionBuff[:,-3] - 0.5*ySectionBuff[:,-4]
        dataSectionBuff[:,-2] = dataSectionBuff[:,-3]

    # Even spacing allows for a fair distance between data and surroundings
    if interpType == 'nearest':
        xSectionBuff[2:-2, 2:-2] = xSection
        ySectionBuff[2:-2, 2:-2] = ySection
        dataSectionBuff[2:-2, 2:-2] = dataSection
        xSectionBuff[1,:] = 2*xSectionBuff[2,:] - xSectionBuff[3,:]
        ySectionBuff[1,:] = 2*ySectionBuff[2,:] - ySectionBuff[3,:]
        dataSectionBuff[1,:] = np.nan
        xSectionBuff[-2,:] = 2*xSectionBuff[-3,:] - xSectionBuff[-4,:]
        ySectionBuff[-2,:] = 2*ySectionBuff[-3,:] - ySectionBuff[-4,:]
        dataSectionBuff[-2,:] = np.nan
        xSectionBuff[:,1] = 2*xSectionBuff[:,2] - xSectionBuff[:,3]
        ySectionBuff[:,1] = 2*ySectionBuff[:,2] - ySectionBuff[:,3]
        dataSectionBuff[:,1] = np.nan
        xSectionBuff[:,-2] = 2*xSectionBuff[:,-3] - xSectionBuff[:,-4]
        ySectionBuff[:,-2] = 2*ySectionBuff[:,-3] - ySectionBuff[:,-4]
        dataSectionBuff[:,-2] = np.nan
    
    # Add an envelope of nans to prevent interpolation artifacts
    xSectionBuff[0,:] = 2.0*xSectionBuff[1,:] - 1.0*xSectionBuff[2,:]
    ySectionBuff[0,:] = 2.0*ySectionBuff[1,:] - 1.0*ySectionBuff[2,:]
    xSectionBuff[-1,:] = 2.0*xSectionBuff[-2,:] - 1.0*xSectionBuff[-3,:]
    ySectionBuff[-1,:] = 2.0*ySectionBuff[-2,:] - 1.0*ySectionBuff[-3,:]
    xSectionBuff[:,0] = 2.0*xSectionBuff[:,1] - 1.0*xSectionBuff[:,2]
    ySectionBuff[:,0] = 2.0*ySectionBuff[:,1] - 1.0*ySectionBuff[:,2]
    xSectionBuff[:,-1] = 2.0*xSectionBuff[:,-2] - 1.0*xSectionBuff[:,-3]
    ySectionBuff[:,-1] = 2.0*ySectionBuff[:,-2] - 1.0*ySectionBuff[:,-3]
                        
    return xSectionBuff, ySectionBuff, dataSectionBuff
    