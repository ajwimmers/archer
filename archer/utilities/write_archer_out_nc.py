
import time, calendar
import numpy as np
import netCDF4
from os import remove
from os.path import isfile

def write_rdf_paths_netcdf(filename, lonPathGrid, latPathGrid, 
                           sourceTime, destinTime, lonArr, latArr, netcdf_format):

    # This is necessary to fix a bug in the Netcdf4 module with file permission settings
    if isfile(filename):
        remove(filename)

    rootgrp = netCDF4.Dataset(filename, 'w' ,format=netcdf_format)
    numRows, numCols = np.shape(lonPathGrid)
    rootgrp.createDimension('lat', numRows)
    rootgrp.createDimension('lon', numCols)
    rootgrp.createDimension('time', None)

    nc_sourceTime = rootgrp.createVariable('sourceTime', 'i8', ('time',))
    nc_sourceTime[:] = sourceTime
    nc_destinTime = rootgrp.createVariable('destinTime', 'i8', ('time',))
    nc_destinTime[:] = destinTime

    nc_lonArr = rootgrp.createVariable('lonArr', 'f8', ('lon',))
    nc_lonArr[:] = lonArr
    nc_latArr = rootgrp.createVariable('latArr', 'f8', ('lat',))
    nc_latArr[:] = latArr
    nc_lonPathGrid = rootgrp.createVariable('lonPathGrid', 'f8', ('lat','lon',), least_significant_digit=2)
    nc_lonPathGrid[:] = lonPathGrid
    nc_latPathGrid = rootgrp.createVariable('latPathGrid', 'f8', ('lat','lon',), least_significant_digit=2)
    nc_latPathGrid[:] = latPathGrid
    rootgrp.close()
    
