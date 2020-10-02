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

import h5py
filename = "GIMGO-SVI01-SVI04-SVI05_j01_d20200208_t0540026_e0545426_b11518_c20200410201727083640_noac_ops.h5"

with h5py.File(filename, "r") as f:
    bt_mx = f['All_Data']['VIIRS-I5-SDR_All']['BrightnessTemperature'][()]
    lat_mx = f['All_Data']['VIIRS-IMG-GEO_All']['Latitude'][()]
    lon_mx = f['All_Data']['VIIRS-IMG-GEO_All']['Longitude'][()]

lat_mx = lat_mx[:-33, :-33]
lon_mx = lon_mx[:-33, :-33]
bt_mx = bt_mx[:-33, :-33] / 256. + 140 # Just a guess

op_lon = 116.8
op_lat = -20.8


"""
fig = plt.figure()

# This explains the nan-related runtime error for the following line: https://github.com/matplotlib/matplotlib/issues/9892
img0 = plt.pcolormesh(lon_mx, lat_mx, bt_mx, vmin=180, vmax=300, cmap='jet_r') 
#img0 = plt.scatter(lon_mx, lat_mx, c=bt_mx, vmin=160, vmax=280, cmap='jet_r', edgecolors='none', marker=',')

aspect_lon = np.cos(np.pi/180*op_lat) # Aspect ratio
plt.xlim(op_lon-2.5/aspect_lon, op_lon+2.5/aspect_lon)
plt.ylim(op_lat-2.5, op_lat+2.5)

fig.set_size_inches(5, 5)
fig.colorbar(img0, orientation='horizontal')
plt.savefig('test.png', dpi=100)
"""

fig = plt.figure()

# This explains the nan-related runtime error for the following line: https://github.com/matplotlib/matplotlib/issues/9892
img0 = plt.imshow(bt_mx[500:2200, 3000:4350], vmin=180, vmax=300, cmap='jet_r') 
#img0 = plt.scatter(lon_mx, lat_mx, c=bt_mx, vmin=160, vmax=280, cmap='jet_r', edgecolors='none', marker=',')

fig.set_size_inches(5, 5)
fig.colorbar(img0, orientation='horizontal')
plt.savefig('test.png', dpi=100)


#
#
#

from satpy import available_readers, Scene
from glob import glob
available_readers()

filenames = glob('./GIMGO-SVI01*')
#scn = Scene(reader='viirs_l1b', filenames=filenames) # 'ValueError: No supported files found'
scn = Scene(reader='viirs_sdr', filenames=filenames)
scn.load(['I05'])
plt.figure()
plt.imshow(scn['I05'][500:2200, 3000:4500], vmin=180, vmax=300, cmap='jet_r')
plt.colorbar()
plt.savefig('test_pysat.png', dpi=100)

#scn.save_datasets(writer='cf', datasets=['I05'], filename='viirs_sdr_i05_damien_test.nc', 
#    exclude_attrs=['raw_metadata'])
scn.save_datasets(writer='cf', datasets=['I05'], filename='viirs_sdr_i05_damien_test.nc', 
    exclude_attrs=['raw_metadata'])

